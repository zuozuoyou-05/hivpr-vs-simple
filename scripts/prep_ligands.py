"""准备配体库：挑分子，生成 3D 构象，转成对接用的 pdbqt。

一次跑完三件事
  1. 解析 DUD-E hivpr 的 .ism 文件，分出阳性、实测阴性、诱饵三组
  2. 尺寸过滤、去掉重复结构、按配置抽样，写出 data/processed/dataset.csv
  3. SMILES → 3D 构象 → meeko → AutoDock pdbqt

三组的分工
  阳性  30 个实测有活性的化合物，用来检验「打分能不能排序」以及当正样本算富集
  阴性  40 个实测无活性的化合物，用来检验假阳性
  诱饵  130 个性质匹配的化合物，模拟真实筛选库

产物
  data/processed/dataset.csv           配体清单，带 p_activity 标签
  data/processed/reference_ligand.sdf  晶体配体 189，定义盒子、也用来做重对接
  data/processed/sdf/<id>.sdf          单分子 3D
  data/processed/ligands_3d.sdf        多分子 3D
  data/processed/pdbqt/<id>.pdbqt      对接输入
"""

import math
import random
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import cfg
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, SDWriter, rdMolDescriptors
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs, ETKDGv3
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMoleculeConfs
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

N_CONFS = 10  # 每个分子生成几个构象取最低能量。跟 Interformer 官方预处理脚本一致


def parse_ism(path):
    """DUD-E 的 .ism 是空格分隔的定长字段。

    列：SMILES 编号 chembl_id 活性类型 关系 数值 单位 ... uniprot 靶点名
    诱饵文件只有前两列有内容。
    """
    rows = []
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.split()
        if len(f) < 2:
            continue
        rec = {"smiles": f[0], "src_id": f[1]}
        if len(f) >= 15:
            rec.update(chembl_id=f[2], activity_type=f[3], relation=f[4],
                       value_nM=f[5], uniprot=f[9])
        rows.append(rec)
    return pd.DataFrame(rows)


def annotate(df):
    """补 RDKit 描述符。解析不出来的分子直接丢掉。"""
    out = []
    for rec in df.to_dict("records"):
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is None:
            continue
        rec["canonical_smiles"] = Chem.MolToSmiles(mol)
        rec["heavy_atoms"] = mol.GetNumHeavyAtoms()
        rec["mw"] = round(Descriptors.MolWt(mol), 2)
        rec["num_rot_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        try:
            rec["murcko_scaffold"] = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        except Exception:
            rec["murcko_scaffold"] = ""
        out.append(rec)
    return pd.DataFrame(out)


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def diversify(pool, n, rng):
    """按 Murcko 骨架贪心去冗余，不够了再放宽。"""
    idx = list(pool.index)
    rng.shuffle(idx)
    picked, seen = [], set()
    for i in idx:
        scaffold = pool.at[i, "murcko_scaffold"]
        if scaffold and scaffold in seen:
            continue
        seen.add(scaffold)
        picked.append(i)
        if len(picked) == n:
            return pool.loc[picked]
    for i in idx:
        if i not in picked:
            picked.append(i)
        if len(picked) == n:
            break
    return pool.loc[picked[:n]]


def pick_dataset(cfgd, raw, rng):
    """从四个 .ism 文件里挑出最终的配体清单。"""
    pos = annotate(parse_ism(raw / "actives_nM_chembl.ism"))
    neg = annotate(parse_ism(raw / "inactives_nM_chembl.ism"))
    dec = annotate(parse_ism(raw / "decoys_final.ism"))
    # 诱饵文件只有 SMILES 和编号两列，其余字段补空，后面拼表时才不会缺列
    for df in (pos, neg, dec):
        for col in ("chembl_id", "activity_type", "relation", "value_nM"):
            if col not in df.columns:
                df[col] = ""
    print(f"  DUD-E 原始条目  actives={len(pos)}  inactives={len(neg)}  decoys={len(dec)}")

    # 阳性只要实测的等号值，且落在配置的 pActivity 区间里
    pos["value_nM_f"] = pos["value_nM"].map(to_float)
    pos = pos[(pos["relation"] == "=") & pos["activity_type"].isin(["Ki", "IC50", "Kd"])
              & (pos["value_nM_f"] > 0)].copy()
    pos["p_activity"] = 9.0 - pos["value_nM_f"].map(math.log10)
    pos = pos[pos["p_activity"].between(cfg.f(cfgd, "min_p_activity"),
                                        cfg.f(cfgd, "max_p_activity"))]

    # 阴性要的是「测过、确实没活性」，也就是 > 40 uM 还没抑制
    neg["value_nM_f"] = neg["value_nM"].map(to_float)
    neg = neg[(neg["relation"] == ">") & (neg["value_nM_f"] >= 40000)].copy()
    neg["p_activity"] = 9.0 - neg["value_nM_f"].map(math.log10)

    # 尺寸过滤 + 跨组去重，优先级 阳性 > 阴性 > 诱饵
    def keep(d):
        return d[(d["heavy_atoms"] <= cfg.i(cfgd, "max_heavy_atoms"))
                 & (d["mw"] <= cfg.f(cfgd, "max_mw"))]

    pos, neg, dec = keep(pos), keep(neg), keep(dec)
    used, uniq = set(), {}
    for name, df in (("阳性", pos), ("阴性", neg), ("诱饵", dec)):
        df = df[~df["canonical_smiles"].isin(used)].drop_duplicates("canonical_smiles")
        used |= set(df["canonical_smiles"])
        uniq[name] = df
    pos, neg, dec = uniq["阳性"], uniq["阴性"], uniq["诱饵"]
    print(f"  过滤去重后      阳性={len(pos)}  阴性={len(neg)}  诱饵={len(dec)}")

    # 阳性按 pActivity 分 6 档，每档抽几个，保证活性跨度而不是全抽到最强的
    bins = pd.cut(pos["p_activity"], bins=6)
    per_bin = max(1, cfg.i(cfgd, "n_positives") // 6)
    chunks = [diversify(grp, min(per_bin, len(grp)), rng)
              for _, grp in pos.groupby(bins, observed=True)]
    pos_sel = pd.concat(chunks)
    if len(pos_sel) < cfg.i(cfgd, "n_positives"):
        rest = pos.drop(index=pos_sel.index)
        pos_sel = pd.concat([pos_sel, diversify(rest, cfg.i(cfgd, "n_positives") - len(pos_sel), rng)])
    pos_sel = pos_sel.head(cfg.i(cfgd, "n_positives")).sort_values("p_activity", ascending=False)

    neg_sel = diversify(neg, min(cfg.i(cfgd, "n_negatives"), len(neg)), rng)
    dec_sel = dec.sample(n=min(cfg.i(cfgd, "n_decoys"), len(dec)),
                         random_state=cfg.i(cfgd, "seed"))

    def pack(df, group, source, prefix):
        d = df.copy().reset_index(drop=True)
        d["group"] = group
        d["source"] = source
        d["ligand_id"] = [d.at[i, "chembl_id"] if str(d.at[i, "chembl_id"]).startswith("CHEMBL")
                          else f"{prefix}{i + 1:04d}" for i in d.index]
        return d

    final = pd.concat([
        pack(pos_sel, "positive", "DUD-E actives_nM_chembl", "POS"),
        pack(neg_sel, "negative", "DUD-E inactives_nM_chembl", "NEG"),
        pack(dec_sel, "decoy", "DUD-E decoys_final", "DEC"),
    ]).reset_index(drop=True)
    cols = ["ligand_id", "smiles", "canonical_smiles", "chembl_id", "group", "source",
            "activity_type", "relation", "value_nM", "p_activity",
            "heavy_atoms", "mw", "num_rot_bonds", "murcko_scaffold"]
    return final[cols]


def smiles_to_3d(task):
    """SMILES → 能量最低的 3D 构象。ETKDGv3(seed=42) 多构象 + UFF 优化。"""
    lid, smiles = task
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return lid, None
    mol = Chem.AddHs(mol)
    ps = ETKDGv3()
    ps.randomSeed = 42
    ps.useRandomCoords = True
    ps.numThreads = 1
    try:
        if len(EmbedMultipleConfs(mol, N_CONFS, ps)) == 0:
            return lid, None
        energies = [e for _, e in UFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=2000)]
    except Exception:
        return lid, None
    out = Chem.Mol(mol, False, int(min(range(len(energies)), key=lambda i: energies[i])))
    out.SetProp("_Name", lid)
    return lid, out


def sdf_to_pdbqt(task):
    """单分子 SDF → AutoDock pdbqt。"""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    sdf_path, out_dir = task
    mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=True)
    if mol is None:
        return sdf_path.stem, False
    try:
        setups = MoleculePreparation().prepare(mol)
    except Exception:
        return sdf_path.stem, False
    if not setups:
        return sdf_path.stem, False
    text, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return sdf_path.stem, False
    (Path(out_dir) / f"{sdf_path.stem}.pdbqt").write_text(text, encoding="utf-8")
    return sdf_path.stem, True


def main():
    cfgd = cfg.load()
    rng = random.Random(cfg.i(cfgd, "seed"))
    tag = cfgd["pdb_tag"]
    proc = cfg.out("data", "processed")
    sdf_dir = cfg.out("data", "processed", "sdf")
    pdbqt_dir = cfg.out("data", "processed", "pdbqt")
    workers = cfg.i(cfgd, "workers")

    print("[prep_ligands] 1/3 挑配体")
    ds = pick_dataset(cfgd, cfg.path("data", "raw", "dude"), rng)
    ds.to_csv(proc / "dataset.csv", index=False)
    print(f"  写出 dataset.csv  {dict(Counter(ds['group']))}，共 {len(ds)} 个")

    # 晶体配体单独存一份：定义盒子、做重对接、也给 Interformer 定位口袋
    ref = Chem.MolFromMol2File(str(cfg.path(cfgd["crystal_ligand"])), removeHs=False, sanitize=True)
    assert ref is not None, "读不到晶体配体"
    ref.SetProp("_Name", f"{tag}_ref")
    with SDWriter(str(proc / "reference_ligand.sdf")) as w:
        w.write(ref)
    print(f"  参考配体 {tag}_ref  {ref.GetNumHeavyAtoms()} 个重原子")

    print("[prep_ligands] 2/3 生成 3D 构象")
    tasks = [(r.ligand_id, r.smiles) for r in ds.itertuples()]
    with Pool(workers) as pool:
        results = pool.map(smiles_to_3d, tasks)
    n_ok = 0
    failed = []
    with SDWriter(str(proc / "ligands_3d.sdf")) as writer:
        for lid, mol in results:
            if mol is None:
                failed.append(lid)
                continue
            # 分子经多进程来回后 _Name 会丢，这里重新写回去，否则后面按名字取分子会错位
            mol.SetProp("_Name", str(lid))
            with SDWriter(str(sdf_dir / f"{lid}.sdf")) as w:
                w.write(mol)
            writer.write(mol)
            n_ok += 1
    print(f"  成功 {n_ok}/{len(results)}" + (f"  失败 {failed[:8]}" if failed else ""))

    ref.SetProp("_Name", f"{tag}_ref")
    with SDWriter(str(sdf_dir / f"{tag}_ref.sdf")) as w:
        w.write(ref)

    print("[prep_ligands] 3/3 转 pdbqt")
    sdfs = sorted(sdf_dir.glob("*.sdf"))
    with Pool(workers) as pool:
        res = pool.map(sdf_to_pdbqt, [(p, pdbqt_dir) for p in sdfs])
    bad = [lid for lid, ok in res if not ok]
    print(f"  成功 {len(res) - len(bad)}/{len(res)}" + (f"  失败 {bad[:8]}" if bad else ""))
    print(f"[prep_ligands] 完成，pdbqt 在 {pdbqt_dir.relative_to(cfg.ROOT)}/")


if __name__ == "__main__":
    sys.exit(main())
