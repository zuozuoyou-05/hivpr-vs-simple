"""准备受体：从晶体结构抽出蛋白、加氢、生成 pdbqt，顺便定出对接盒子。

受体用 RCSB 1XL2（1.50 A，HIV-1 蛋白酶 + 抑制剂 189 的复合物）。走四步
  1. 抽出 A/B 两条链，丢掉所有 HETATM（配体 189、甘油、氯离子、水），只保留蛋白质
  2. reduce 加氢，按局部环境翻转 Asn/Gln/His
  3. meeko 生成 AutoDock pdbqt，极性氢 + Gasteiger 电荷
  4. 由晶体配体的包围盒推出盒子，写 data/processed/box.json

为什么不直接用 DUD-E 自带的 receptor.pdb：那份文件把残基名改过（GLY→GLZ、ASP→ASQ、
质子化的 ASP→ASH），而且 element 列是空的，meeko 匹配模板时会失败。用 1XL2 原件更省事。

产物
  data/processed/receptor.pdbqt  对接受体（刚性）
  data/processed/box.txt         盒子的中心与尺寸，人看的
  data/processed/box.json        同上，脚本读的
"""

import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import cfg

# reduce 和 obabel 在 Interformer 那个环境里，本机没有系统级安装
REDUCE = Path(os.environ.get("REDUCE_BIN", Path.home() / "miniconda3/envs/interformer/bin/reduce"))
OBABEL = Path(os.environ.get("OBABEL_BIN", Path.home() / "miniconda3/envs/interformer/bin/obabel"))


def guess_element(atom_name):
    """PDB 原子名推元素符号。蛋白里只有 C/N/O/S/H。"""
    core = atom_name.strip().lstrip("0123456789")
    if not core:
        return ""
    return "H" if core[0] == "H" else core[0].upper()


def extract_protein(src, dst, chains):
    """抽出指定链的蛋白原子，补上 element 列（P7-P8，PDB 里经常是空的）。"""
    n = 0
    lines = []
    for line in open(src, encoding="utf-8", errors="replace"):
        if not line.startswith("ATOM"):
            continue
        if line[21] not in chains or line[16] not in (" ", "A"):
            continue
        line = line.rstrip("\n").ljust(80)
        elem = line[76:78].strip() or guess_element(line[12:16])
        if elem == "H":
            continue
        lines.append(line[:16] + " " + line[17:76] + f"{elem:>2}" + "\n")
        n += 1
    dst.write_text("".join(lines), encoding="utf-8")
    return n


def fix_pdb_elements(src, dst):
    """把 element 列写回去。reduce 的输出偶尔缺这一列。"""
    out = []
    for line in open(src, encoding="utf-8", errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            line = line.rstrip("\n").ljust(80)
            elem = line[76:78].strip() or guess_element(line[12:16])
            line = line[:76] + f"{elem:>2}" + line[78:]
        out.append(line if line.endswith("\n") else line + "\n")
    dst.write_text("".join(out), encoding="utf-8")


def complete_histidine(path):
    """给咪唑环上两个氮都没氢的组氨酸补一个 HE2。

    reduce 碰到没有明确氢键环境的组氨酸会保持未质子化，meeko 于是在 HID/HIE/HIP
    之间判不出来（报 missing 和 excess 并列）。这里统一按 HIE 补 HE2。新原子必须插在
    该残基原子块的最后一行后面，否则 meeko 会说 residues interrupted。
    """
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    his = {}
    for idx, line in enumerate(lines):
        if not line.startswith("ATOM") or line[17:20].strip() != "HIS":
            continue
        rec = his.setdefault((line[21], line[22:27].strip()), {"atoms": {}, "last": idx})
        rec["atoms"][line[12:16].strip()] = [float(line[30:38]), float(line[39:46]), float(line[46:54])]
        rec["last"] = idx

    inserts = []
    for (chain, resnum), rec in his.items():
        a = rec["atoms"]
        if "HD1" in a or "HE2" in a or not {"NE2", "CD2", "CE1"} <= set(a):
            continue
        # 氢沿着 NE2 背向咪唑环心伸出去，键长 1.01 A
        v = [a["CD2"][k] + a["CE1"][k] - 2 * a["NE2"][k] for k in range(3)]
        norm = math.sqrt(sum(c * c for c in v)) or 1.0
        h = [a["NE2"][k] - 1.01 * v[k] / norm for k in range(3)]
        newline = ("ATOM  " + "    0" + " " + " HE2" + " " + "HIS" + " " + chain
                   + f"{int(resnum):4d}" + " " + "   "
                   + f"{h[0]:8.3f}{h[1]:8.3f}{h[2]:8.3f}"
                   + "  1.00" + "  0.00" + " " * 10 + " H\n")
        inserts.append((rec["last"], newline))

    for idx, newline in sorted(inserts, reverse=True):
        lines.insert(idx + 1, newline)
    if inserts:
        path.write_text("".join(lines), encoding="utf-8")
    return len(inserts)


def make_box(mol2_path, proc, pad, min_size):
    """盒子 = 晶体配体坐标的包围盒，双侧各加 pad。"""
    coords, in_atom = [], False
    for line in open(mol2_path, encoding="utf-8", errors="replace"):
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_atom = False
            continue
        if in_atom and line.strip():
            f = line.split()
            coords.append((float(f[2]), float(f[3]), float(f[4])))
    assert coords, f"从 {mol2_path} 读不到配体坐标"
    xs, ys, zs = zip(*coords)
    center = [round((min(v) + max(v)) / 2, 3) for v in (xs, ys, zs)]
    size = [round(max(max(v) - min(v) + 2 * pad, min_size), 3) for v in (xs, ys, zs)]
    box = {"center": center, "size": size}
    (proc / "box.json").write_text(json.dumps(box, indent=2), encoding="utf-8")
    (proc / "box.txt").write_text(
        f"center_x = {center[0]}\ncenter_y = {center[1]}\ncenter_z = {center[2]}\n"
        f"size_x = {size[0]}\nsize_y = {size[1]}\nsize_z = {size[2]}\n", encoding="utf-8")
    return box


def run_meeko(rec_h, proc, center, size, max_retry=30):
    """跑 mk_prepare_receptor。组氨酸互变异构判不出来时，报错信息里会带 residue_key，
    解析出来显式指定成 HIE 再重试。"""
    fixes = {}
    res = None
    for _ in range(max_retry):
        cmd = ["mk_prepare_receptor.py", "--read_pdb", str(rec_h),
               "-o", str(proc / "receptor"),
               "-p", str(proc / "receptor.pdbqt"),
               "-j", str(proc / "receptor.json"),
               "-v", str(proc / "box_vina.txt"),
               "--box_center", *[str(c) for c in center],
               "--box_size", *[str(s) for s in size],
               "--charge_model", "gasteiger"]
        for key, tmpl in fixes.items():
            cmd += ["-n", f"{key}={tmpl}"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0 and (proc / "receptor.pdbqt").exists():
            return True, res, fixes
        m = re.search(r"residue_key='([^']+)'", res.stderr)
        if not m or m.group(1) in fixes:
            return False, res, fixes
        fixes[m.group(1)] = "HIE"
    return False, res, fixes


def main():
    cfgd = cfg.load()
    proc = cfg.out("data", "processed")
    rec_src = cfg.path(cfgd["receptor_pdb"])
    rec_noH = proc / "receptor_noH.pdb"
    rec_h = proc / "receptor_h.pdb"

    print("[prep_receptor] 1/4 抽出蛋白")
    n = extract_protein(rec_src, rec_noH, cfg.lst(cfgd, "receptor_chains"))
    print(f"  {rec_src.name} 的 {'/'.join(cfg.lst(cfgd, 'receptor_chains'))} 链，{n} 个原子")

    print("[prep_receptor] 2/4 reduce 加氢")
    assert REDUCE.exists(), f"找不到 reduce：{REDUCE}（用 REDUCE_BIN 环境变量指路）"
    with open(proc / "reduce_tmp.pdb", "w", encoding="utf-8") as fh, open(proc / "reduce.log", "w") as log:
        subprocess.run([str(REDUCE), str(rec_noH)], stdout=fh, stderr=log, check=True)
    fix_pdb_elements(proc / "reduce_tmp.pdb", rec_h)
    (proc / "reduce_tmp.pdb").unlink(missing_ok=True)
    n_his = complete_histidine(rec_h)
    n_h = sum(1 for l in open(rec_h) if l.startswith("ATOM") and l[76:78].strip() == "H")
    print(f"  加氢完成，共 {sum(1 for l in open(rec_h) if l.startswith('ATOM'))} 个原子"
          f"（氢 {n_h}）" + (f"，补了 {n_his} 个组氨酸的 HE2" if n_his else ""))

    print("[prep_receptor] 3/4 定盒子")
    box = make_box(cfg.path(cfgd["crystal_ligand"]), proc,
                   cfg.f(cfgd, "box_padding"), cfg.f(cfgd, "box_min_size"))
    print(f"  center={box['center']}  size={box['size']}")

    print("[prep_receptor] 4/4 meeko 生成 pdbqt")
    (proc / "receptor.pdbqt").unlink(missing_ok=True)
    ok, res, fixes = run_meeko(rec_h, proc, box["center"], box["size"])
    if fixes:
        print(f"  组氨酸互变异构体显式指定了 {len(fixes)} 处 " + ", ".join(f"{k}={v}" for k, v in fixes.items()))
    if not ok:
        print("  meeko 失败，退回 openbabel -xr。stderr 末尾：\n" + res.stderr[-800:])
        subprocess.run([str(OBABEL), str(rec_h), "-O", str(proc / "receptor.pdbqt"), "-xr"],
                       check=True, capture_output=True)
    else:
        print("  meeko 成功")
    n_atom = sum(1 for l in open(proc / "receptor.pdbqt") if l.startswith(("ATOM", "HETATM")))
    print(f"[prep_receptor] 完成，receptor.pdbqt {n_atom} 个原子")


if __name__ == "__main__":
    sys.exit(main())
