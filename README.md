# HIV-1 蛋白酶批量对接
## 它做什么

拿 DUD-E 的 hivpr 当靶点。200 个分子加一个晶体配体全部丢给 Vina 跑一遍，然后回答三个问题。

**位姿复现。** 把晶体配体重新对一次，看它能不能回到晶体里那个位置。判据是重原子 RMSD 小于 2 Å。

**判别富集。** 阳性、实测阴性、诱饵混在一个库里，看打分能不能把真活性物排到前面。指标用 ROC-AUC、
AUPRC、EF 和 BEDROC。

**打分与活性的相关性。** 只在 30 个阳性内部算，看对接打分跟实测 pActivity 对得上多少。另外算一份按重原子数
归一化的版本，用来判断这个正相关里有多少是「分子越大分越好」撑起来的。

## 数据

| 组 | 数量 | 来源 | 用途 |
| --- | --- | --- | --- |
| 阳性 | 30 | DUD-E `actives_nM_chembl`，只取实测的 Ki / IC50 / Kd | 当正样本，也当回归样本 |
| 实测阴性 | 40 | DUD-E `inactives_nM_chembl`，原文写的是 > 40 µM 仍无抑制 | 查假阳性 |
| 诱饵 | 130 | DUD-E `decoys_final`，分子量、logP、氢键数跟阳性配过对 | 模拟真实筛选库 |
| 参考配体 | 1 | 1XL2 晶体里的配体 189 | 定盒子，做重对接 |

阳性按 pActivity 分六档抽样，每档都抽，免得清一色全是强抑制剂。三组之间再按 Murcko 骨架去冗余。
pActivity 取 9 − log10(数值 nM)，所以 Ki = 1 nM 就是 9。

诱饵和实测阴性得分开报。诱饵是「看着像药但大概不结合」，实测阴性是「看着像药、测过、确实不结合」。
一个只在阳性对诱饵上好看、一换成实测阴性就塌掉的打分函数，说明它认的只是理化性质，不是相互作用。
这是虚筛里最常见的伪成功，所以两个数都要有。

受体用 RCSB 的 [1XL2](https://www.rcsb.org/structure/1XL2)（1.50 Å），也就是 DUD-E hivpr 的原始来源。
自己的 `receptor.pdb` 我用不了，它把残基名改过（GLY→GLZ、ASP→ASQ），element 列还是空的，
meeko 匹配模板会失败。直接从 1XL2 抽 A/B 链更省事。

## 跑起来

需要一个能 import vina、rdkit、meeko 的环境

```bash
conda activate vsbench
bash scripts/run_all.sh
```


```bash
bash scripts/run_all.sh 20
```

| 脚本 | 干什么 |
| --- | --- |
| `prep_ligands.py` | 解析 DUD-E 的 .ism，挑分子，SMILES 生成 3D 构象，转 pdbqt |
| `prep_receptor.py` | 从 1XL2 抽 A/B 链、加氢、生成受体 pdbqt，由晶体配体算出盒子 |
| `batchdock.sh` | 批量对接。格点算一次，然后逐个配体并行跑，已有的跳过 |
| `dock_one.py` | 单个配体对接，被 `batchdock.sh` 调起来 |
| `extract_scores.py` | 把每个配体的 log 汇总成 `results/vina_scores.csv` |
| `analyze.py` | 位姿、富集、相关性三件事，出图，写 `metrics.json` |
| `run_interformer.sh` | 可选，同一批配体再用 Interformer 跑一遍做对照 |

## 看结果

跑完的东西都落在 `results/`。想直接看结论，打开 `vs_results.ipynb`，三件事的表和图都摆好了，不用重跑对接。
想从零复现，`bash scripts/run_all.sh`，本机 8 进程大约 45 分钟。

## 流程

```
prep_ligands.py   DUD-E .ism → 过滤去重分层抽样 → dataset.csv
                  SMILES → ETKDGv3(seed=42) + UFF → 3D SDF → meeko → pdbqt
       │
prep_receptor.py  1XL2.pdb → 抽 A/B 链 → reduce 加氢 → meeko → receptor.pdbqt
       │                                                      └ box.txt
       ▼
batchdock.sh      affinity maps 算一次 → 8 路并行，逐配体 load_maps → vina_outs/<id>.log
       │
extract_scores.py vina_outs/*.log → results/vina_scores.csv
       │
analyze.py        位姿 RMSD + 富集指标 + 相关性 → metrics.json + figures/
```


```bash
printf '%s\n' $ligs | xargs -P "$WORKERS" -I {} bash -c '
  id=$(basename "$1" .pdbqt)
  [ -s "$OUT/$id.pdbqt" ] && exit 0        # 已经有了就跳过，中断后重跑不会白干
  $PY scripts/dock_one.py --ligand "$1" --maps "$MAPS" --out "$OUT/$id.pdbqt" ...
' _ {}
```


## 来源

- 受体 [RCSB 1XL2](https://www.rcsb.org/structure/1XL2)
- 活性与诱饵 [DUD-E](http://dude.docking.org/) 的 hivpr
- AutoDock Vina 1.2.x，Eberhardt et al., *J. Chem. Inf. Model.* 2021
- Interformer，Lai et al., *Nat. Commun.* 2024
- RDKit、Meeko、reduce (MolProbity)

MIT
