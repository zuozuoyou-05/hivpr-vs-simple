#!/usr/bin/env bash
# 可选：同一批配体再让 Interformer 跑一遍，跟 Vina 做对照。
#
# Interformer 官方仓库自带四步（生成查询表 → 能量模型 → Monte Carlo 采样 → 亲和力打分），
# 这个脚本只负责把四步串起来，再把每个配体的最优构象汇总到 results/interformer_scores.csv。
# 之后 analyze.py 会自己发现它，并把它一起算进对比。
#
# 用法
#   IF_REPO=/mnt/d/Interformer bash scripts/run_interformer.sh
#   IF_REPO=/mnt/d/Interformer IF_WORK=/home/zuoyou/ifwork bash scripts/run_interformer.sh 40
#     最后一个参数是只跑前 N 个配体（想快点见效就先用它）
#
# 注：这个分支在本次简化里没有重跑，四步的命令行参数是从完整版搬过来的（那边跑通过）。
#     results/ 里现有的 interformer_scores.csv 也是完整版跑出来的，配体集跟这里一致。
#
# 两个坑，都踩过
#   1. 第 2、4 步必须带 -debug。官方 parser 里 n_jobs 默认 70、DataLoader worker 默认 5，
#      那是给大内存机器的。本机 7 GB 的 WSL 上不加必被内核 OOM-killer 收走，而且它优先杀
#      systemd，表现出来是整个发行版重启、后台任务一起消失，不报错，很难查。
#   2. -ensemble 的通配符不能加引号。parser 里它是 nargs='+'，靠 shell 展开成 4 个模型目录；
#      一旦加引号就只剩字面量 model*，实际只跑 1 个模型，同样不报错。
#
# 耗时：Monte Carlo 采样是单进程 CPU，约 1.5 分钟/配体。201 个配体要 5 小时以上。

set -euo pipefail

cd "$(dirname "$0")/.." || exit 1

PY=${PYTHON:-python}
IF_REPO=${IF_REPO:-/mnt/d/Interformer}
IF_ENV=${IF_ENV:-interformer}
WORK=${IF_WORK:-work/interformer}
LIMIT=${1:-0}
TAG=$($PY scripts/cfg.py pdb_tag)

[ -d "$IF_REPO" ] || { echo "找不到 Interformer 仓库：$IF_REPO（用 IF_REPO=... 指路）"; exit 1; }

run() {
  conda run -n "$IF_ENV" --no-capture-output bash -c "cd '$IF_REPO' && $*"
}

echo "[interformer] 0/5 准备输入"
mkdir -p "$WORK"/{ligand,pocket,uff,infer}
cp "data/processed/sdf/${TAG}_ref.sdf" "$WORK/ligand/${TAG}_docked.sdf"   # 晶体构象，当作口袋的对照
cp data/processed/receptor_h.pdb "$WORK/pocket/${TAG}_pocket.pdb"        # 加氢受体
# SDF 本身是纯文本，直接拼就行。参考配体放第一条：MC 采样按索引顺序跑，它最先出结果
cat "data/processed/sdf/${TAG}_ref.sdf" data/processed/ligands_3d.sdf > "$WORK/uff/${TAG}_uff.sdf"
if [ "$LIMIT" -gt 0 ] 2>/dev/null; then
  awk -v n="$LIMIT" '/^\$\$\$\$/{c++; print; if (c>=n) exit; next} {print}' \
      "$WORK/uff/${TAG}_uff.sdf" > "$WORK/uff/trim.sdf"
  mv "$WORK/uff/trim.sdf" "$WORK/uff/${TAG}_uff.sdf"
fi
echo "  配体数 $(( $(grep -c '^\$\$\$\$' "$WORK/uff/${TAG}_uff.sdf") ))"

echo "[interformer] 1/5 生成查询表"
run "python tools/inference/inter_sdf2csv.py \"$WORK/uff/${TAG}_uff.sdf\" 1"

echo "[interformer] 2/5 能量模型预测原子对交互"
run "PYTHONPATH=interformer/ python inference.py -test_csv \"$WORK/uff/${TAG}_uff_infer.csv\" \
     -work_path \"$WORK/\" -ensemble checkpoints/v0.2_energy_model -gpus 1 -batch_size 1 \
     -posfix '*val_loss*' -energy_output_folder \"$WORK/energy_VS\" -uff_as_ligand -debug -reload"

echo "[interformer] 3/5 Monte Carlo 采样生成构象"
run "cp \"$WORK/uff/${TAG}_uff.sdf\" \"$WORK/energy_VS/uff/\" && \
     OMP_NUM_THREADS=8,8 python docking/reconstruct_ligands.py -y --cwd \"$WORK/energy_VS\" \
       --find_all --uff_folder uff find && \
     python docking/reconstruct_ligands.py -y --cwd \"$WORK/energy_VS\" --find_all stat"

echo "[interformer] 4/5 PoseScore + 亲和力打分"
run "rm -rf \"$WORK/infer\" && mkdir -p \"$WORK/infer\" && \
     cp -r \"$WORK/energy_VS/ligand_reconstructing/\"*.sdf \"$WORK/infer/\" && \
     python tools/inference/inter_sdf2csv.py \"$WORK/infer/${TAG}_docked.sdf\" 0 && \
     PYTHONPATH=interformer/ python inference.py -test_csv \"$WORK/infer/${TAG}_docked_infer.csv\" \
       -work_path \"$WORK/\" -ligand_folder infer/ \
       -ensemble checkpoints/v0.2_affinity_model/model* -use_ff_ligands '' -vs -gpus 1 \
       -batch_size 4 -posfix '*val_loss*' --pose_sel True -debug"

echo "[interformer] 5/5 汇总打分"
IF_REPO="$IF_REPO" TAG="$TAG" $PY - <<'PY'
import os
from pathlib import Path

import pandas as pd

repo = Path(os.environ["IF_REPO"])
tag = os.environ["TAG"]
src = sorted(repo.glob(f"result/{tag}_docked_infer_ensemble.csv")) or \
      sorted((repo / "result").glob("*docked_infer_ensemble.csv"))
if not src:
    raise SystemExit(f"没找到打分结果，检查 {repo}/result/")
df = pd.read_csv(src[-1])
print(f"  读到 {src[-1].name}，{len(df)} 行")

# 每个配体留下最好的那条构象：pred_pIC50 最高，同分再看 PoseScore
best = df.sort_values(["pred_pIC50", "pred_pose"], ascending=[False, False]) \
         .groupby("Molecule ID", as_index=False).first() \
         .rename(columns={"Molecule ID": "ligand_id"})
cols = ["ligand_id"] + [c for c in ["energy", "pred_pIC50", "pred_pose", "num_torsions", "rmsd"]
                        if c in best.columns]
out = Path("results/interformer_scores.csv")
best[cols].to_csv(out, index=False)
print(f"  写出 {out}  {len(best)} 个配体")
PY

echo "[interformer] 完成。下一步 python scripts/analyze.py"
