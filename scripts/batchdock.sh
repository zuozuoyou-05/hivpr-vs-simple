#!/usr/bin/env bash
# 批量对接。结构照 ChameleonBinder 的 batchdock.sh：先把格点备好，然后逐个配体跑 vina，
# 已经在 vina_outs/ 里有结果的直接跳过，所以中断后重跑不会从头再来。
#
# 用法
#   bash scripts/batchdock.sh          跑全部
#   bash scripts/batchdock.sh 20       先跑 20 个试水
#
# 跑之前先激活一个能 import vina 的环境（本机是 vsbench）。用别的解释器就
#   PYTHON=/path/to/python bash scripts/batchdock.sh

set -uo pipefail

cd "$(dirname "$0")/.." || exit 1

PY=${PYTHON:-python}
OUT=vina_outs
LIG_DIR=data/processed/pdbqt
MAPS="$OUT/receptor"

EXH=$($PY scripts/cfg.py exhaustiveness)
NPOSES=$($PY scripts/cfg.py n_poses)
WORKERS=$($PY scripts/cfg.py workers)
SEED=$($PY scripts/cfg.py vina_seed)
TAG=$($PY scripts/cfg.py pdb_tag)

mkdir -p "$OUT"

# 格点只算一次：它只跟受体和盒子有关，跟配体无关
if [ ! -f "$MAPS.maps" ]; then
  echo "[batchdock] 算 affinity maps"
  $PY scripts/dock_one.py --make-maps --maps "$MAPS" || exit 1
fi

# 配体清单。参考配体排第一个，因为最先要拿到的就是它的重对接 RMSD
ligs=$(ls "$LIG_DIR/${TAG}_ref.pdbqt" "$LIG_DIR"/*.pdbqt 2>/dev/null | sort -u)
limit=${1:-0}
if [ "$limit" -gt 0 ] 2>/dev/null; then
  ligs=$(printf '%s\n' $ligs | head -n "$limit")
fi
total=$(printf '%s\n' $ligs | grep -c .)
echo "[batchdock] 待跑 $total 个配体，并行 $WORKERS 路，exhaustiveness=$EXH"

export OUT MAPS EXH NPOSES SEED PY
printf '%s\n' $ligs | xargs -P "$WORKERS" -I {} bash -c '
  lig="$1"
  id=$(basename "$lig" .pdbqt)
  [ -s "$OUT/$id.pdbqt" ] && exit 0
  if $PY scripts/dock_one.py --ligand "$lig" --maps "$MAPS" --out "$OUT/$id.pdbqt" \
        --log "$OUT/$id.log" --exh "$EXH" --poses "$NPOSES" --seed "$SEED" >/dev/null 2>&1; then
    echo "  $id 完成"
  else
    echo "  $id 失败" >&2
  fi
' _ {}

echo "[batchdock] 结束，vina_outs/ 里现在有 $(ls "$OUT"/*.pdbqt 2>/dev/null | wc -l) 个结果"
