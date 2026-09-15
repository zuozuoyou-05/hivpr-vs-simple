#!/usr/bin/env bash
# 全流程。数据已经在 data/raw/ 里了，直接跑。
#
#   conda activate vsbench
#   bash scripts/run_all.sh
#
# 想省事只跑前 20 个配体试水：bash scripts/run_all.sh 20

set -euo pipefail

cd "$(dirname "$0")/.." || exit 1

PY=${PYTHON:-python}
LIMIT=${1:-0}

$PY scripts/prep_ligands.py
$PY scripts/prep_receptor.py
bash scripts/batchdock.sh "$LIMIT"
$PY scripts/extract_scores.py
$PY scripts/analyze.py
