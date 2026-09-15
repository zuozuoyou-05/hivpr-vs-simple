"""把 vina_outs/*.log 里的打分抠出来，合成一张表。

跟 ChameleonBinder 的 extract_score_and_pose_confidence.py 一个路子：先找到
'-----+...' 那行表头，下面每一行是一个构象，第 0 列是序号、第 1 列是亲和力。
这里顺手做了阈值统计，看每个组里有多少配体打到了 -7 kcal/mol 以内。

产物 results/vina_scores.csv
  ligand_id, vina_affinity, vina_affinity_2, n_poses, seconds, hit_7
  vina_affinity     最优构象的打分
  vina_affinity_2   次优构象的打分，用来判断能量面平不平
  hit_7             打分 <= -7 kcal/mol 记为 True
"""

import csv
import re
import sys
from pathlib import Path

import cfg

SEPARATOR = re.compile(r"^[-+]{5,}$")


def parse_log(path):
    """从一份 vina log 里取打分。返回 (各构象打分列表, 耗时秒)。"""
    scores, seconds, in_table = [], None, False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if SEPARATOR.match(line.strip()):
            in_table = True
            continue
        if line.startswith("# seconds"):
            seconds = float(line.split()[-1])
            continue
        if not in_table or not line.strip() or not line.strip()[0].isdigit():
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                scores.append(float(parts[1]))
            except ValueError:
                pass
    return scores, seconds


def main():
    cfgd = cfg.load()
    tag = cfgd["pdb_tag"]
    vina_outs = cfg.path("vina_outs")
    res = cfg.out("results")

    rows = []
    for log in sorted(vina_outs.glob("*.log")):
        scores, seconds = parse_log(log)
        if not scores:
            continue
        rows.append({
            "ligand_id": log.stem,
            "vina_affinity": round(min(scores), 3),
            "vina_affinity_2": round(sorted(scores)[1], 3) if len(scores) > 1 else "",
            "n_poses": len(scores),
            "seconds": seconds if seconds is not None else "",
            "hit_7": min(scores) <= -7.0,
        })
    if not rows:
        raise SystemExit(f"{vina_outs} 里没有可解析的 log，先把 batchdock.sh 跑完")

    rows.sort(key=lambda r: (r["ligand_id"] != f"{tag}_ref", r["ligand_id"]))
    out_csv = res / "vina_scores.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[extract_scores] {len(rows)} 个配体的打分写到 {out_csv.relative_to(cfg.ROOT)}")

    # 阈值统计。诱饵那一组最值得看：它们跟阳性性质匹配，最像真实筛选库
    ds_path = cfg.path("data", "processed", "dataset.csv")
    group = {}
    if ds_path.exists():
        import pandas as pd
        ds = pd.read_csv(ds_path)
        group = dict(zip(ds["ligand_id"], ds["group"]))
    names = {"positive": "阳性", "negative": "实测阴性", "decoy": "诱饵"}
    print("  打分 <= -7 kcal/mol 的比例")
    for key, label in names.items():
        sub = [r for r in rows if group.get(r["ligand_id"]) == key]
        if sub:
            n_hit = sum(1 for r in sub if r["hit_7"])
            print(f"    {label:<6s} {n_hit}/{len(sub)}  = {n_hit / len(sub):.0%}")
    ref = [r for r in rows if r["ligand_id"] == f"{tag}_ref"]
    if ref:
        print(f"    参考配体 {tag}_ref  {ref[0]['vina_affinity']:.2f} kcal/mol")


if __name__ == "__main__":
    sys.exit(main())
