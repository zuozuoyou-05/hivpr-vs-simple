"""方法学验证：位姿、富集、相关性三件事，配上图和 metrics.json。

A. 位姿复现  晶体配体重新对接一遍，对称校正后的重原子 RMSD，判据 2.0 A
B. 判别富集  阳性 vs 实测阴性 / 诱饵 / 两者合起来，ROC-AUC、AUPRC、EF、BEDROC
C. 相关性    只在阳性内部，打分与实测 pActivity 的 Spearman / Pearson / R² / RMSE
              另外算一份按重原子数归一化的，用来判断正相关是不是「分子越大分越好」

默认只有 Vina。如果 results/interformer_scores.csv 也在（用 run_interformer.sh 跑出来的），
会把它一并算进来，而且对比时只取两个引擎都有打分的分子，保证是配对的。

产物 results/metrics.json、results/figures/*.png
"""

import json
import sys

import cfg
import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem import rdMolAlign  # noqa: E402
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve  # noqa: E402

RDLogger.DisableLog("rdApp.*")
CJK = cfg.cjk_font()
RED, GREY = "#A32D2D", "#5F5E5A"
# 引擎 -> (打分列, 方向, 颜色)。-1 表示越负越强
ENGINES = {
    "Vina": ("vina_affinity", -1, "#185FA5"),
    "Interformer": ("pred_pIC50", +1, "#0F6E56"),
}
GROUP_LABEL = {"positive": "阳性", "negative": "实测阴性", "decoy": "诱饵"}


def enrichment_factor(y, score, frac):
    """EF@frac：打分最高的前 frac 比例里，富集了阳性多少倍。"""
    n, n_act = len(y), int(y.sum())
    if n_act == 0:
        return float("nan")
    k = max(1, int(round(n * frac)))
    hits = int(y[np.argsort(-score)[:k]].sum())
    return (hits / n_act) / (k / n)


def bedroc(y, score, alpha=20.0):
    """BEDROC：给排在前面的命中按指数加权，越靠前权重越大。"""
    n, n_act = len(y), int(y.sum())
    if n_act == 0 or n_act == n:
        return float("nan")
    hit_ranks = np.arange(1, n + 1)[y[np.argsort(-score)] == 1]
    rie = float(np.sum(np.exp(-alpha * hit_ranks / n)) / n_act)
    best = np.arange(1, n_act + 1)
    return rie / float(np.sum(np.exp(-alpha * best / n)) / n_act)


def discrimination(y, score):
    """score 已统一成越大结合越强。"""
    return {
        "n": int(len(y)), "n_positive": int(y.sum()),
        "roc_auc": round(float(roc_auc_score(y, score)), 3),
        "auprc": round(float(average_precision_score(y, score)), 3),
        "ef1": round(enrichment_factor(y, score, 0.01), 3),
        "ef5": round(enrichment_factor(y, score, 0.05), 3),
        "ef10": round(enrichment_factor(y, score, 0.10), 3),
        "bedroc20": round(bedroc(y, score, 20.0), 3),
    }


def correlation(x, y):
    """x 已统一成越大越强，y 是实测 pActivity。"""
    slope, intercept, r, _, _ = stats.linregress(x, y)
    rmse = float(np.sqrt(np.mean((y - (slope * x + intercept)) ** 2)))
    return {
        "n": int(len(x)),
        "spearman_rho": round(float(stats.spearmanr(x, y).statistic), 3),
        "spearman_p": round(float(stats.spearmanr(x, y).pvalue), 4),
        "pearson_r": round(float(stats.pearsonr(x, y).statistic), 3),
        "r2": round(float(r ** 2), 3),
        "rmse_pactivity": round(rmse, 3),
    }


def redock_rmsd(pose_pdbqt, ref_sdf):
    """重对接位姿与晶体构象的对称校正重原子 RMSD。"""
    if not (pose_pdbqt.exists() and ref_sdf.exists()):
        return None
    from meeko import PDBQTMolecule, RDKitMolCreate

    pm = PDBQTMolecule.from_file(str(pose_pdbqt), skip_typing=True)
    ref = Chem.RemoveHs(Chem.MolFromMolFile(str(ref_sdf), removeHs=False, sanitize=True))
    got = Chem.RemoveHs(Chem.Mol(RDKitMolCreate.from_pdbqt_mol(pm)[0]))
    return round(float(rdMolAlign.GetBestRMS(got, ref)), 3)


def plot_distributions(df, engines, out):
    fig, axes = plt.subplots(1, len(engines), figsize=(4.4 * len(engines), 4.2), squeeze=False)
    for ax, (name, (col, sign, color)) in zip(axes[0], engines.items()):
        data, labels = [], []
        for key, label in GROUP_LABEL.items():
            sub = df.loc[df.group == key, col].dropna().values * sign
            if len(sub):
                data.append(sub)
                labels.append(f"{label}\nn={len(sub)}")
        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, widths=0.55)
        for patch, c in zip(bp["boxes"], ["#E6F1FB", "#F1EFE8", "#D3D1C7"]):
            patch.set_facecolor(c)
            patch.set_edgecolor(GREY)
        for med in bp["medians"]:
            med.set_color(RED)
        ax.set_title(name, fontsize=11)
        ax.set_ylabel("打分（越大结合越强）" if CJK else "score (higher = stronger)")
        ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    fig.suptitle("三组分子的打分分布" if CJK else "Score by group", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "score_distributions.png", dpi=160)
    plt.close(fig)


def plot_roc_efficiency(df, engines, out):
    """左 ROC，右富集因子曲线。每个引擎两组对比。"""
    pairs = [("阳性 vs 实测阴性", ["positive", "negative"]),
             ("阳性 vs 诱饵", ["positive", "decoy"])]
    styles = ["-", "--"]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
    fracs = np.arange(0.01, 0.51, 0.01)
    for name, (col, sign, color) in engines.items():
        for (label, groups), ls in zip(pairs, styles):
            sub = df[df.group.isin(groups)].dropna(subset=[col])
            y = (sub.group == "positive").astype(int).values
            # 某个引擎可能压根没跑那一组分子，或者只剩一类，这种对比算不出 AUC，跳过
            if len(sub) < 10 or not 0 < y.sum() < len(y):
                continue
            score = sub[col].values * sign
            fpr, tpr, _ = roc_curve(y, score)
            tag = f"{name} {label}"
            axes[0].plot(fpr, tpr, linewidth=1.5, color=color, linestyle=ls,
                         label=f"{tag}  {roc_auc_score(y, score):.3f}")
            axes[1].plot(fracs * 100, [enrichment_factor(y, score, f) for f in fracs],
                         linewidth=1.5, color=color, linestyle=ls, label=tag)
    axes[0].plot([0, 1], [0, 1], "--", color="#B4B2A9", linewidth=0.8)
    axes[0].set_xlabel("假阳性率" if CJK else "FPR")
    axes[0].set_ylabel("真阳性率" if CJK else "TPR")
    axes[0].set_title("ROC 曲线" if CJK else "ROC", fontsize=11)
    axes[1].axhline(1.0, linestyle="--", color="#B4B2A9", linewidth=0.8)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("取打分最高的前 x%" if CJK else "top x%")
    axes[1].set_ylabel("富集因子" if CJK else "EF")
    axes[1].set_title("富集曲线" if CJK else "Enrichment", fontsize=11)
    for ax in axes:
        ax.legend(fontsize=7.5, frameon=False)
        ax.grid(alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(out / "roc_enrichment.png", dpi=160)
    plt.close(fig)


def plot_correlation(pos, engines, out):
    fig, axes = plt.subplots(1, len(engines), figsize=(4.6 * len(engines), 4.4), squeeze=False)
    for ax, (name, (col, sign, color)) in zip(axes[0], engines.items()):
        sub = pos.dropna(subset=[col])
        x, y = sub[col].values * sign, sub["p_activity"].values
        ax.scatter(x, y, s=28, color=color, alpha=0.75, edgecolor="white", linewidth=0.5)
        slope, intercept, r, _, _ = stats.linregress(x, y)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, slope * xs + intercept, color=RED, linewidth=1.3)
        rho = stats.spearmanr(x, y)
        ax.text(0.04, 0.94, f"Spearman ρ={rho.statistic:.2f}\nR²={r ** 2:.2f}\nn={len(x)}",
                transform=ax.transAxes, va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.35", facecolor="#F1EFE8",
                          edgecolor="#D3D1C7", linewidth=0.5))
        ax.set_title(name, fontsize=11)
        ax.set_xlabel("打分（越大结合越强）" if CJK else "score (higher = stronger)")
        ax.set_ylabel("实测 pActivity")
        ax.grid(alpha=0.25, linewidth=0.5)
    fig.suptitle("阳性内部的打分-活性关系" if CJK else "Score vs activity", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "correlation.png", dpi=160)
    plt.close(fig)


def main():
    cfgd = cfg.load()
    tag = cfgd["pdb_tag"]
    proc = cfg.path("data", "processed")
    res = cfg.out("results")
    figs = cfg.out("results", "figures")

    merged = pd.read_csv(proc / "dataset.csv").merge(
        pd.read_csv(res / "vina_scores.csv"), on="ligand_id", how="left")

    # Interformer 是可选的第二把尺子，跑了才并进来
    if_csv = res / "interformer_scores.csv"
    if if_csv.exists():
        merged = merged.merge(pd.read_csv(if_csv), on="ligand_id", how="left")
    engines = {n: spec for n, spec in ENGINES.items()
               if spec[0] in merged.columns and merged[spec[0]].notna().any()}
    if not engines:
        raise SystemExit("没有可用打分，先跑 batchdock.sh + extract_scores.py")
    if len(engines) > 1:
        coverage = merged[[spec[0] for spec in engines.values()]]
        print(f"[analyze] 两个引擎都覆盖到的分子 {int(coverage.notna().all(axis=1).sum())} 个")

    metrics = {
        "target": cfgd["target_name"],
        "engines": list(engines),
        "n_molecules": {g: int((merged.group == g).sum()) for g in GROUP_LABEL},
    }

    # ---- A 位姿复现 ----
    pose = {"Vina": {
        "reference_ligand": f"{tag}_ref",
        "redock_rmsd": redock_rmsd(cfg.path("vina_outs", f"{tag}_ref.pdbqt"),
                                   proc / "sdf" / f"{tag}_ref.sdf"),
    }}
    if pose["Vina"]["redock_rmsd"] is not None:
        pose["Vina"]["pass_2A"] = pose["Vina"]["redock_rmsd"] <= 2.0
    if_json = res / "pose_reproduction.json"
    if if_json.exists():
        p = json.loads(if_json.read_text(encoding="utf-8"))
        # Monte Carlo 采样会把输入构象原样留成一条候选位姿，而这里输入的就是晶体构象，
        # 所以「最小 RMSD」恒等于 0，没有信息量。要看它剔掉起始构象后自己生成了什么
        pose["Interformer"] = {
            "best_sampled_rmsd": p.get("best_sampled_rmsd"),
            "median_sampled_rmsd": p.get("median_sampled_rmsd"),
            "n_within_2A": p.get("n_sampled_within_2A"),
            "n_poses_generated": p.get("n_poses_generated"),
        }
    metrics["pose_reproduction"] = pose

    # ---- B 判别富集。每个引擎按自己覆盖到的分子算，n 单独标出来 ----
    # Interformer 这边没跑诱饵，所以它进不了「阳性 vs 诱饵」那一行，n 会露馅，不用假装对齐
    comparisons = {
        "pos_vs_neg": merged[merged.group.isin(["positive", "negative"])],
        "pos_vs_decoy": merged[merged.group.isin(["positive", "decoy"])],
        "pos_vs_all": merged,
    }
    disc = {}
    for comp, sub in comparisons.items():
        per_engine = {}
        for n, (col, sign, _) in engines.items():
            s = sub.dropna(subset=[col])
            if s.empty:
                continue
            y = (s.group == "positive").astype(int).values
            if not 0 < y.sum() < len(y):
                continue
            per_engine[n] = discrimination(y, s[col].values * sign)
        if per_engine:
            disc[comp] = per_engine
    metrics["discrimination"] = disc

    # ---- C 相关性，只在阳性内部 ----
    positives = merged[merged.group == "positive"]
    corr = {}
    for n, (col, sign, _) in engines.items():
        s = positives.dropna(subset=[col])
        score = s[col].values * sign
        c = correlation(score, s["p_activity"].values)
        # 除以重原子数：如果正相关是「分子越大分越好」撑起来的，这个数会塌下来
        c["spearman_rho_size_normalised"] = round(float(
            stats.spearmanr(score / s["heavy_atoms"].values, s["p_activity"].values).statistic), 3)
        corr[n] = c
    metrics["correlation"] = corr

    # ---- 三组的尺码。尺码不匹配时 AUC 会被尺寸效应污染，摆出来自己看 ----
    metrics["group_properties"] = {
        g: {"n": int((merged.group == g).sum()),
            "mean_heavy_atoms": round(float(merged.loc[merged.group == g, "heavy_atoms"].mean()), 1),
            "mean_mw": round(float(merged.loc[merged.group == g, "mw"].mean()), 1)}
        for g in GROUP_LABEL}
    metrics["note"] = "打分统一成数值越大表示结合越强之后再算指标；富集对比只用共同覆盖的分子"

    (res / "metrics.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")
    merged.to_csv(res / "scores_merged.csv", index=False)

    plot_distributions(merged, engines, figs)
    plot_roc_efficiency(merged, engines, figs)
    plot_correlation(positives, engines, figs)

    # ---- 控制台摘要 ----
    print("\n=== A 位姿复现（重对接 RMSD，判据 2.0 A）===")
    r = pose["Vina"]["redock_rmsd"]
    print(f"  Vina  {tag}_ref  RMSD = {r if r is not None else 'n/a'} A  "
          f"{'通过' if pose['Vina'].get('pass_2A') else '未通过'}")
    if "Interformer" in pose:
        v = pose["Interformer"]
        print(f"  Interformer 剔掉起始构象后最好 {v['best_sampled_rmsd']} A、"
              f"中位 {v['median_sampled_rmsd']} A，{v['n_poses_generated']} 条里 "
              f"{v['n_within_2A']} 条在 2 A 以内")

    print("\n=== B 判别富集 ===")
    for comp, per_engine in disc.items():
        print(f"  {comp}")
        print(f"    {'引擎':<14s}{'n':>6s}{'ROC-AUC':>9s}{'AUPRC':>8s}{'EF1%':>7s}{'EF5%':>7s}{'BEDROC':>9s}")
        for name, m in per_engine.items():
            print(f"    {name:<14s}{m['n']:>6d}{m['roc_auc']:>9.3f}{m['auprc']:>8.3f}"
                  f"{m['ef1']:>7.2f}{m['ef5']:>7.2f}{m['bedroc20']:>9.3f}")

    print("\n=== C 相关性（阳性）===")
    for name, m in corr.items():
        print(f"  {name:<14s} n={m['n']:<3d} Spearman ρ={m['spearman_rho']:>6.3f}"
              f" (p={m['spearman_p']:.3f})  R²={m['r2']:.3f}  RMSE={m['rmse_pactivity']:.3f}"
              f"  按重原子归一后 ρ={m['spearman_rho_size_normalised']:.3f}")

    print("\n=== 三组的尺码（判断尺寸混杂用）===")
    for g, v in metrics["group_properties"].items():
        print(f"  {GROUP_LABEL[g]:<6s} n={v['n']:<4d} 平均重原子 {v['mean_heavy_atoms']:>5.1f}"
              f"  平均 MW {v['mean_mw']:>6.1f}")
    print(f"\n[analyze] metrics.json 和 {len(list(figs.glob('*.png')))} 张图写到 results/")


if __name__ == "__main__":
    sys.exit(main())
