"""
Plot IMP (delta_dG) vs AAR and IMP vs RMSD scatter plots for design structures.

IMP per sample = dG_gt - dG_pred (positive = improved binding over native).
dG comes from relaxed PDBs, RMSD/AAR from unrelaxed evaluation.
Merged on target code (PDB_Name).

Reads:
  r1_eval/dg_results/{model}_dg.csv
  {DESIGN_BASE}/{model}/predictions/results.csv
  evaluate/AbX_eval/imp_results/test_ground_binding.csv

Writes:
  r1_eval/plots/imp_vs_aar_all.png
  r1_eval/plots/imp_vs_rmsd_all.png
  r1_eval/plots/imp_vs_aar_{model}.png   (3 individual)
  r1_eval/plots/imp_vs_rmsd_{model}.png  (3 individual)
  stdout: correlation coefficients
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN_BASE = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_exp_data/r1_design_structures"
GT_CSV = os.path.join(PROJECT_DIR, "evaluate/AbX_eval/imp_results/test_ground_binding.csv")

MODELS = [
    "gen_trial_0",
    "gen_trial_1",
    "mfdesign_test_bagel_protein_boltz_base",
]

SHORT_NAMES = {
    "gen_trial_0": "gt0",
    "gen_trial_1": "gt1",
    "mfdesign_test_bagel_protein_boltz_base": "mfdesign_bagel",
}

CDR_PAIRS = [
    ("heavy_cdr1", "H-CDR1"),
    ("heavy_cdr2", "H-CDR2"),
    ("heavy_cdr3", "H-CDR3"),
    ("light_cdr1", "L-CDR1"),
    ("light_cdr2", "L-CDR2"),
    ("light_cdr3", "L-CDR3"),
]

COLORS = {
    "H-CDR1": "#1f77b4",
    "H-CDR2": "#ff7f0e",
    "H-CDR3": "#d62728",
    "L-CDR1": "#9467bd",
    "L-CDR2": "#8c564b",
    "L-CDR3": "#2ca02c",
}


def load_merged(model, script_dir, gt_dg):
    """Load and merge dG with RMSD/AAR data, compute delta_dG."""
    dg_csv = os.path.join(script_dir, "dg_results", f"{model}_dg.csv")
    results_csv = os.path.join(DESIGN_BASE, model, "predictions", "results.csv")

    if not os.path.exists(dg_csv) or not os.path.exists(results_csv):
        return None

    dg = pd.read_csv(dg_csv)
    results = pd.read_csv(results_csv)

    merged = results.merge(dg, left_on="code", right_on="PDB_Name", how="inner")
    merged = merged.dropna(subset=["dG"])

    # Compute delta_dG = dG_gt - dG_pred (positive = improved)
    merged["dG_gt"] = merged["code"].map(gt_dg)
    merged = merged.dropna(subset=["dG_gt"])
    merged["delta_dG"] = merged["dG_gt"] - merged["dG"]
    merged["imp"] = (merged["dG"] < merged["dG_gt"]).astype(int)

    return merged


def scatter_plot(ax, df, x_col_suffix, y_col, title, xlabel, ylabel):
    """Scatter plot with per-CDR coloring."""
    all_x, all_y = [], []

    for cdr_prefix, cdr_label in CDR_PAIRS:
        col = f"{cdr_prefix}_{x_col_suffix}"
        if col not in df.columns:
            continue
        mask = df[col].notna() & df[y_col].notna()
        x = df.loc[mask, col].values
        y = df.loc[mask, y_col].values
        if len(x) == 0:
            continue
        ax.scatter(x, y, s=12, alpha=0.5, label=cdr_label, color=COLORS[cdr_label])
        all_x.extend(x)
        all_y.extend(y)

    all_x, all_y = np.array(all_x), np.array(all_y)
    if len(all_x) > 2:
        rp, pp = stats.pearsonr(all_x, all_y)
        rs, ps = stats.spearmanr(all_x, all_y)
        ax.text(0.02, 0.98,
                f"Pearson r={rp:.3f} (p={pp:.1e})\nSpearman r={rs:.3f} (p={ps:.1e})",
                transform=ax.transAxes, fontsize=7, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=6, loc="best", markerscale=1.5)

    return all_x, all_y


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    plot_dir = os.path.join(script_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    # Load ground truth dG
    gt = pd.read_csv(GT_CSV)
    gt_dg = dict(zip(gt["PDB_Name"], gt["dG"]))

    # Print header
    print(f"{'Setting':<25} {'Metric':<16} {'Pearson r':>10} {'p':>12} {'Spearman r':>10} {'p':>12} {'N':>6}")
    print("-" * 95)

    # Combined figures
    fig_aar, axes_aar = plt.subplots(1, 3, figsize=(15, 4.5))
    fig_rmsd, axes_rmsd = plt.subplots(1, 3, figsize=(15, 4.5))
    fig_aar.suptitle("IMP (delta_dG) vs AAR — positive = improved binding", fontsize=12)
    fig_rmsd.suptitle("IMP (delta_dG) vs RMSD — positive = improved binding", fontsize=12)

    for idx, model in enumerate(MODELS):
        short = SHORT_NAMES[model]
        df = load_merged(model, script_dir, gt_dg)
        if df is None:
            print(f"  {short}: data not found, skipping")
            continue

        imp_rate = df["imp"].mean()
        n_targets = len(df)
        print(f"  {short}: {n_targets} targets, IMP rate = {imp_rate:.3f}")

        # --- delta_dG vs AAR ---
        ax = axes_aar[idx]
        all_x, all_y = scatter_plot(ax, df, "AAR", "delta_dG",
                                    f"{short} (IMP={imp_rate:.2f})", "AAR",
                                    "delta_dG (dG_gt - dG_pred)")
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)

        fig_s, ax_s = plt.subplots(figsize=(6, 5))
        scatter_plot(ax_s, df, "AAR", "delta_dG",
                     f"{short}: IMP vs AAR", "AAR", "delta_dG (dG_gt - dG_pred)")
        ax_s.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
        fig_s.tight_layout()
        fig_s.savefig(os.path.join(plot_dir, f"imp_vs_aar_{short}.png"), dpi=150)
        plt.close(fig_s)

        if len(all_x) > 2:
            rp, pp = stats.pearsonr(all_x, all_y)
            rs, ps = stats.spearmanr(all_x, all_y)
            print(f"{short:<25} {'IMP vs AAR':<16} {rp:>10.4f} {pp:>12.2e} {rs:>10.4f} {ps:>12.2e} {len(all_x):>6}")

        # --- delta_dG vs RMSD ---
        ax = axes_rmsd[idx]
        all_x, all_y = scatter_plot(ax, df, "RMSD", "delta_dG",
                                    f"{short} (IMP={imp_rate:.2f})", "RMSD",
                                    "delta_dG (dG_gt - dG_pred)")
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)

        fig_s, ax_s = plt.subplots(figsize=(6, 5))
        scatter_plot(ax_s, df, "RMSD", "delta_dG",
                     f"{short}: IMP vs RMSD", "RMSD", "delta_dG (dG_gt - dG_pred)")
        ax_s.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
        fig_s.tight_layout()
        fig_s.savefig(os.path.join(plot_dir, f"imp_vs_rmsd_{short}.png"), dpi=150)
        plt.close(fig_s)

        if len(all_x) > 2:
            rp, pp = stats.pearsonr(all_x, all_y)
            rs, ps = stats.spearmanr(all_x, all_y)
            print(f"{short:<25} {'IMP vs RMSD':<16} {rp:>10.4f} {pp:>12.2e} {rs:>10.4f} {ps:>12.2e} {len(all_x):>6}")

    fig_aar.tight_layout(rect=[0, 0, 1, 0.93])
    fig_aar.savefig(os.path.join(plot_dir, "imp_vs_aar_all.png"), dpi=150)
    plt.close(fig_aar)

    fig_rmsd.tight_layout(rect=[0, 0, 1, 0.93])
    fig_rmsd.savefig(os.path.join(plot_dir, "imp_vs_rmsd_all.png"), dpi=150)
    plt.close(fig_rmsd)

    print(f"\nPlots saved to {plot_dir}/")


if __name__ == "__main__":
    main()
