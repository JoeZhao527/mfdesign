"""
Plot dG vs AAR and dG vs RMSD scatter plots for design structures.

dG comes from relaxed PDBs, RMSD/AAR from unrelaxed evaluation.
Merged on target code (PDB_Name).

Reads:
  r1_eval/dg_results/{model}_dg.csv
  {DESIGN_BASE}/{model}/predictions/results.csv

Writes:
  r1_eval/plots/dg_vs_aar_all.png
  r1_eval/plots/dg_vs_rmsd_all.png
  r1_eval/plots/dg_vs_aar_{model}.png   (3 individual)
  r1_eval/plots/dg_vs_rmsd_{model}.png  (3 individual)
  stdout: correlation coefficients
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

DESIGN_BASE = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_exp_data/r1_design_structures"

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


def load_merged(model, script_dir):
    """Load and merge dG with RMSD/AAR data."""
    dg_csv = os.path.join(script_dir, "dg_results", f"{model}_dg.csv")
    results_csv = os.path.join(DESIGN_BASE, model, "predictions", "results.csv")

    if not os.path.exists(dg_csv) or not os.path.exists(results_csv):
        return None

    dg = pd.read_csv(dg_csv)
    results = pd.read_csv(results_csv)

    merged = results.merge(dg, left_on="code", right_on="PDB_Name", how="inner")
    merged = merged.dropna(subset=["dG"])
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

    # Print header
    print(f"{'Setting':<25} {'Metric':<12} {'Pearson r':>10} {'p':>12} {'Spearman r':>10} {'p':>12} {'N':>6}")
    print("-" * 90)

    # Combined figures
    fig_aar, axes_aar = plt.subplots(1, 3, figsize=(15, 4.5))
    fig_rmsd, axes_rmsd = plt.subplots(1, 3, figsize=(15, 4.5))
    fig_aar.suptitle("dG vs AAR (per CDR, per target)", fontsize=13)
    fig_rmsd.suptitle("dG vs RMSD (per CDR, per target)", fontsize=13)

    for idx, model in enumerate(MODELS):
        short = SHORT_NAMES[model]
        df = load_merged(model, script_dir)
        if df is None:
            print(f"  {short}: data not found, skipping")
            continue

        # --- dG vs AAR ---
        ax = axes_aar[idx]
        all_x, all_y = scatter_plot(ax, df, "AAR", "dG", f"{short} design", "AAR", "dG (kcal/mol)")

        # Individual plot
        fig_s, ax_s = plt.subplots(figsize=(6, 5))
        scatter_plot(ax_s, df, "AAR", "dG", f"{short} design: dG vs AAR", "AAR", "dG (kcal/mol)")
        fig_s.tight_layout()
        fig_s.savefig(os.path.join(plot_dir, f"dg_vs_aar_{short}.png"), dpi=150)
        plt.close(fig_s)

        if len(all_x) > 2:
            rp, pp = stats.pearsonr(all_x, all_y)
            rs, ps = stats.spearmanr(all_x, all_y)
            print(f"{short:<25} {'dG vs AAR':<12} {rp:>10.4f} {pp:>12.2e} {rs:>10.4f} {ps:>12.2e} {len(all_x):>6}")

        # --- dG vs RMSD ---
        ax = axes_rmsd[idx]
        all_x, all_y = scatter_plot(ax, df, "RMSD", "dG", f"{short} design", "RMSD", "dG (kcal/mol)")

        # Individual plot
        fig_s, ax_s = plt.subplots(figsize=(6, 5))
        scatter_plot(ax_s, df, "RMSD", "dG", f"{short} design: dG vs RMSD", "RMSD", "dG (kcal/mol)")
        fig_s.tight_layout()
        fig_s.savefig(os.path.join(plot_dir, f"dg_vs_rmsd_{short}.png"), dpi=150)
        plt.close(fig_s)

        if len(all_x) > 2:
            rp, pp = stats.pearsonr(all_x, all_y)
            rs, ps = stats.spearmanr(all_x, all_y)
            print(f"{short:<25} {'dG vs RMSD':<12} {rp:>10.4f} {pp:>12.2e} {rs:>10.4f} {ps:>12.2e} {len(all_x):>6}")

    fig_aar.tight_layout(rect=[0, 0, 1, 0.93])
    fig_aar.savefig(os.path.join(plot_dir, "dg_vs_aar_all.png"), dpi=150)
    plt.close(fig_aar)

    fig_rmsd.tight_layout(rect=[0, 0, 1, 0.93])
    fig_rmsd.savefig(os.path.join(plot_dir, "dg_vs_rmsd_all.png"), dpi=150)
    plt.close(fig_rmsd)

    print(f"\nPlots saved to {plot_dir}/")


if __name__ == "__main__":
    main()
