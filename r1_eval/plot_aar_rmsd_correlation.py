"""
Plot AAR vs RMSD scatter plots for each setting, with Pearson/Spearman correlation.

Reads:
  {DESIGN_BASE}/{model}/predictions/results.csv   (3 design results)
  {PROTENIX_PDB}/{model}/predictions/results.csv   (3 protenix results)

Writes:
  r1_eval/plots/aar_vs_rmsd_{setting}.png          (6 scatter plots)
  r1_eval/plots/aar_vs_rmsd_all.png                (combined 2x3 figure)
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
PROTENIX_PDB = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_protenix_pdb"

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
    ("heavy_cdr1_AAR", "heavy_cdr1_RMSD", "H-CDR1"),
    ("heavy_cdr2_AAR", "heavy_cdr2_RMSD", "H-CDR2"),
    ("heavy_cdr3_AAR", "heavy_cdr3_RMSD", "H-CDR3"),
    ("light_cdr1_AAR", "light_cdr1_RMSD", "L-CDR1"),
    ("light_cdr2_AAR", "light_cdr2_RMSD", "L-CDR2"),
    ("light_cdr3_AAR", "light_cdr3_RMSD", "L-CDR3"),
]

COLORS = {
    "H-CDR1": "#1f77b4",
    "H-CDR2": "#ff7f0e",
    "H-CDR3": "#d62728",
    "L-CDR1": "#9467bd",
    "L-CDR2": "#8c564b",
    "L-CDR3": "#2ca02c",
}


def load_setting(csv_path):
    if not os.path.exists(csv_path):
        return None
    return pd.read_csv(csv_path)


def plot_single_setting(df, title, ax):
    """Plot all 6 CDR AAR-vs-RMSD pairs on one axis."""
    all_aar, all_rmsd = [], []

    for aar_col, rmsd_col, cdr_label in CDR_PAIRS:
        if aar_col not in df.columns or rmsd_col not in df.columns:
            continue
        mask = df[aar_col].notna() & df[rmsd_col].notna()
        x = df.loc[mask, aar_col].values
        y = df.loc[mask, rmsd_col].values
        if len(x) == 0:
            continue
        ax.scatter(x, y, s=12, alpha=0.5, label=cdr_label, color=COLORS[cdr_label])
        all_aar.extend(x)
        all_rmsd.extend(y)

    all_aar = np.array(all_aar)
    all_rmsd = np.array(all_rmsd)

    # Correlation
    corr_text = ""
    if len(all_aar) > 2:
        r_pearson, p_pearson = stats.pearsonr(all_aar, all_rmsd)
        r_spearman, p_spearman = stats.spearmanr(all_aar, all_rmsd)
        corr_text = f"Pearson r={r_pearson:.3f} (p={p_pearson:.1e})\nSpearman r={r_spearman:.3f} (p={p_spearman:.1e})"
        ax.text(0.02, 0.98, corr_text, transform=ax.transAxes,
                fontsize=7, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    ax.set_xlabel("AAR")
    ax.set_ylabel("RMSD")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=6, loc="upper right", markerscale=1.5)

    return corr_text


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    plot_dir = os.path.join(script_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    settings = []
    for model in MODELS:
        short = SHORT_NAMES[model]
        design_csv = os.path.join(DESIGN_BASE, model, "predictions", "results.csv")
        protenix_csv = os.path.join(PROTENIX_PDB, model, "predictions", "results.csv")
        settings.append((f"{short}_design", design_csv))
        settings.append((f"{short}_protenix", protenix_csv))

    # Combined figure: 2 rows (design/protenix) x 3 cols (models)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("AAR vs RMSD Correlation (per CDR, per target)", fontsize=13)

    print(f"{'Setting':<30} {'Pearson r':>12} {'p-value':>12} {'Spearman r':>12} {'p-value':>12} {'N':>6}")
    print("-" * 90)

    for idx, (label, csv_path) in enumerate(settings):
        df = load_setting(csv_path)
        if df is None:
            print(f"{label:<30} -- file not found --")
            continue

        row, col = divmod(idx, 3)
        ax = axes[row, col]
        plot_single_setting(df, label, ax)

        # Also save individual plot
        fig_single, ax_single = plt.subplots(figsize=(6, 5))
        plot_single_setting(df, label, ax_single)
        fig_single.tight_layout()
        fig_single.savefig(os.path.join(plot_dir, f"aar_vs_rmsd_{label}.png"), dpi=150)
        plt.close(fig_single)

        # Print correlation summary
        all_aar, all_rmsd = [], []
        for aar_col, rmsd_col, _ in CDR_PAIRS:
            if aar_col in df.columns and rmsd_col in df.columns:
                mask = df[aar_col].notna() & df[rmsd_col].notna()
                all_aar.extend(df.loc[mask, aar_col].values)
                all_rmsd.extend(df.loc[mask, rmsd_col].values)
        all_aar, all_rmsd = np.array(all_aar), np.array(all_rmsd)
        if len(all_aar) > 2:
            rp, pp = stats.pearsonr(all_aar, all_rmsd)
            rs, ps = stats.spearmanr(all_aar, all_rmsd)
            print(f"{label:<30} {rp:>12.4f} {pp:>12.2e} {rs:>12.4f} {ps:>12.2e} {len(all_aar):>6}")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    combined_path = os.path.join(plot_dir, "aar_vs_rmsd_all.png")
    fig.savefig(combined_path, dpi=150)
    plt.close(fig)

    print(f"\nPlots saved to {plot_dir}/")


if __name__ == "__main__":
    main()
