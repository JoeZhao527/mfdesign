"""
Aggregate results from all 6 evaluation runs into a comparison table.

Reads:
  {DESIGN_BASE}/{model}/predictions/results.csv   (3 design results)
  {PROTENIX_PDB}/{model}/predictions/results.csv   (3 protenix results)
  r1_eval/tmscore_results/{model}_{type}_tmscore.csv (TM-score results, optional)

Writes:
  r1_eval/r1_comparison.csv
  stdout (formatted table)
"""

import os
import pandas as pd
import numpy as np

DESIGN_BASE = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_exp_data/r1_design_structures"
PROTENIX_PDB = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_protenix_pdb"

MODELS = [
    "gen_trial_0",
    "gen_trial_1",
    "mfdesign_test_bagel_protein_boltz_base",
]

KEY_METRICS = [
    "heavy_cdr1_RMSD", "heavy_cdr1_AAR",
    "heavy_cdr2_RMSD", "heavy_cdr2_AAR",
    "heavy_cdr3_RMSD", "heavy_cdr3_AAR",
    "heavy_cdr3_Loop_RMSD", "heavy_cdr3_Loop_AAR",
    "light_cdr1_RMSD", "light_cdr1_AAR",
    "light_cdr2_RMSD", "light_cdr2_AAR",
    "light_cdr3_RMSD", "light_cdr3_AAR",
    "full_RMSD",
]

TMSCORE_METRICS = [
    "antibody_TMscore",
    "heavy_TMscore",
    "light_TMscore",
]


def load_and_summarize(csv_path, label):
    """Load results.csv, compute mean +/- SEM grouped by code."""
    if not os.path.exists(csv_path):
        print(f"  Warning: {csv_path} not found, skipping")
        return None

    df = pd.read_csv(csv_path)
    n_total = len(df)

    # Group by code (target), take mean per target, then mean across targets
    grouped = df.groupby("code").mean(numeric_only=True)
    means = grouped.mean(numeric_only=True)
    sems = grouped.sem(numeric_only=True)
    n_targets = len(grouped)

    row = {"label": label, "n_samples": n_total, "n_targets": n_targets}
    for metric in KEY_METRICS:
        if metric in means.index:
            row[metric] = f"{means[metric]:.4f} +/- {sems[metric]:.4f}"
            row[f"{metric}_mean"] = means[metric]
        else:
            row[metric] = "N/A"
            row[f"{metric}_mean"] = np.nan

    return row


def load_tmscore(tmscore_csv, label):
    """Load TM-score CSV and compute mean +/- SEM grouped by code."""
    if not os.path.exists(tmscore_csv):
        return {}

    df = pd.read_csv(tmscore_csv)
    grouped = df.groupby("code").mean(numeric_only=True)
    means = grouped.mean(numeric_only=True)
    sems = grouped.sem(numeric_only=True)

    result = {}
    for metric in TMSCORE_METRICS:
        if metric in means.index:
            result[metric] = f"{means[metric]:.4f} +/- {sems[metric]:.4f}"
            result[f"{metric}_mean"] = means[metric]
        else:
            result[metric] = "N/A"
            result[f"{metric}_mean"] = np.nan
    return result


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    tmscore_dir = os.path.join(script_dir, "tmscore_results")

    rows = []
    for model in MODELS:
        short_name = model.replace("mfdesign_test_bagel_protein_boltz_base", "mfdesign_bagel")

        # Design results
        design_csv = os.path.join(DESIGN_BASE, model, "predictions", "results.csv")
        print(f"Loading design: {model}")
        row = load_and_summarize(design_csv, f"{short_name} (design)")
        if row:
            # Merge TM-score if available
            tm_csv = os.path.join(tmscore_dir, f"{model}_design_tmscore.csv")
            row.update(load_tmscore(tm_csv, f"{short_name} (design)"))
            rows.append(row)

        # Protenix results
        protenix_csv = os.path.join(PROTENIX_PDB, model, "predictions", "results.csv")
        print(f"Loading protenix: {model}")
        row = load_and_summarize(protenix_csv, f"{short_name} (protenix)")
        if row:
            # Merge TM-score if available
            tm_csv = os.path.join(tmscore_dir, f"{model}_protenix_tmscore.csv")
            row.update(load_tmscore(tm_csv, f"{short_name} (protenix)"))
            rows.append(row)

    if not rows:
        print("No results found!")
        return

    # Build comparison table
    summary_df = pd.DataFrame(rows)

    # Print formatted table
    print("\n" + "=" * 80)
    print("  R1 Evaluation Comparison (mean +/- SEM, grouped by code)")
    print("=" * 80)

    display_cols = ["label", "n_targets"] + KEY_METRICS + TMSCORE_METRICS
    available_cols = [c for c in display_cols if c in summary_df.columns]
    print(summary_df[available_cols].to_string(index=False))

    # Save CSV with numeric means for easy downstream use
    all_metrics = KEY_METRICS + TMSCORE_METRICS
    numeric_cols = ["label", "n_samples", "n_targets"] + [
        f"{m}_mean" for m in all_metrics if f"{m}_mean" in summary_df.columns
    ]
    out_path = os.path.join(script_dir, "r1_comparison.csv")
    summary_df[numeric_cols].to_csv(out_path, index=False)
    print(f"\nNumeric results saved to {out_path}")

    # Also save the full table with +/- format
    full_path = os.path.join(script_dir, "r1_comparison_full.csv")
    summary_df[available_cols].to_csv(full_path, index=False)
    print(f"Full results saved to {full_path}")


if __name__ == "__main__":
    main()
