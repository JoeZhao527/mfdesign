"""
Compute TM-score for design and protenix structures against reference.

Uses tmtools (Python binding for TM-align) to align antibody CA coordinates.
Computes antibody-level TM-score (heavy + light combined), plus per-chain scores.

Reads:
  {DESIGN_BASE}/{model}/predictions/{ID}/{ID}{suffix}        (design structures)
  {PROTENIX_PDB}/{model}/predictions/{ID}/{ID}_sample_0.pdb  (protenix structures)
  data/test_entry_pdb_files/{ID}_reference.pdb                (references)
  data/test_entry.json                                        (target list)

Writes:
  r1_eval/tmscore_results/{model}_design_tmscore.csv   (3 design CSVs)
  r1_eval/tmscore_results/{model}_protenix_tmscore.csv (3 protenix CSVs)
  stdout: comparison table
"""

import os
import sys
import json
import argparse
import multiprocessing as mp
from functools import partial

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import tmtools

import warnings
warnings.filterwarnings("ignore")

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN_BASE = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_exp_data/r1_design_structures"
PROTENIX_PDB = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_protenix_pdb"
REF_DIR = os.path.join(PROJECT_DIR, "data/test_entry_pdb_files")
TEST_JSON = os.path.join(PROJECT_DIR, "data/test_entry.json")

MODELS = [
    "gen_trial_0",
    "gen_trial_1",
    "mfdesign_test_bagel_protein_boltz_base",
]


def extract_chain_ca(structure_model, chain_id):
    """Extract CA coordinates and sequence from a single chain."""
    if chain_id == '' or chain_id not in [c.id for c in structure_model]:
        return None, ""
    chain = structure_model[chain_id]
    residues = [r for r in chain.get_residues() if r.id[0] == ' ']
    coords = []
    seq_chars = []
    for r in residues:
        if 'CA' in r:
            coords.append(r['CA'].get_vector().get_array())
            seq_chars.append(seq1(r.get_resname()))
    if len(coords) == 0:
        return None, ""
    return np.array(coords), ''.join(seq_chars)


def extract_antibody_ca(pdb_path, heavy_chain_id, light_chain_id):
    """Extract CA coords and sequence for antibody (heavy + light)."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('s', pdb_path)
    except Exception:
        return None, None, None, "", "", ""
    model = structure[0]

    heavy_coords, heavy_seq = extract_chain_ca(model, heavy_chain_id)
    light_coords, light_seq = extract_chain_ca(model, light_chain_id)

    # Combined antibody
    parts_coords = [c for c in [heavy_coords, light_coords] if c is not None]
    parts_seq = [s for s in [heavy_seq, light_seq] if s]
    if not parts_coords:
        return None, None, None, "", "", ""

    ab_coords = np.concatenate(parts_coords, axis=0)
    ab_seq = ''.join(parts_seq)

    return ab_coords, heavy_coords, light_coords, ab_seq, heavy_seq, light_seq


def compute_tmscore(ref_coords, ref_seq, pred_coords, pred_seq):
    """Compute TM-score using tmtools. Returns TM-score normalized by reference length."""
    if ref_coords is None or pred_coords is None:
        return np.nan
    if len(ref_coords) == 0 or len(pred_coords) == 0:
        return np.nan
    try:
        result = tmtools.tm_align(ref_coords, pred_coords, ref_seq, pred_seq)
        # tm_norm_chain1 = normalized by length of chain1 (reference)
        return result.tm_norm_chain1
    except Exception:
        return np.nan


def compute_single(args_tuple, ref_data, use_boltz_chains=True):
    """Compute TM-scores for a single prediction."""
    pred_path, pdb_name = args_tuple

    if pdb_name not in ref_data:
        return None

    _, orig_heavy, orig_light, _ = pdb_name.split('_')

    # Chain IDs for prediction (boltz model: A=heavy, B=light)
    if use_boltz_chains:
        pred_heavy = 'A' if orig_heavy != '' else ''
        pred_light = 'B' if orig_light != '' else ''
    else:
        pred_heavy = orig_heavy
        pred_light = orig_light

    pred_ab, pred_h, pred_l, pred_ab_seq, pred_h_seq, pred_l_seq = \
        extract_antibody_ca(pred_path, pred_heavy, pred_light)

    if pred_ab is None:
        return None

    ref = ref_data[pdb_name]

    # Antibody TM-score (heavy + light combined)
    ab_tm = compute_tmscore(ref['ab_coords'], ref['ab_seq'], pred_ab, pred_ab_seq)

    # Per-chain TM-scores
    heavy_tm = compute_tmscore(ref['heavy_coords'], ref['heavy_seq'], pred_h, pred_h_seq)
    light_tm = compute_tmscore(ref['light_coords'], ref['light_seq'], pred_l, pred_l_seq)

    return {
        'code': pdb_name,
        'antibody_TMscore': ab_tm,
        'heavy_TMscore': heavy_tm,
        'light_TMscore': light_tm,
        'file_path': pred_path,
    }


def load_reference_data(test_names):
    """Pre-load reference structures for all test targets."""
    ref_data = {}
    for pdb_name in test_names:
        ref_path = os.path.join(REF_DIR, f"{pdb_name}_reference.pdb")
        if not os.path.exists(ref_path):
            continue
        _, heavy_id, light_id, _ = pdb_name.split('_')
        ab_coords, heavy_coords, light_coords, ab_seq, heavy_seq, light_seq = \
            extract_antibody_ca(ref_path, heavy_id, light_id)
        if ab_coords is None:
            continue
        ref_data[pdb_name] = {
            'ab_coords': ab_coords,
            'heavy_coords': heavy_coords,
            'light_coords': light_coords,
            'ab_seq': ab_seq,
            'heavy_seq': heavy_seq,
            'light_seq': light_seq,
        }
    return ref_data


def collect_predictions(data_dir, pdb_names_set, suffix):
    """Collect prediction PDB paths matching the target list."""
    tasks = []
    if not os.path.isdir(data_dir):
        return tasks
    for item in sorted(os.listdir(data_dir)):
        item_path = os.path.join(data_dir, item)
        if not os.path.isdir(item_path):
            continue
        if item not in pdb_names_set:
            continue
        # Find PDB with matching suffix
        for fname in os.listdir(item_path):
            if fname.endswith(suffix) and not fname.endswith('_relaxed.pdb'):
                tasks.append((os.path.join(item_path, fname), item))
                break  # one per target dir
    return tasks


def run_evaluation(tasks, ref_data, cpus, use_boltz_chains=True):
    """Run TM-score computation in parallel."""
    func = partial(compute_single, ref_data=ref_data, use_boltz_chains=use_boltz_chains)
    if cpus > 1:
        with mp.Pool(cpus) as pool:
            results = pool.map(func, tasks)
    else:
        results = [func(t) for t in tasks]
    return [r for r in results if r is not None]


def summarize(results, label):
    """Compute mean +/- SEM for TM-score metrics, grouped by code."""
    if not results:
        return None
    df = pd.DataFrame(results)
    grouped = df.groupby('code').mean(numeric_only=True)
    means = grouped.mean(numeric_only=True)
    sems = grouped.sem(numeric_only=True)
    n_targets = len(grouped)

    row = {'label': label, 'n_targets': n_targets}
    for metric in ['antibody_TMscore', 'heavy_TMscore', 'light_TMscore']:
        if metric in means.index:
            m, s = means[metric], sems[metric]
            row[metric] = f"{m:.4f} +/- {s:.4f}"
            row[f"{metric}_mean"] = m
        else:
            row[metric] = "N/A"
            row[f"{metric}_mean"] = np.nan
    return row


def main():
    parser = argparse.ArgumentParser(description="Compute TM-score for antibody structures")
    parser.add_argument("--cpus", type=int, default=40)
    parser.add_argument("--design-suffix", type=str, default="_model_0.pdb",
                        help="Suffix for design structure PDB files")
    parser.add_argument("--protenix-suffix", type=str, default="_sample_0.pdb",
                        help="Suffix for protenix structure PDB files")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, "tmscore_results")
    os.makedirs(out_dir, exist_ok=True)

    # Load test targets
    with open(TEST_JSON) as f:
        test_names = json.load(f)
    test_names_set = set(test_names)
    print(f"Loaded {len(test_names)} test targets")

    # Load reference structures
    print("Loading reference structures...")
    ref_data = load_reference_data(test_names)
    print(f"  Loaded {len(ref_data)} references")

    summary_rows = []

    for model in MODELS:
        short_name = model.replace("mfdesign_test_bagel_protein_boltz_base", "mfdesign_bagel")

        # --- Design structures ---
        design_dir = os.path.join(DESIGN_BASE, model, "predictions")
        design_tasks = collect_predictions(design_dir, test_names_set, args.design_suffix)
        if design_tasks:
            print(f"\n[{short_name}] Design: {len(design_tasks)} structures")
            design_results = run_evaluation(design_tasks, ref_data, args.cpus, use_boltz_chains=True)
            out_csv = os.path.join(out_dir, f"{model}_design_tmscore.csv")
            pd.DataFrame(design_results).to_csv(out_csv, index=False)
            print(f"  Saved {len(design_results)} results to {out_csv}")
            row = summarize(design_results, f"{short_name} (design)")
            if row:
                summary_rows.append(row)
        else:
            print(f"\n[{short_name}] Design: no structures found in {design_dir}")

        # --- Protenix structures ---
        protenix_dir = os.path.join(PROTENIX_PDB, model, "predictions")
        protenix_tasks = collect_predictions(protenix_dir, test_names_set, args.protenix_suffix)
        if protenix_tasks:
            print(f"[{short_name}] Protenix: {len(protenix_tasks)} structures")
            protenix_results = run_evaluation(protenix_tasks, ref_data, args.cpus, use_boltz_chains=True)
            out_csv = os.path.join(out_dir, f"{model}_protenix_tmscore.csv")
            pd.DataFrame(protenix_results).to_csv(out_csv, index=False)
            print(f"  Saved {len(protenix_results)} results to {out_csv}")
            row = summarize(protenix_results, f"{short_name} (protenix)")
            if row:
                summary_rows.append(row)
        else:
            print(f"[{short_name}] Protenix: no structures found in {protenix_dir}")

    # Print comparison table
    if summary_rows:
        print("\n" + "=" * 90)
        print("  TM-score Comparison (mean +/- SEM, grouped by code)")
        print("=" * 90)
        summary_df = pd.DataFrame(summary_rows)
        display_cols = ['label', 'n_targets', 'antibody_TMscore', 'heavy_TMscore', 'light_TMscore']
        available = [c for c in display_cols if c in summary_df.columns]
        print(summary_df[available].to_string(index=False))

        # Save numeric summary
        numeric_cols = ['label', 'n_targets'] + [
            f"{m}_mean" for m in ['antibody_TMscore', 'heavy_TMscore', 'light_TMscore']
            if f"{m}_mean" in summary_df.columns
        ]
        summary_csv = os.path.join(out_dir, "tmscore_comparison.csv")
        summary_df[numeric_cols].to_csv(summary_csv, index=False)
        print(f"\nSummary saved to {summary_csv}")
    else:
        print("\nNo results found!")


if __name__ == "__main__":
    main()
