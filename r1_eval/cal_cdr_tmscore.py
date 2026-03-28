"""
Compute per-CDR and whole-chain TM-scores for antibody design evaluation.

Reads prediction PDBs + reference PDBs + YAML CDR masks, computes TM-score at:
  - Whole antibody (heavy + light)
  - Per chain (heavy, light)
  - Per CDR (H1, H2, H3, L1, L2, L3)
  - CDR3 loop (H3 core excluding flanking residues)

Usage:
  python cal_cdr_tmscore.py \
      --pred_dir  inference_outputs/rabd_h3_only_trial0/predictions \
      --ref_dir   mfdesign_rabd_npz/reference_pdbs \
      --yaml_dir  mfdesign_rabd_npz/test_yaml_dir/h3_only \
      --target_json mfdesign_rabd_npz/rabd_split.json \
      --suffix _model_0.pdb \
      --cpus 40 \
      --out results_tmscore.csv
"""

import os
import json
import argparse
import multiprocessing as mp
from functools import partial

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from ruamel.yaml import YAML
import tmtools

import warnings
warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# CDR numbering from YAML spec_mask  (mirrors boltz_cdr_numbering in abx)
# ---------------------------------------------------------------------------

def read_yaml_masks(yaml_path, heavy_id, light_id):
    """Read spec_mask strings from a YAML file."""
    yaml = YAML()
    with open(yaml_path) as f:
        data = yaml.load(f)
    heavy_mask = light_mask = ""
    for entry in data["sequences"]:
        p = entry["protein"]
        if p["id"] == heavy_id:
            heavy_mask = p["spec_mask"]
        elif p["id"] == light_id:
            light_mask = p["spec_mask"]
    return heavy_mask, light_mask


def mask_to_region_def(mask_str, chain_type):
    """Convert a binary mask string to region IDs.

    Region IDs follow the convention in abx:
      Heavy: FR1=0, CDR1=1, FR2=2, CDR2=3, FR3=4, CDR3=5, FR4=6
      Light: offset by 7  (FR1=7, CDR1=8, ..., FR4=13)
    """
    n = len(mask_str)
    region_def = np.full(n, -1, dtype=np.int32)
    last = None
    idx = -1
    for i, ch in enumerate(mask_str):
        if ch != last:
            idx += 1
        region_def[i] = idx
        last = ch
    if chain_type == "L":
        region_def += 7
    return region_def


# ---------------------------------------------------------------------------
# PDB helpers
# ---------------------------------------------------------------------------

def extract_chain_ca(model, chain_id):
    """Extract CA coords (Nx3) and 1-letter sequence from a single chain."""
    if chain_id not in [c.id for c in model]:
        return None, ""
    chain = model[chain_id]
    coords, seq_chars = [], []
    for r in chain.get_residues():
        if r.id[0] != " ":
            continue
        if "CA" in r:
            coords.append(r["CA"].get_vector().get_array())
            seq_chars.append(seq1(r.get_resname()))
    if not coords:
        return None, ""
    return np.array(coords, dtype=np.float64), "".join(seq_chars)


# ---------------------------------------------------------------------------
# TM-score helpers
# ---------------------------------------------------------------------------

def _tm_score(ref_coords, ref_seq, pred_coords, pred_seq):
    """Compute TM-score normalised by reference length. Returns NaN on failure."""
    if ref_coords is None or pred_coords is None:
        return np.nan
    if len(ref_coords) == 0 or len(pred_coords) == 0:
        return np.nan
    if len(ref_coords) < 3 or len(pred_coords) < 3:
        return np.nan
    try:
        res = tmtools.tm_align(ref_coords, pred_coords, ref_seq, pred_seq)
        return res.tm_norm_chain1
    except Exception:
        return np.nan


# CDR region ID -> label mapping
CDR_REGIONS = {
    1: "heavy_cdr1",
    3: "heavy_cdr2",
    5: "heavy_cdr3",
    8: "light_cdr1",
    10: "light_cdr2",
    12: "light_cdr3",
}


def compute_cdr_tmscores(ref_coords, ref_seq, pred_coords, pred_seq, cdr_def):
    """Compute TM-score for each CDR region defined in cdr_def."""
    results = {}
    for region_id, label in CDR_REGIONS.items():
        mask = cdr_def == region_id
        if mask.sum() < 3:
            results[f"{label}_TMscore"] = np.nan
            continue
        rc = ref_coords[mask]
        pc = pred_coords[mask]
        rs = "".join(c for c, m in zip(ref_seq, mask) if m)
        ps = "".join(c for c, m in zip(pred_seq, mask) if m)
        results[f"{label}_TMscore"] = _tm_score(rc, rs, pc, ps)

    # CDR3 loop (heavy_cdr3 core: skip first 4 and last 2, IMGT convention)
    h3_mask = cdr_def == 5
    h3_indices = np.where(h3_mask)[0]
    if len(h3_indices) > 6:
        loop_idx = h3_indices[4:-2]
        rc = ref_coords[loop_idx]
        pc = pred_coords[loop_idx]
        rs = "".join(ref_seq[i] for i in loop_idx)
        ps = "".join(pred_seq[i] for i in loop_idx)
        results["heavy_cdr3_loop_TMscore"] = _tm_score(rc, rs, pc, ps)
    else:
        results["heavy_cdr3_loop_TMscore"] = np.nan

    return results


# ---------------------------------------------------------------------------
# Single-target evaluation
# ---------------------------------------------------------------------------

def compute_single(args_tuple, ref_cache, yaml_dir, use_boltz_chains=True):
    """Evaluate one prediction: whole-chain + per-CDR TM-scores."""
    pred_path, pdb_name = args_tuple

    if pdb_name not in ref_cache:
        return None

    ref = ref_cache[pdb_name]
    parts = pdb_name.split("_")
    orig_heavy, orig_light = parts[1], parts[2]

    # --- prediction chain IDs ---
    if use_boltz_chains:
        pred_heavy = "A"
        pred_light = "B"
    else:
        pred_heavy = orig_heavy
        pred_light = orig_light

    parser = PDBParser(QUIET=True)
    try:
        model = parser.get_structure("s", pred_path)[0]
    except Exception:
        return None

    pred_h_coords, pred_h_seq = extract_chain_ca(model, pred_heavy)
    pred_l_coords, pred_l_seq = extract_chain_ca(model, pred_light)

    if pred_h_coords is None and pred_l_coords is None:
        return None

    row = {"code": pdb_name, "file_path": pred_path}

    # --- whole antibody TM-score ---
    parts_c = [c for c in [pred_h_coords, pred_l_coords] if c is not None]
    parts_s = [s for s in [pred_h_seq, pred_l_seq] if s]
    pred_ab_coords = np.concatenate(parts_c)
    pred_ab_seq = "".join(parts_s)
    row["antibody_TMscore"] = _tm_score(
        ref["ab_coords"], ref["ab_seq"], pred_ab_coords, pred_ab_seq
    )

    # --- per-chain TM-scores ---
    row["heavy_TMscore"] = _tm_score(
        ref["heavy_coords"], ref["heavy_seq"], pred_h_coords, pred_h_seq
    )
    row["light_TMscore"] = _tm_score(
        ref["light_coords"], ref["light_seq"], pred_l_coords, pred_l_seq
    )

    # --- per-CDR TM-scores ---
    # Build full cdr_def for prediction (same region structure as reference)
    # Lengths must match; if not, skip CDR-level
    ref_cdr_def = ref["cdr_def"]
    pred_full_coords = np.concatenate(
        [c for c in [pred_h_coords, pred_l_coords] if c is not None]
    )
    pred_full_seq = "".join([s for s in [pred_h_seq, pred_l_seq] if s])

    if len(ref_cdr_def) == len(ref["ab_coords"]) == len(pred_full_coords):
        cdr_scores = compute_cdr_tmscores(
            ref["ab_coords"], ref["ab_seq"],
            pred_full_coords, pred_full_seq,
            ref_cdr_def,
        )
        row.update(cdr_scores)
    else:
        for label in CDR_REGIONS.values():
            row[f"{label}_TMscore"] = np.nan
        row["heavy_cdr3_loop_TMscore"] = np.nan

    return row


# ---------------------------------------------------------------------------
# Reference loading
# ---------------------------------------------------------------------------

def load_references(target_names, ref_dir, yaml_dir):
    """Pre-load reference structures and CDR definitions."""
    parser = PDBParser(QUIET=True)
    cache = {}
    for name in target_names:
        ref_path = os.path.join(ref_dir, f"{name}_reference.pdb")
        yaml_path = os.path.join(yaml_dir, f"{name}.yaml")
        if not os.path.exists(ref_path) or not os.path.exists(yaml_path):
            continue

        parts = name.split("_")
        heavy_id, light_id = parts[1], parts[2]

        try:
            model = parser.get_structure("r", ref_path)[0]
        except Exception:
            continue

        h_coords, h_seq = extract_chain_ca(model, heavy_id)
        l_coords, l_seq = extract_chain_ca(model, light_id)
        if h_coords is None and l_coords is None:
            continue

        # CDR mask from YAML
        h_mask, l_mask = read_yaml_masks(yaml_path, heavy_id, light_id)
        if len(h_mask) != len(h_seq) or len(l_mask) != len(l_seq):
            continue
        h_region = mask_to_region_def(h_mask, "H")
        l_region = mask_to_region_def(l_mask, "L")
        cdr_def = np.concatenate([h_region, l_region])

        ab_parts_c = [c for c in [h_coords, l_coords] if c is not None]
        ab_parts_s = [s for s in [h_seq, l_seq] if s]

        cache[name] = {
            "ab_coords": np.concatenate(ab_parts_c),
            "ab_seq": "".join(ab_parts_s),
            "heavy_coords": h_coords,
            "heavy_seq": h_seq,
            "light_coords": l_coords,
            "light_seq": l_seq,
            "cdr_def": cdr_def,
        }
    return cache


def collect_predictions(pred_dir, target_set, suffix):
    """Collect prediction PDB paths."""
    tasks = []
    if not os.path.isdir(pred_dir):
        return tasks
    for item in sorted(os.listdir(pred_dir)):
        item_path = os.path.join(pred_dir, item)
        if not os.path.isdir(item_path):
            continue
        if item not in target_set:
            continue
        for fname in os.listdir(item_path):
            if fname.endswith(suffix) and "_relaxed" not in fname:
                tasks.append((os.path.join(item_path, fname), item))
                break
    return tasks


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Per-CDR and whole-chain TM-score evaluation"
    )
    parser.add_argument("--pred_dir", required=True,
                        help="Directory with predictions/{ID}/{ID}{suffix}")
    parser.add_argument("--ref_dir", required=True,
                        help="Directory with {ID}_reference.pdb files")
    parser.add_argument("--yaml_dir", required=True,
                        help="Directory with {ID}.yaml CDR mask files")
    parser.add_argument("--target_json", required=True,
                        help="JSON file with list of target names")
    parser.add_argument("--suffix", default="_model_0.pdb",
                        help="PDB filename suffix to match")
    parser.add_argument("--boltz_chains", action="store_true", default=True,
                        help="Predictions use boltz convention A=heavy B=light (default)")
    parser.add_argument("--original_chains", action="store_true", default=False,
                        help="Predictions use original chain IDs from filename")
    parser.add_argument("--cpus", type=int, default=40)
    parser.add_argument("--out", default=None,
                        help="Output CSV path (default: {pred_dir}/tmscore_cdr.csv)")
    args = parser.parse_args()

    use_boltz = not args.original_chains

    # Load targets
    with open(args.target_json) as f:
        target_names = json.load(f)
    target_set = set(target_names)
    print(f"Targets: {len(target_names)}")

    # Load references + CDR defs
    print("Loading references ...")
    ref_cache = load_references(target_names, args.ref_dir, args.yaml_dir)
    print(f"  Loaded {len(ref_cache)} references with CDR masks")

    # Collect predictions
    tasks = collect_predictions(args.pred_dir, target_set, args.suffix)
    print(f"Predictions: {len(tasks)}")
    if not tasks:
        print("No predictions found!")
        return

    # Run
    func = partial(compute_single, ref_cache=ref_cache,
                   yaml_dir=args.yaml_dir, use_boltz_chains=use_boltz)
    if args.cpus > 1:
        with mp.Pool(args.cpus) as pool:
            results = pool.map(func, tasks)
    else:
        results = [func(t) for t in tasks]

    results = [r for r in results if r is not None]
    print(f"Valid results: {len(results)}")

    if not results:
        print("No valid results!")
        return

    df = pd.DataFrame(results)

    # --- Summary: mean ± SEM grouped by code ---
    metric_cols = [c for c in df.columns if c.endswith("_TMscore")]
    grouped = df.groupby("code")[metric_cols].mean()
    means = grouped.mean()
    sems = grouped.sem()

    print(f"\n{'='*70}")
    print(f"  TM-score Results (mean ± SEM, n_targets={len(grouped)})")
    print(f"{'='*70}")
    for m in metric_cols:
        print(f"  {m:30s}: {means[m]:.4f} ± {sems[m]:.4f}")

    # Save
    out_path = args.out or os.path.join(args.pred_dir, "tmscore_cdr.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
