"""
Compute lDDT-Cα for antibody design evaluation.

Reports three granularities:
  - Whole antibody (heavy + light chains)
  - All CDR regions combined (H1+H2+H3+L1+L2+L3)
  - CDR-H3 only

lDDT is superposition-free: it compares pairwise CA distance patterns in
reference vs prediction within an inclusion radius (default 15 Å), averaged
over four thresholds (0.5, 1, 2, 4 Å).

Usage:
  python cal_lddt.py \
      --pred_dir  r1_abl/inference_outputs/rabd_h3_only_trial0/predictions \
      --ref_dir   r1_abl/mfdesign_rabd_npz/reference_pdbs \
      --yaml_dir  r1_abl/mfdesign_rabd_npz/test_yaml_dir/h3_only \
      --target_json r1_abl/mfdesign_rabd_npz/rabd_split.json \
      --suffix _model_0.pdb \
      --cpus 20
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

import warnings
warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# lDDT-Cα implementation
# ---------------------------------------------------------------------------

LDDT_THRESHOLDS = [0.5, 1.0, 2.0, 4.0]
INCLUSION_RADIUS = 15.0


def compute_lddt(ref_coords, pred_coords,
                 inclusion_radius=INCLUSION_RADIUS,
                 thresholds=LDDT_THRESHOLDS):
    """Compute lDDT-Cα.

    Args:
        ref_coords:  (N, 3) reference CA coordinates
        pred_coords: (N, 3) predicted CA coordinates

    Returns:
        global_lddt: scalar, mean over all residues
        per_residue: (N,) per-residue lDDT scores
    """
    N = len(ref_coords)
    if N == 0:
        return np.nan, np.array([])

    # Pairwise distance matrices  (N, N)
    ref_dists = np.linalg.norm(
        ref_coords[:, None, :] - ref_coords[None, :, :], axis=-1)
    pred_dists = np.linalg.norm(
        pred_coords[:, None, :] - pred_coords[None, :, :], axis=-1)

    diff = np.abs(ref_dists - pred_dists)

    # Inclusion mask: within radius in reference, exclude self-pairs
    mask = ref_dists < inclusion_radius
    np.fill_diagonal(mask, False)

    n_contacts = mask.sum(axis=1).astype(np.float64)
    n_contacts = np.clip(n_contacts, 1, None)  # avoid division by zero

    per_residue = np.zeros(N, dtype=np.float64)
    for t in thresholds:
        preserved = (diff < t) & mask
        per_residue += preserved.sum(axis=1) / n_contacts
    per_residue /= len(thresholds)

    return float(np.mean(per_residue)), per_residue


# ---------------------------------------------------------------------------
# CDR mask helpers (same as cal_cdr_tmscore.py)
# ---------------------------------------------------------------------------

CDR_REGION_IDS = {1, 3, 5, 8, 10, 12}        # all CDRs
H3_REGION_ID = 5                               # heavy CDR3


def read_yaml_masks(yaml_path, heavy_id, light_id):
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
    if chain_id not in [c.id for c in model]:
        return None, ""
    coords, seq_chars = [], []
    for r in model[chain_id].get_residues():
        if r.id[0] != " ":
            continue
        if "CA" in r:
            coords.append(r["CA"].get_vector().get_array())
            seq_chars.append(seq1(r.get_resname()))
    if not coords:
        return None, ""
    return np.array(coords, dtype=np.float64), "".join(seq_chars)


# ---------------------------------------------------------------------------
# Single-target evaluation
# ---------------------------------------------------------------------------

def compute_single(args_tuple, ref_cache, use_boltz_chains=True):
    pred_path, pdb_name = args_tuple
    if pdb_name not in ref_cache:
        return None

    ref = ref_cache[pdb_name]
    parts = pdb_name.split("_")
    orig_heavy, orig_light = parts[1], parts[2]

    pred_heavy = "A" if use_boltz_chains else orig_heavy
    pred_light = "B" if use_boltz_chains else orig_light

    parser = PDBParser(QUIET=True)
    try:
        model = parser.get_structure("s", pred_path)[0]
    except Exception:
        return None

    ph_coords, _ = extract_chain_ca(model, pred_heavy)
    pl_coords, _ = extract_chain_ca(model, pred_light)
    if ph_coords is None and pl_coords is None:
        return None

    pred_ab = np.concatenate([c for c in [ph_coords, pl_coords] if c is not None])
    ref_ab = ref["ab_coords"]
    cdr_def = ref["cdr_def"]

    if len(pred_ab) != len(ref_ab):
        return None

    # 1. Whole antibody lDDT
    ab_lddt, per_res = compute_lddt(ref_ab, pred_ab)

    # 2. Heavy chain only lDDT
    n_heavy = len(ph_coords) if ph_coords is not None else 0
    if n_heavy > 0:
        heavy_lddt = float(np.mean(per_res[:n_heavy]))
    else:
        heavy_lddt = np.nan

    # 3. All-CDR lDDT
    cdr_mask = np.isin(cdr_def, list(CDR_REGION_IDS))
    if cdr_mask.sum() > 0:
        all_cdr_lddt = float(np.mean(per_res[cdr_mask]))
    else:
        all_cdr_lddt = np.nan

    # 4. CDR-H3 lDDT
    h3_mask = cdr_def == H3_REGION_ID
    if h3_mask.sum() > 0:
        h3_lddt = float(np.mean(per_res[h3_mask]))
    else:
        h3_lddt = np.nan

    return {
        "code": pdb_name,
        "antibody_lDDT": ab_lddt,
        "heavy_lDDT": heavy_lddt,
        "all_CDR_lDDT": all_cdr_lddt,
        "heavy_cdr3_lDDT": h3_lddt,
        "file_path": pred_path,
    }


# ---------------------------------------------------------------------------
# Reference loading
# ---------------------------------------------------------------------------

def load_references(target_names, ref_dir, yaml_dir):
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

        h_mask, l_mask = read_yaml_masks(yaml_path, heavy_id, light_id)
        if len(h_mask) != len(h_seq) or len(l_mask) != len(l_seq):
            continue
        h_region = mask_to_region_def(h_mask, "H")
        l_region = mask_to_region_def(l_mask, "L")

        cache[name] = {
            "ab_coords": np.concatenate([c for c in [h_coords, l_coords] if c is not None]),
            "cdr_def": np.concatenate([h_region, l_region]),
        }
    return cache


def collect_predictions(pred_dir, target_set, suffix):
    tasks = []
    if not os.path.isdir(pred_dir):
        return tasks
    for item in sorted(os.listdir(pred_dir)):
        item_path = os.path.join(pred_dir, item)
        if not os.path.isdir(item_path) or item not in target_set:
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
    ap = argparse.ArgumentParser(description="lDDT-Cα evaluation")
    ap.add_argument("--pred_dir", required=True)
    ap.add_argument("--ref_dir", required=True)
    ap.add_argument("--yaml_dir", required=True)
    ap.add_argument("--target_json", required=True)
    ap.add_argument("--suffix", default="_model_0.pdb")
    ap.add_argument("--boltz_chains", action="store_true", default=True)
    ap.add_argument("--original_chains", action="store_true", default=False)
    ap.add_argument("--cpus", type=int, default=40)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    use_boltz = not args.original_chains

    with open(args.target_json) as f:
        target_names = json.load(f)
    target_set = set(target_names)
    print(f"Targets: {len(target_names)}")

    print("Loading references ...")
    ref_cache = load_references(target_names, args.ref_dir, args.yaml_dir)
    print(f"  Loaded {len(ref_cache)} references")

    tasks = collect_predictions(args.pred_dir, target_set, args.suffix)
    print(f"Predictions: {len(tasks)}")
    if not tasks:
        return

    func = partial(compute_single, ref_cache=ref_cache, use_boltz_chains=use_boltz)
    if args.cpus > 1:
        with mp.Pool(args.cpus) as pool:
            results = pool.map(func, tasks)
    else:
        results = [func(t) for t in tasks]

    results = [r for r in results if r is not None]
    print(f"Valid results: {len(results)}")
    if not results:
        return

    df = pd.DataFrame(results)
    metric_cols = [c for c in df.columns if "lDDT" in c]
    grouped = df.groupby("code")[metric_cols].mean()
    means = grouped.mean()
    sems = grouped.sem()

    print(f"\n{'='*60}")
    print(f"  lDDT Results (mean ± SEM, n_targets={len(grouped)})")
    print(f"{'='*60}")
    for m in metric_cols:
        print(f"  {m:25s}: {means[m]:.4f} ± {sems[m]:.4f}")

    out_path = args.out or os.path.join(args.pred_dir, "lddt.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
