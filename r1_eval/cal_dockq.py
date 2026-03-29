"""
Compute DockQ for antibody-antigen complexes.

Reports per-interface DockQ (heavy-antigen, light-antigen, heavy-light)
and an antigen-interface average.

Requires: pip install DockQ

Usage:
  python cal_dockq.py \
      --pred_dir  r1_abl/inference_outputs/rabd_h3_only_trial0/predictions \
      --ref_dir   r1_abl/mfdesign_rabd_npz/reference_pdbs \
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
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

import warnings
warnings.filterwarnings("ignore")


def build_chain_map(pdb_name):
    """Build {native_chain: model_chain} mapping from pdb_name.

    pdb_name format: {code}_{heavy}_{light}_{antigen_chars}
    Boltz prediction: A=heavy, B=light, C,D,...=antigen chains
    """
    parts = pdb_name.split("_")
    heavy_id = parts[1]
    light_id = parts[2]
    antigen_ids = list(parts[3])  # e.g. "AC" -> ['A','C']

    chain_map = {heavy_id: "A", light_id: "B"}
    for i, ag_id in enumerate(antigen_ids):
        chain_map[ag_id] = chr(ord("C") + i)
    return chain_map, heavy_id, light_id, antigen_ids


def compute_single(args_tuple):
    """Compute DockQ for a single prediction."""
    pred_path, pdb_name, ref_dir = args_tuple

    ref_path = os.path.join(ref_dir, f"{pdb_name}_reference.pdb")
    if not os.path.exists(ref_path):
        return None

    try:
        chain_map, heavy_id, light_id, antigen_ids = build_chain_map(pdb_name)
    except Exception:
        return None

    try:
        model = load_PDB(pred_path)
        native = load_PDB(ref_path)
        result_dict, _ = run_on_all_native_interfaces(
            model, native, chain_map=chain_map)
    except Exception as e:
        print(f"  [WARN] {pdb_name}: {e}")
        return None

    row = {"code": pdb_name, "file_path": pred_path}

    # Collect per-interface DockQ
    antigen_dockqs = []
    for interface_key, metrics in result_dict.items():
        dq = metrics["DockQ"]
        fnat = metrics["fnat"]
        irmsd = metrics["iRMSD"]
        lrmsd = metrics["LRMSD"]

        # Determine interface type from native chain IDs
        c1, c2 = list(interface_key)
        is_antigen_interface = (c1 in antigen_ids or c2 in antigen_ids)
        is_heavy_antigen = (heavy_id in interface_key and
                            any(a in interface_key for a in antigen_ids))

        row[f"DockQ_{interface_key}"] = dq
        row[f"Fnat_{interface_key}"] = fnat
        row[f"iRMSD_{interface_key}"] = irmsd
        row[f"LRMSD_{interface_key}"] = lrmsd

        if is_antigen_interface:
            antigen_dockqs.append(dq)

        if is_heavy_antigen:
            row["DockQ_heavy_antigen"] = dq

    # Average DockQ over antigen-containing interfaces
    if antigen_dockqs:
        row["DockQ_antigen_avg"] = float(np.mean(antigen_dockqs))
    else:
        row["DockQ_antigen_avg"] = np.nan

    return row


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


def main():
    ap = argparse.ArgumentParser(description="DockQ evaluation")
    ap.add_argument("--pred_dir", required=True)
    ap.add_argument("--ref_dir", required=True)
    ap.add_argument("--target_json", required=True)
    ap.add_argument("--suffix", default="_model_0.pdb")
    ap.add_argument("--cpus", type=int, default=20)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    with open(args.target_json) as f:
        target_names = json.load(f)
    target_set = set(target_names)
    print(f"Targets: {len(target_names)}")

    tasks_raw = collect_predictions(args.pred_dir, target_set, args.suffix)
    # Attach ref_dir to each task for the worker
    tasks = [(path, name, args.ref_dir) for path, name in tasks_raw]
    print(f"Predictions: {len(tasks)}")
    if not tasks:
        return

    if args.cpus > 1:
        with mp.Pool(args.cpus) as pool:
            results = pool.map(compute_single, tasks)
    else:
        results = [compute_single(t) for t in tasks]

    results = [r for r in results if r is not None]
    print(f"Valid results: {len(results)}")
    if not results:
        return

    df = pd.DataFrame(results)

    # Summary
    summary_cols = ["DockQ_antigen_avg"]
    if "DockQ_heavy_antigen" in df.columns:
        summary_cols.append("DockQ_heavy_antigen")
    available = [c for c in summary_cols if c in df.columns]

    grouped = df.groupby("code")[available].mean()
    means = grouped.mean()
    sems = grouped.sem()

    print(f"\n{'='*60}")
    print(f"  DockQ Results (mean ± SEM, n_targets={len(grouped)})")
    print(f"{'='*60}")
    for m in available:
        print(f"  {m:30s}: {means[m]:.4f} ± {sems[m]:.4f}")

    out_path = args.out or os.path.join(args.pred_dir, "dockq.csv")
    df.to_csv(out_path, index=False)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
