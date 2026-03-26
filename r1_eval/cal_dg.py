"""
Compute dG (binding free energy) for design structures using PyRosetta.

For gen_trial_0 and gen_trial_1, reuses existing dG from binding_res/.
For mfdesign_test_bagel_protein_boltz_base, computes from relaxed PDBs.

Reads:
  {DESIGN_BASE}/{model}/predictions/{ID}/{ID}_model_0_relaxed.pdb
  evaluate/AbX_eval/imp_results/test_ground_binding.csv   (ground truth dG)
  evaluate/AbX_eval/regular_list.json / nano_list.json    (antibody classification)
  binding_res/gen_trial_{0,1}/mfdesign_uniform_imp_energy_validation.csv  (reuse)

Writes:
  r1_eval/dg_results/{model}_dg.csv   (3 per-sample dG CSVs)
  stdout: IMP values per setting
"""

import os
import sys
import json
import argparse
import multiprocessing as mp
import shutil
import traceback
import logging

import pandas as pd
import numpy as np
from Bio.PDB import PDBParser

# PyRosetta imports
from pyrosetta import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

import warnings
warnings.filterwarnings("ignore")

init('-use_input_sc -input_ab_scheme AHo_Scheme -ignore_unrecognized_res '
     '-ignore_zero_occupancy false -load_PDB_components true -relax:default_repeats 2 -no_fconfig')

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DESIGN_BASE = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_exp_data/r1_design_structures"
BINDING_RES = os.path.join(PROJECT_DIR, "binding_res")
GT_CSV = os.path.join(PROJECT_DIR, "evaluate/AbX_eval/imp_results/test_ground_binding.csv")
ANTI_JSON = os.path.join(PROJECT_DIR, "evaluate/AbX_eval/regular_list.json")
NANO_JSON = os.path.join(PROJECT_DIR, "evaluate/AbX_eval/nano_list.json")

MODELS = [
    "gen_trial_0",
    "gen_trial_1",
    "mfdesign_test_bagel_protein_boltz_base",
]

# Models that already have dG computed in binding_res/
REUSE_MODELS = {
    "gen_trial_0": os.path.join(BINDING_RES, "gen_trial_0", "mfdesign_uniform_imp_energy_validation.csv"),
    "gen_trial_1": os.path.join(BINDING_RES, "gen_trial_1", "mfdesign_uniform_imp_energy_validation.csv"),
}


def remap_chain_ids(pdb_core_name):
    """Map original chain IDs to boltz ordering (A, B, C, D, ...) for Rosetta interface."""
    parts = pdb_core_name.split('_')
    code = parts[0]
    order_name = ['A', 'B', 'C', 'D', 'E', 'F']
    sum_chain = len(''.join(parts[1:]))
    select_name = order_name[:sum_chain]
    anti_count = len(parts[1] + parts[2])

    if anti_count == 1:  # nanobody
        heavy_name = select_name[0]
        light_name = ''
        ag_name = ''.join(select_name[1:])
    elif anti_count == 2:  # full antibody
        heavy_name = select_name[0]
        light_name = select_name[1]
        ag_name = ''.join(select_name[2:])
    else:
        heavy_name = select_name[0]
        light_name = ''
        ag_name = ''.join(select_name[1:])

    new_name = f'{code}_{heavy_name}_{light_name}_{ag_name}'
    return new_name


def pyrosetta_interface_energy(pdb_path, interface):
    pose = pose_from_pdb(pdb_path)
    mover = InterfaceAnalyzerMover()
    mover.set_interface(interface)
    mover.set_scorefunction(create_score_function('ref2015'))
    mover.apply(pose)
    return pose.scores['dG_separated']


def compute_dg_single(args_tuple):
    """Compute dG for a single PDB file."""
    pdb_path, pdb_name, original_name = args_tuple
    try:
        code, heavy, light, antigen = pdb_name.split('_')
        antibody_chains = heavy + light
        antigen_chains = antigen
        interface = f"{antibody_chains}_{antigen_chains}"
        dg = pyrosetta_interface_energy(pdb_path, interface)
        return original_name, dg
    except Exception as e:
        print(f"  Error computing dG for {original_name}: {e}")
        return original_name, np.nan


def collect_pdb_tasks(pdb_dir, suffix='_model_0_relaxed.pdb'):
    """Collect PDB files and remap chain IDs, matching cal_interface_energy.py logic."""
    tasks = []
    items = [d for d in os.listdir(pdb_dir)
             if not d.endswith('.log') and not d.endswith('.csv')]

    for item_dir in sorted(items):
        item_path = os.path.join(pdb_dir, item_dir)
        if not os.path.isdir(item_path):
            continue
        for fname in os.listdir(item_path):
            if fname.endswith(suffix) and not fname.endswith('_relaxed_relaxed.pdb'):
                pdb_path = os.path.join(item_path, fname)
                # Extract original name from filename
                fname_parts = fname.split('_')
                original_name = f'{fname_parts[0]}_{fname_parts[1]}_{fname_parts[2]}_{fname_parts[3]}'
                # Remap to boltz chain ordering for Rosetta
                remapped_name = remap_chain_ids(original_name)
                tasks.append((pdb_path, remapped_name, original_name))
    return tasks


def compute_imp(dg_csv_path, gt_csv_path, anti_json_path, nano_json_path):
    """Compute IMP metric: fraction of samples with dG < ground truth dG."""
    gt = pd.read_csv(gt_csv_path)
    samples = pd.read_csv(dg_csv_path)
    samples = samples.dropna(subset=['dG'])

    with open(anti_json_path) as f:
        anti_list = json.load(f)
    with open(nano_json_path) as f:
        nano_list = json.load(f)

    gt_dg = dict(zip(gt['PDB_Name'], gt['dG']))

    results = {}
    for label, subset_list in [('anti', anti_list), ('nano', nano_list), ('all', None)]:
        if subset_list is not None:
            df = samples[samples['PDB_Name'].isin(subset_list)]
        else:
            df = samples

        succ = 0
        total = 0
        for _, row in df.iterrows():
            if row['PDB_Name'] in gt_dg:
                total += 1
                if row['dG'] < gt_dg[row['PDB_Name']]:
                    succ += 1

        if total > 0:
            p = succ / total
            n_targets = df['PDB_Name'].nunique()
            se = np.sqrt(p * (1 - p) / n_targets) if n_targets > 0 else 0
            results[label] = (p, se, total)
        else:
            results[label] = (0, 0, 0)

    return results


def main():
    parser = argparse.ArgumentParser(description="Compute dG for design structures")
    parser.add_argument("--cpus", type=int, default=20)
    parser.add_argument("--recompute-all", action="store_true",
                        help="Recompute dG for all models instead of reusing existing")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(script_dir, "dg_results")
    os.makedirs(out_dir, exist_ok=True)

    for model in MODELS:
        out_csv = os.path.join(out_dir, f"{model}_dg.csv")

        # Reuse existing dG if available
        if not args.recompute_all and model in REUSE_MODELS:
            src = REUSE_MODELS[model]
            if os.path.exists(src):
                print(f"[{model}] Reusing existing dG from {src}")
                shutil.copy2(src, out_csv)
                continue

        # Compute dG
        pdb_dir = os.path.join(DESIGN_BASE, model, "predictions")
        print(f"[{model}] Computing dG from {pdb_dir}")
        tasks = collect_pdb_tasks(pdb_dir, suffix='_model_0_relaxed.pdb')
        print(f"  Found {len(tasks)} PDB files")

        with mp.Pool(args.cpus) as pool:
            results = pool.map(compute_dg_single, tasks)

        df = pd.DataFrame(results, columns=["PDB_Name", "dG"])
        df.to_csv(out_csv, index=False)
        print(f"  Saved to {out_csv} ({len(df)} entries)")

    # Compute IMP for all 3
    print("\n" + "=" * 60)
    print("  IMP Results (dG < ground truth dG)")
    print("=" * 60)
    print(f"{'Model':<45} {'anti':>18} {'nano':>18} {'all':>18}")
    print("-" * 100)

    for model in MODELS:
        out_csv = os.path.join(out_dir, f"{model}_dg.csv")
        if not os.path.exists(out_csv):
            continue
        imp = compute_imp(out_csv, GT_CSV, ANTI_JSON, NANO_JSON)
        parts = []
        for label in ['anti', 'nano', 'all']:
            p, se, n = imp[label]
            parts.append(f"{p:.4f}±{se:.4f} (n={n})")
        short = model.replace("mfdesign_test_bagel_protein_boltz_base", "mfdesign_bagel")
        print(f"{short:<45} {parts[0]:>18} {parts[1]:>18} {parts[2]:>18}")


if __name__ == "__main__":
    main()
