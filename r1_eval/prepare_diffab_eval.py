#!/usr/bin/env python3
"""
Reorganise DiffAb (ABX) output into the directory layout that
evaluate/AbX_eval/eval_metric.py expects (--model boltz).

eval_metric.py --model boltz assumes chain A = heavy, B = light.
DiffAb keeps original PDB chain IDs, so this script remaps them.

DiffAb output layout
---------------------
  {diffab_dir}/
    0000/                          # sample 0
      {pdb_name}.pdb               # unrelaxed design  (original chain IDs)
      {pdb_name}_relaxed.pdb       # relaxed design     (may be absent)
    0001/                          # sample 1  (optional)
      ...
    reference/
      {pdb_name}.pdb               # reference structure (original chain IDs)

Target layout  (--model boltz)
------------------------------
  {output_dir}/predictions/
    {pdb_name}/
      {pdb_name}_model_0.pdb       # chains remapped: heavy→A, light→B, ag→C,D,...
      {pdb_name}_model_0_relaxed.pdb
      ...
  {output_dir}/references/
    {pdb_name}_reference.pdb       # original chain IDs (eval_metric uses them directly)

Usage
-----
  python r1_eval/prepare_diffab_eval.py \\
      --diffab_dir /hai/scratch/fangwu97/abx/output/DiffAb_design/design \\
      --output_dir /hai/scratch/fangwu97/abx/output/DiffAb_design/eval_ready
"""

import argparse
import json
import os
import re
import shutil


def find_sample_dirs(diffab_dir):
    """Return sorted list of numeric sample directories (0000, 0001, …)."""
    dirs = []
    for name in os.listdir(diffab_dir):
        full = os.path.join(diffab_dir, name)
        if os.path.isdir(full) and re.match(r"^\d{4}$", name):
            dirs.append(name)
    return sorted(dirs)


def extract_pdb_names(sample_dir):
    """Return set of pdb_names (without .pdb) from non-relaxed PDB files."""
    names = set()
    for f in os.listdir(sample_dir):
        if f.endswith(".pdb") and not f.endswith("_relaxed.pdb"):
            names.add(f[:-4])
    return names


def build_chain_mapping(pdb_name):
    """Build chain ID mapping: original → boltz convention.

    boltz convention: heavy=A, light=B, antigen chains=C,D,E,...
    For nanobodies (no light chain): heavy=A, antigen=B,C,...
    """
    parts = pdb_name.split("_")
    heavy = parts[1]       # e.g. 'H', 'B', 'A'
    light = parts[2]       # e.g. 'L', 'C', '' (empty for nanobodies)
    antigen = parts[3]     # e.g. 'A', 'AB', 'EI', 'BD'

    mapping = {}
    next_target = ord("A")

    if heavy:
        mapping[heavy] = chr(next_target)
        next_target += 1

    if light:
        mapping[light] = chr(next_target)
        next_target += 1

    for ag_char in antigen:
        if ag_char not in mapping:
            mapping[ag_char] = chr(next_target)
            next_target += 1

    return mapping


def remap_pdb_chains(src_path, dst_path, chain_mapping):
    """Rewrite PDB file with remapped chain IDs (column 22, 0-indexed pos 21)."""
    with open(src_path) as f:
        lines = f.readlines()

    with open(dst_path, "w") as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM", "TER", "ANISOU")):
                if len(line) > 21:
                    old_chain = line[21]
                    if old_chain in chain_mapping:
                        line = line[:21] + chain_mapping[old_chain] + line[22:]
            f.write(line)


def link_or_copy(src, dst, use_symlinks=True):
    """Create a symlink or copy, overwriting dst if it exists."""
    src = os.path.abspath(src)
    if os.path.exists(dst) or os.path.islink(dst):
        os.remove(dst)
    if use_symlinks:
        os.symlink(src, dst)
    else:
        shutil.copy2(src, dst)


def main():
    parser = argparse.ArgumentParser(
        description="Convert DiffAb output to eval_metric.py boltz format"
    )
    parser.add_argument(
        "--diffab_dir", type=str, required=True,
        help="DiffAb design directory containing 0000/, reference/, etc.",
    )
    parser.add_argument(
        "--output_dir", type=str, required=True,
        help="Output directory for reorganised files",
    )
    parser.add_argument(
        "--ref_dir", type=str, default=None,
        help="External reference dir (e.g. data/test_entry_pdb_files). "
             "If not set, uses {diffab_dir}/reference/.",
    )
    parser.add_argument(
        "--test_json", type=str, default=None,
        help="test_entry.json to filter targets. If not set, all targets included.",
    )
    args = parser.parse_args()

    diffab_dir = args.diffab_dir
    pred_out = os.path.join(args.output_dir, "predictions")
    ref_out = os.path.join(args.output_dir, "references")
    os.makedirs(pred_out, exist_ok=True)
    os.makedirs(ref_out, exist_ok=True)

    # ---- discover samples ------------------------------------------------
    sample_dirs = find_sample_dirs(diffab_dir)
    if not sample_dirs:
        raise RuntimeError(f"No sample directories (0000, 0001, …) in {diffab_dir}")
    print(f"Found {len(sample_dirs)} sample dir(s): {sample_dirs}")

    # ---- collect all target pdb_names ------------------------------------
    all_targets = set()
    for sd in sample_dirs:
        all_targets |= extract_pdb_names(os.path.join(diffab_dir, sd))

    if args.test_json:
        with open(args.test_json) as f:
            test_set = set(json.load(f))
        before = len(all_targets)
        all_targets &= test_set
        print(f"Filtered to {len(all_targets)} targets (was {before}, "
              f"test_json has {len(test_set)})")

    print(f"Total targets: {len(all_targets)}")

    # ---- reorganise predictions (with chain remapping) -------------------
    pred_count = 0
    for sample_idx, sd_name in enumerate(sample_dirs):
        sd_path = os.path.join(diffab_dir, sd_name)
        for pdb_name in sorted(all_targets):
            target_dir = os.path.join(pred_out, pdb_name)
            os.makedirs(target_dir, exist_ok=True)
            chain_map = build_chain_mapping(pdb_name)

            # unrelaxed
            src = os.path.join(sd_path, f"{pdb_name}.pdb")
            if os.path.isfile(src):
                dst = os.path.join(target_dir, f"{pdb_name}_model_{sample_idx}.pdb")
                remap_pdb_chains(src, dst, chain_map)
                pred_count += 1

            # relaxed
            src_relax = os.path.join(sd_path, f"{pdb_name}_relaxed.pdb")
            if os.path.isfile(src_relax):
                dst_relax = os.path.join(
                    target_dir, f"{pdb_name}_model_{sample_idx}_relaxed.pdb"
                )
                remap_pdb_chains(src_relax, dst_relax, chain_map)

    print(f"Wrote {pred_count} prediction PDBs (chain-remapped) "
          f"into {len(os.listdir(pred_out))} target dirs")

    # ---- reorganise references (keep original chain IDs) -----------------
    #   eval_metric.py reads references with original chain IDs from pdb_name
    if args.ref_dir:
        ref_count = 0
        for pdb_name in sorted(all_targets):
            src = os.path.join(args.ref_dir, f"{pdb_name}_reference.pdb")
            if os.path.isfile(src):
                dst = os.path.join(ref_out, f"{pdb_name}_reference.pdb")
                link_or_copy(src, dst)
                ref_count += 1
            else:
                print(f"  WARNING: reference not found for {pdb_name}")
    else:
        ref_src_dir = os.path.join(diffab_dir, "reference")
        ref_count = 0
        if os.path.isdir(ref_src_dir):
            for pdb_name in sorted(all_targets):
                src = os.path.join(ref_src_dir, f"{pdb_name}.pdb")
                if os.path.isfile(src):
                    dst = os.path.join(ref_out, f"{pdb_name}_reference.pdb")
                    link_or_copy(src, dst)
                    ref_count += 1
                else:
                    print(f"  WARNING: reference not found for {pdb_name}")
        else:
            print(f"  WARNING: no reference directory at {ref_src_dir}")

    print(f"Linked {ref_count} reference PDBs")

    # ---- write target list JSON ------------------------------------------
    target_json_path = os.path.join(args.output_dir, "diffab_targets.json")
    with open(target_json_path, "w") as f:
        json.dump(sorted(all_targets), f, indent=2)
    print(f"Wrote target list ({len(all_targets)} entries) → {target_json_path}")

    # ---- print eval commands ---------------------------------------------
    print("\n" + "=" * 60)
    print("Ready! Run evaluation with:")
    print("=" * 60)
    print(f"""
# Unrelaxed designs
python evaluate/AbX_eval/eval_metric.py \\
    --data_dir {pred_out} \\
    --ref_dir  {ref_out} \\
    --model boltz \\
    --suffix _model_0.pdb \\
    --test_yaml_dir data/test_yaml_dir/ab \\
    -test_json_fpath {target_json_path} \\
    -c 40

# Relaxed designs
python evaluate/AbX_eval/eval_metric.py \\
    --data_dir {pred_out} \\
    --ref_dir  {ref_out} \\
    --model boltz \\
    --suffix _model_0_relaxed.pdb \\
    --test_yaml_dir data/test_yaml_dir/ab \\
    -test_json_fpath {target_json_path} \\
    -c 40
""")


if __name__ == "__main__":
    main()
