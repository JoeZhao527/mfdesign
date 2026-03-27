#!/usr/bin/env python3
"""
Convert BoltzGen test_hcdr3 PDB bundle to eval_metric.py boltz format.

BoltzGen bundle layout
----------------------
  {bundle_dir}/{sample_id}/
    {sample_id}_predicted_inverse_folded.pdb   (chains A,B,C = antigen,heavy,light)
    {sample_id}_ground_truth_native.pdb        (original chain IDs)

Naming conventions
------------------
  BoltzGen:      {code}_{antigen}_{heavy}_{light}    e.g. 1a14_N_H_L
  eval_metric:   {code}_{heavy}_{light}_{antigen}    e.g. 1a14_H_L_N

Chain ID mapping
----------------
  BoltzGen prediction PDB: A=1st(antigen), B=2nd(heavy), C=3rd(light)
  eval_metric --model boltz expects: A=heavy, B=light
  Remap: B→A, C→B, A→C

Output
------
  {output_dir}/
    predictions/{new_name}/{new_name}_model_0.pdb   (chain-remapped)
    references/{new_name}_reference.pdb              (symlink, original chain IDs)
    yaml/{new_name}.yaml                             (CDR masks via ANARCI)
    boltzgen_targets.json

Usage
-----
  python r1_eval/prepare_boltzgen_eval.py \\
      --bundle_dir test_hcdr3_num1_pdb_bundle_20260327 \\
      --output_dir test_hcdr3_num1_pdb_bundle_20260327/eval_ready
"""

import argparse
import csv
import json
import os
import sys
import io
import numpy as np

# abx lives under evaluate/AbX_eval/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "evaluate", "AbX_eval"))

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from abx.preprocess.numbering import renumber_ab_seq, get_ab_regions


# ---------------------------------------------------------------------------
# Name conversion
# ---------------------------------------------------------------------------

def swap_name(sample_id):
    """BoltzGen {code}_{ag}_{hv}_{lt} → eval {code}_{hv}_{lt}_{ag}."""
    parts = sample_id.split("_")
    # parts: [code, antigen, heavy, light]
    return f"{parts[0]}_{parts[2]}_{parts[3]}_{parts[1]}"


# ---------------------------------------------------------------------------
# PDB chain remapping
# ---------------------------------------------------------------------------

# BoltzGen prediction: A=antigen, B=heavy, C=light
# eval_metric boltz:   A=heavy,   B=light, C=antigen
_PRED_CHAIN_REMAP = {"A": "C", "B": "A", "C": "B"}


def remap_pdb_chains(src_path, dst_path):
    """Rewrite prediction PDB with chain IDs remapped for eval_metric."""
    with open(src_path) as f:
        lines = f.readlines()
    with open(dst_path, "w") as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM", "TER", "ANISOU")) and len(line) > 21:
                old = line[21]
                if old in _PRED_CHAIN_REMAP:
                    line = line[:21] + _PRED_CHAIN_REMAP[old] + line[22:]
            f.write(line)


# ---------------------------------------------------------------------------
# Sequence extraction
# ---------------------------------------------------------------------------

def extract_chain_seq(pdb_path, chain_id):
    """Extract single-letter amino-acid sequence for one chain."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    residues = list(structure[0][chain_id].get_residues())
    return "".join(seq1(r.get_resname()) for r in residues)


# ---------------------------------------------------------------------------
# YAML CDR mask generation
# ---------------------------------------------------------------------------

def _spec_mask_for_chain(sequence, chain_type):
    """Generate binary CDR spec_mask via ANARCI IMGT numbering.

    Returns a string of '0'/'1' with len == len(sequence), or None on failure.
    chain_type: 'H' or 'L'.
    """
    allow = ["H"] if chain_type == "H" else ["K", "L"]

    # renumber_ab_seq has noisy prints; suppress them
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    try:
        anarci_res = renumber_ab_seq(sequence, allow=allow, scheme="imgt")
    finally:
        sys.stdout = old_stdout

    domain_numbering = anarci_res.get("domain_numbering")
    start = anarci_res.get("start")
    if domain_numbering is None:
        return None

    region_def = get_ab_regions(domain_numbering, chain_id=chain_type)

    # Heavy CDR values: 1(H1), 3(H2), 5(H3)
    # Light CDR values: 8(L1), 10(L2), 12(L3)
    cdr_values = {1, 3, 5} if chain_type == "H" else {8, 10, 12}

    mask = ["0"] * len(sequence)
    for i, val in enumerate(region_def):
        if val in cdr_values:
            mask[start + i] = "1"

    return "".join(mask)


def write_yaml(path, heavy_id, light_id, heavy_seq, light_seq, heavy_mask, light_mask):
    """Write a minimal YAML CDR mask file compatible with read_yaml_file()."""
    with open(path, "w") as f:
        f.write("version: 1\nsequences:\n")
        f.write("  - protein:\n")
        f.write(f"      id: {heavy_id}\n")
        f.write(f"      sequence: {heavy_seq}\n")
        f.write(f"      spec_mask: '{heavy_mask}'\n")
        f.write(f"      ground_truth: {heavy_seq}\n")
        if light_id:
            f.write("  - protein:\n")
            f.write(f"      id: {light_id}\n")
            f.write(f"      sequence: {light_seq}\n")
            f.write(f"      spec_mask: '{light_mask}'\n")
            f.write(f"      ground_truth: {light_seq}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert BoltzGen PDB bundle to eval_metric.py format"
    )
    parser.add_argument("--bundle_dir", type=str, required=True,
                        help="Path to the PDB bundle directory")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Output directory for eval-ready files")
    args = parser.parse_args()

    # ---- parse manifest --------------------------------------------------
    manifest = os.path.join(args.bundle_dir, "manifest.csv")
    targets = []
    with open(manifest) as f:
        for row in csv.DictReader(f):
            targets.append(row)
    print(f"Found {len(targets)} targets in manifest")

    # ---- create output dirs ----------------------------------------------
    pred_out = os.path.join(args.output_dir, "predictions")
    ref_out = os.path.join(args.output_dir, "references")
    yaml_out = os.path.join(args.output_dir, "yaml")
    for d in [pred_out, ref_out, yaml_out]:
        os.makedirs(d, exist_ok=True)

    new_names = []
    errors = []

    for target in targets:
        sample_id = target["sample_id"]
        new_name = swap_name(sample_id)
        parts = new_name.split("_")
        heavy_id, light_id = parts[1], parts[2]  # original chain IDs

        pred_src = os.path.join(args.bundle_dir, target["predicted_pdb"])
        gt_src = os.path.join(args.bundle_dir, target["ground_truth_pdb"])

        if not os.path.isfile(pred_src) or not os.path.isfile(gt_src):
            errors.append(f"{sample_id}: PDB file(s) missing")
            continue

        # 1. Prediction: remap chains
        tgt_dir = os.path.join(pred_out, new_name)
        os.makedirs(tgt_dir, exist_ok=True)
        remap_pdb_chains(pred_src, os.path.join(tgt_dir, f"{new_name}_model_0.pdb"))

        # 2. Reference: symlink GT (original chain IDs match new_name)
        ref_dst = os.path.join(ref_out, f"{new_name}_reference.pdb")
        if os.path.exists(ref_dst) or os.path.islink(ref_dst):
            os.remove(ref_dst)
        os.symlink(os.path.abspath(gt_src), ref_dst)

        # 3. YAML: generate CDR masks from GT sequence via ANARCI
        #    On failure, clean up prediction/reference created above.
        def _cleanup_failed():
            import shutil
            shutil.rmtree(tgt_dir, ignore_errors=True)
            if os.path.exists(ref_dst) or os.path.islink(ref_dst):
                os.remove(ref_dst)

        try:
            heavy_seq = extract_chain_seq(gt_src, heavy_id)
            heavy_mask = _spec_mask_for_chain(heavy_seq, "H")
            if heavy_mask is None:
                errors.append(f"{sample_id}: ANARCI failed on heavy chain ({heavy_id})")
                _cleanup_failed()
                continue

            if light_id:
                light_seq = extract_chain_seq(gt_src, light_id)
                light_mask = _spec_mask_for_chain(light_seq, "L")
                if light_mask is None:
                    errors.append(f"{sample_id}: ANARCI failed on light chain ({light_id})")
                    _cleanup_failed()
                    continue
            else:
                light_seq, light_mask = "", ""

            write_yaml(
                os.path.join(yaml_out, f"{new_name}.yaml"),
                heavy_id, light_id, heavy_seq, light_seq, heavy_mask, light_mask,
            )
        except Exception as e:
            errors.append(f"{sample_id}: {e}")
            continue

        new_names.append(new_name)

    # ---- write target list -----------------------------------------------
    json_path = os.path.join(args.output_dir, "boltzgen_targets.json")
    with open(json_path, "w") as f:
        json.dump(sorted(new_names), f, indent=2)

    # ---- summary ---------------------------------------------------------
    print(f"Processed {len(new_names)} targets OK, {len(errors)} errors")
    for e in errors:
        print(f"  ERROR: {e}")

    print(f"\n{'=' * 60}")
    print("Run evaluation:")
    print("=" * 60)
    print(f"""
python evaluate/AbX_eval/eval_metric.py \\
    --data_dir {pred_out} \\
    --ref_dir  {ref_out} \\
    --test_yaml_dir {yaml_out} \\
    -test_json_fpath {json_path} \\
    --model boltz \\
    --suffix _model_0.pdb \\
    -c 40
""")


if __name__ == "__main__":
    main()
