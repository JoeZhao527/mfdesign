#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Prepare RAbD benchmark data for the mfdesign preprocessing pipeline.

Reads RAbD MMAP binary data and generates all inputs required by
``scripts/process/antibody.py`` and ``scripts/process/convert_msa.py``:

  <output_dir>/rabd_processed.csv        -- Processed CSV (antibody.py input)
  <output_dir>/summary.json              -- Metadata JSON (convert_msa.py input)
  <output_dir>/yamls/{entry}.yaml        -- YAML entity-chain mapping files
  <output_dir>/pdbs/{pdb_id}.pdb         -- PDB structure files (from GT CIF)

Usage:
    python rabd_prepare/prepare_rabd_for_mfdesign.py \\
        --data_dir rabd/ \\
        --rabd_eval_dir rabd_eval/ \\
        --output_dir rabd_eval/mfdesign_input/
"""

import argparse
import io
import gzip
import json
import mmap
import os

import gemmi
import numpy as np
import pandas as pd
from ruamel.yaml import YAML


# ============================================================
#  Constants
# ============================================================

AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
STANDARD_AA_3 = set(AA_3TO1.keys())


# ============================================================
#  MMAP I/O  (reused from convert_rabd_to_jsonl.py)
# ============================================================

def decompress(compressed_x: bytes):
    buf = io.BytesIO(compressed_x)
    with gzip.GzipFile(fileobj=buf, mode="rb") as f:
        serialized = f.read().decode()
    return json.loads(serialized)


def read_test_index(data_dir: str):
    idx_path = os.path.join(data_dir, "test_index.txt")
    entries = []
    with open(idx_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            sample_id, start, end = parts[0], int(parts[1]), int(parts[2])
            props = json.loads(parts[3])
            entries.append((sample_id, start, end, props))
    return entries


def load_complex_raw(data_dir: str, byte_start: int, byte_end: int):
    """Load and decompress raw complex tuple from the MMAP binary file."""
    data_bin = os.path.join(data_dir, "data.bin")
    with open(data_bin, "rb") as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        raw = decompress(mm[byte_start:byte_end])
        mm.close()
    return raw


def extract_chain_sequences(raw_complex) -> dict:
    """Extract 1-letter sequences for all chains from raw complex tuple.

    raw_complex layout: (name, molecules, bonds, properties)
    molecule layout:    (name, blocks, chain_id, properties)
    block layout:       (residue_name_3letter, atoms, id, properties)
    """
    sequences = {}
    for mol in raw_complex[1]:
        chain_id = mol[2]
        seq = "".join(AA_3TO1[b[0]] for b in mol[1] if b[0] in STANDARD_AA_3)
        if seq:
            sequences[chain_id] = seq
    return sequences


# ============================================================
#  Helpers
# ============================================================

def build_masked_seq(sequence: str, mark: str) -> str:
    """Replace CDR positions (mark digits != '0') with 'X'."""
    if len(sequence) != len(mark):
        raise ValueError(
            f"Sequence/mark length mismatch: {len(sequence)} vs {len(mark)}"
        )
    return "".join("X" if m != "0" else s for s, m in zip(sequence, mark))


def build_entry_name(sample_id: str, props: dict) -> str:
    """Build mfdesign-style entry name: {pdb}_{heavy}_{light}_{antigen}."""
    pdb_id = sample_id.split("_", 1)[0]
    heavy = props["heavy_chain_id"]
    light = props.get("light_chain_id") or ""
    antigen = "".join(props["target_chain_ids"])
    return f"{pdb_id}_{heavy}_{light}_{antigen}"


# ============================================================
#  YAML generation
# ============================================================

def generate_yaml(entry_name, heavy_id, heavy_seq, light_id, light_seq,
                  antigen_chain_ids, antigen_seq, output_dir):
    """Generate a Boltz-compatible YAML file for one entry."""
    yaml_content = {"version": 1, "sequences": []}

    yaml_content["sequences"].append(
        {"protein": {"id": heavy_id, "sequence": heavy_seq}}
    )
    if light_id:
        yaml_content["sequences"].append(
            {"protein": {"id": light_id, "sequence": light_seq}}
        )
    for ag_id in antigen_chain_ids:
        if ag_id in antigen_seq:
            yaml_content["sequences"].append(
                {"protein": {"id": ag_id, "sequence": antigen_seq[ag_id]}}
            )

    yaml = YAML()
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.width = 4096
    yaml_path = os.path.join(output_dir, f"{entry_name}.yaml")
    with open(yaml_path, "w") as f:
        yaml.dump(yaml_content, f)


# ============================================================
#  CIF -> PDB conversion
# ============================================================

def convert_cif_to_pdb(cif_path: str, pdb_path: str):
    """Convert a CIF file to PDB format using gemmi."""
    structure = gemmi.read_structure(cif_path)
    structure.write_pdb(pdb_path)


# ============================================================
#  Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Prepare RAbD data for mfdesign preprocessing"
    )
    parser.add_argument(
        "--data_dir", type=str, required=True,
        help="RAbD data directory (containing data.bin, test_index.txt)",
    )
    parser.add_argument(
        "--rabd_eval_dir", type=str, required=True,
        help="rabd_eval directory (containing cifs/)",
    )
    parser.add_argument(
        "--output_dir", type=str, default="rabd_eval/mfdesign_input/",
        help="Output directory for all generated files",
    )
    args = parser.parse_args()

    # Create output directories
    yaml_dir = os.path.join(args.output_dir, "yamls")
    pdb_dir = os.path.join(args.output_dir, "pdbs")
    os.makedirs(yaml_dir, exist_ok=True)
    os.makedirs(pdb_dir, exist_ok=True)

    # Read test index
    entries = read_test_index(args.data_dir)
    print(f"[INFO] Read {len(entries)} test entries from {args.data_dir}")

    csv_rows = []
    summary_dict = {}
    converted_pdbs = set()

    for i, (sample_id, byte_start, byte_end, props) in enumerate(entries):
        pdb_id = sample_id.split("_", 1)[0]
        entry_name = build_entry_name(sample_id, props)

        # ---- Sequences from test_index.txt ----
        heavy_id = props["heavy_chain_id"]
        heavy_seq = props["heavy_chain_sequence"]
        heavy_mark = props["heavy_chain_mark"]

        light_id = props.get("light_chain_id")
        light_seq = props.get("light_chain_sequence")
        light_mark = props.get("light_chain_mark")

        antigen_chain_ids = props["target_chain_ids"]

        # ---- Antigen sequences from MMAP Complex ----
        raw_cplx = load_complex_raw(args.data_dir, byte_start, byte_end)
        chain_seqs = extract_chain_sequences(raw_cplx)

        antigen_seq = {}
        for ag_id in antigen_chain_ids:
            if ag_id in chain_seqs:
                antigen_seq[ag_id] = chain_seqs[ag_id]
            else:
                print(f"  [WARN] Antigen chain {ag_id} not found in {sample_id}")

        # ---- Masked sequences (all CDR positions -> X) ----
        heavy_masked = build_masked_seq(heavy_seq, heavy_mark)
        light_masked = None
        if light_seq and light_mark:
            light_masked = build_masked_seq(light_seq, light_mark)

        # ---- CSV row (matches antibody.py fetch() expectations) ----
        csv_rows.append({
            "file_name": entry_name,
            "pdb": pdb_id,
            "H_chain_id": heavy_id,
            "L_chain_id": light_id if light_id else np.nan,
            "H_chain_seq": heavy_seq,
            "L_chain_seq": light_seq if light_seq else np.nan,
            "H_chain_masked_seq": heavy_masked,
            "L_chain_masked_seq": light_masked if light_masked else np.nan,
            "antigen_chain_id": str(antigen_chain_ids),
            "antigen_seq": str(antigen_seq),
            "resolution": 0.0,
            "scfv": False,
        })

        # ---- YAML (entity-chain mapping for Boltz) ----
        generate_yaml(
            entry_name, heavy_id, heavy_seq,
            light_id, light_seq,
            antigen_chain_ids, antigen_seq,
            yaml_dir,
        )

        # ---- CIF -> PDB (once per pdb_id, from GT structures) ----
        if pdb_id not in converted_pdbs:
            cif_path = os.path.join(args.rabd_eval_dir, "cifs", f"{pdb_id}.cif")
            pdb_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
            if os.path.exists(cif_path):
                convert_cif_to_pdb(cif_path, pdb_path)
                converted_pdbs.add(pdb_id)
            else:
                print(f"  [WARN] GT CIF not found: {cif_path}")

        # ---- summary.json entry (for convert_msa.py) ----
        summary_dict[entry_name] = {
            "file_name": entry_name,
            "pdb": pdb_id,
            "H_chain_id": heavy_id,
            "L_chain_id": light_id,
            "H_chain_seq": heavy_seq,
            "L_chain_seq": light_seq,
            "H_chain_masked_seq": heavy_masked,
            "L_chain_masked_seq": light_masked,
            "antigen_chain_id": antigen_chain_ids,
            "antigen_seq": antigen_seq,
        }

        print(f"  [{i+1}/{len(entries)}] {entry_name}")

    # ---- Write CSV ----
    df = pd.DataFrame(csv_rows).set_index("file_name")
    csv_path = os.path.join(args.output_dir, "rabd_processed.csv")
    df.to_csv(csv_path)

    # ---- Write summary.json ----
    summary_path = os.path.join(args.output_dir, "summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary_dict, f, indent=2)

    print(f"\n[DONE] CSV:     {csv_path} ({len(df)} entries)")
    print(f"[DONE] Summary: {summary_path}")
    print(f"[DONE] YAMLs:   {yaml_dir}/ ({len(os.listdir(yaml_dir))} files)")
    print(f"[DONE] PDBs:    {pdb_dir}/ ({len(os.listdir(pdb_dir))} files)")


if __name__ == "__main__":
    main()
