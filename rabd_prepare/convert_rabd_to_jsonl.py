#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Convert RabD (Rosetta Antibody Design) benchmark data from UniMoMo's MMAP binary
format into the CIF + JSONL + cdr_info format expected by the cdr_eval pipeline.

Two modes:

  --prepare (default):
    Reads data.bin/test_index.txt, produces:
      <output_dir>/cifs/{pdb_id}.cif              -- GT CIF files (one per complex)
      <output_dir>/cdr_info/*_cdr_info.json        -- CDR info JSON files
      <output_dir>/rabd_cdr_eval_split.json        -- Split JSON for cdr_masking
      <output_dir>/rabd_h3_test.jsonl              -- JSONL with source_cif = {pdb_id}.cif

  --finalize:
    After Protenix prediction + remap, produces:
      <output_dir>/rabd_h3_test.jsonl              -- JSONL with source_cif = best remap CIF

Usage:
    # Phase 1: Prepare
    python -m lmms_engine.pyscripts.convert_rabd_to_jsonl \
        --data_dir /path/to/sabdab_processed/processed \
        --output_dir data/rabd_eval/

    # Phase 4: Finalize (after protenix predict + remap)
    python -m lmms_engine.pyscripts.convert_rabd_to_jsonl \
        --data_dir /path/to/sabdab_processed/processed \
        --output_dir data/rabd_eval/ \
        --finalize \
        --predictions_dir data/rabd_eval/predictions/ \
        --remap_dir data/rabd_eval/predictions_remap/ \
        --seed 101
"""

import argparse
import io
import gzip
import json
import mmap
import os
from copy import copy
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Tuple

import numpy as np

# [ref] docs/claude/implements/rabd-h3-eval-pipeline-plan.md#ADR-5
# Standard amino acid 3-letter codes for hetero detection
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}


# ============================================================
#  Data classes (minimal subset from UniMoMo bioparse/hierarchy.py)
# ============================================================
# [src] proteinfm_joint_train/scripts/read_unimomo_antibody_test.py:L132-L241


class BondType(Enum):
    NONE = 0
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 4


@dataclass
class Bond:
    index1: Tuple[int, int, int]
    index2: Tuple[int, int, int]
    bond_type: BondType

    @classmethod
    def from_tuple(cls, tup):
        return Bond(tuple(tup[0]), tuple(tup[1]), BondType(tup[2]))


class Atom:
    def __init__(self, name, coordinate, element, id, properties=None):
        self.name = name
        self.coordinate = coordinate
        self.element = element
        self.id = id
        self.properties = properties or {}

    def get_element(self):
        return self.element

    def get_coord(self):
        return copy(self.coordinate)

    @classmethod
    def from_tuple(cls, data):
        return Atom(data[0], data[1], data[2], data[3], data[4])


class Block:
    def __init__(self, name, atoms, id, properties=None):
        self.name = name
        self.atoms = atoms
        self.id = tuple(id)
        self.properties = properties or {}

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    @classmethod
    def from_tuple(cls, data):
        return Block(data[0], [Atom.from_tuple(a) for a in data[1]], tuple(data[2]), data[3])


class Molecule:
    def __init__(self, name, blocks, id, properties=None):
        self.name = name
        self.blocks = blocks
        self.id = id
        self.properties = properties or {}

    def __len__(self):
        return len(self.blocks)

    def __iter__(self):
        return iter(self.blocks)

    @classmethod
    def from_tuple(cls, data):
        return Molecule(data[0], [Block.from_tuple(b) for b in data[1]], data[2], data[3])


class Complex:
    def __init__(self, name, molecules, bonds, properties=None):
        self.name = name
        self.molecules = molecules
        self.bonds = bonds
        self.properties = properties or {}

    def __len__(self):
        return len(self.molecules)

    def __iter__(self):
        return iter(self.molecules)

    @classmethod
    def from_tuple(cls, data):
        return Complex(
            data[0],
            [Molecule.from_tuple(m) for m in data[1]],
            [Bond.from_tuple(b) for b in data[2]],
            data[3] if len(data) > 3 else {},
        )


# ============================================================
#  MMAP I/O
# ============================================================
# [src] proteinfm_joint_train/scripts/read_unimomo_antibody_test.py:L269-L273


def decompress(compressed_x: bytes):
    """Decompress gzipped JSON bytes into a Python object."""
    buf = io.BytesIO(compressed_x)
    with gzip.GzipFile(fileobj=buf, mode="rb") as f:
        serialized = f.read().decode()
    return json.loads(serialized)


def read_test_index(data_dir: str) -> List[Tuple[str, int, int, dict]]:
    """Read test_index.txt and return list of (sample_id, byte_start, byte_end, properties)."""
    idx_path = os.path.join(data_dir, "test_index.txt")
    entries = []
    with open(idx_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            sample_id, start, end = parts[0], int(parts[1]), int(parts[2])
            props = json.loads(parts[3])
            entries.append((sample_id, start, end, props))
    return entries


def load_complex(data_dir: str, byte_start: int, byte_end: int) -> Complex:
    """Load and decompress a single Complex from the MMAP binary file."""
    data_bin = os.path.join(data_dir, "data.bin")
    with open(data_bin, "rb") as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        raw = decompress(mm[byte_start:byte_end])
        mm.close()
    return Complex.from_tuple(raw)


# ============================================================
#  Complex -> Biotite AtomArray -> CIF via Protenix CIFWriter
# ============================================================
# [ref] docs/claude/implements/rabd-h3-eval-pipeline-plan.md#CIF-Export-Strategy


def complex_to_atom_array(cplx: Complex):
    """Convert a Complex object to a biotite AtomArray.

    Iterates molecules (chains) -> blocks (residues) -> atoms.
    Skips hydrogen atoms. Sets sequential res_id from 1 per chain.
    """
    import biotite.structure as struc
    from biotite.structure import Atom as BiotiteAtom, AtomArray

    atoms = []
    for mol in cplx:
        res_id = 0
        for block in mol:
            # Skip non-standard residues (metal ions like [Mg], ligands, etc.)
            # These cause Protenix CCD parsing failures downstream
            if block.name not in STANDARD_AA:
                continue
            res_id += 1
            for atom in block:
                if atom.get_element() == "H":
                    continue
                coord = atom.get_coord()
                ba = BiotiteAtom(
                    coord=coord,
                    chain_id=mol.id,
                    res_id=res_id,
                    res_name=block.name,
                    atom_name=atom.name,
                    element=atom.get_element(),
                    hetero=False,
                )
                atoms.append(ba)

    atom_array = struc.array(atoms)

    # Attach an empty BondList (required by CIFWriter's bond filtering logic)
    atom_array.bonds = struc.BondList(len(atom_array))

    # Assign label_entity_id: one entity per unique chain
    entity_map = {}
    entity_counter = 0
    label_entity_id = np.empty(len(atom_array), dtype="<U4")
    for i, cid in enumerate(atom_array.chain_id):
        if cid not in entity_map:
            entity_counter += 1
            entity_map[cid] = str(entity_counter)
        label_entity_id[i] = entity_map[cid]
    atom_array.set_annotation("label_entity_id", label_entity_id)

    # Set label_asym_id = chain_id, label_atom_id = atom_name, label_seq_id = res_id
    atom_array.set_annotation("label_asym_id", atom_array.chain_id.copy())
    atom_array.set_annotation("label_atom_id", atom_array.atom_name.copy())
    atom_array.set_annotation("label_seq_id", atom_array.res_id.copy())

    return atom_array


def build_entity_poly_type(atom_array) -> Dict[str, str]:
    """Build entity_poly_type dict: entity_id -> 'polypeptide(L)' for protein chains."""
    import biotite.structure as struc

    entity_poly_type = {}
    entity_ids = np.unique(atom_array.label_entity_id)
    for eid in entity_ids:
        mask = atom_array.label_entity_id == eid
        chain_array = atom_array[mask]
        res_starts = struc.get_residue_starts(chain_array, add_exclusive_stop=False)
        res_names = chain_array[res_starts].res_name
        protein_count = sum(1 for name in res_names if name in STANDARD_AA)
        if protein_count >= 2:
            entity_poly_type[eid] = "polypeptide(L)"
    return entity_poly_type


def complex_to_protenix_cif(cplx: Complex, output_path: str, entry_id: str):
    """Export a Complex to a Protenix-compatible CIF file using CIFWriter.

    Falls back to basic biotite CIF writing if CIFWriter is unavailable.
    """
    atom_array = complex_to_atom_array(cplx)
    entity_poly_type = build_entity_poly_type(atom_array)

    try:
        import _setup_protenix  # noqa: F401
        from protenix.data.utils import CIFWriter

        writer = CIFWriter(atom_array=atom_array, entity_poly_type=entity_poly_type)
        writer.save_to_cif(output_path, entry_id=entry_id)
    except ImportError:
        import biotite.structure.io.pdbx as pdbx

        cif_file = pdbx.CIFFile()
        pdbx.set_structure(cif_file, atom_array, data_block=entry_id)
        cif_file.write(output_path)
        print("[WARN] CIFWriter unavailable, using basic biotite CIF output (may lack entity_poly metadata)")


# ============================================================
#  JSONL generation
# ============================================================
# [ref] docs/claude/implements/rabd-h3-eval-pipeline-plan.md#ADR-3


def build_jsonl_entry(sample_id: str, props: dict, cdr_type: str, source_cif: Optional[str] = None) -> dict:
    """Build one JSONL entry for a RabD test sample.

    Args:
        sample_id: e.g. "4cmh_A_B_C"
        props: properties dict from test_index.txt
        cdr_type: CDR type string, e.g. "H3"
        source_cif: Override for source_cif filename (default: "{pdb_id}.cif")
    """
    # Extract CDR sequence from the mark
    mark = props["heavy_chain_mark"]
    cdr_digit = cdr_type[-1]  # "3" from "H3"
    h3_indices = [i for i, c in enumerate(mark) if c == cdr_digit]
    h3_start = h3_indices[0]
    h3_end = h3_indices[-1] + 1
    h3_seq = props["heavy_chain_sequence"][h3_start:h3_end]
    h3_len = h3_end - h3_start

    # input_text: no design_points (ADR-3)
    input_text = (
        "<TASK=AB_CDR_REDESIGN_SFT_V2>\n"
        "You are redesigning masked CDR regions of an antibody to improve binding to the antigen.\n"
        "Output JSON with keys: task, thinking, answer."
    )

    # target_text with H3-only (ADR-4)
    cdr_tag = f"<HCDR{cdr_digit}>"
    cdr_close_tag = f"</HCDR{cdr_digit}>"
    filled_positions = [{"pos": i + 1, "aa": h3_seq[i]} for i in range(h3_len)]
    target_obj = {
        "task": "AB_CDR_REDESIGN_SFT_V2",
        "thinking": "",
        "answer": {
            "cdrs_present": [cdr_type],
            "cdr_sequences": {
                cdr_type: {
                    "len": h3_len,
                    "seq": f"{cdr_tag}{h3_seq}{cdr_close_tag}",
                    "filled_positions": filled_positions,
                }
            },
        },
    }
    target_text = json.dumps(target_obj)

    # meta
    pdb_id = sample_id.split("_", 1)[0]
    meta = {
        "pdb_id": pdb_id,
        "entry": sample_id,
        "heavy_chain_id": props["heavy_chain_id"],
        "light_chain_id": props["light_chain_id"],
        "antigen_chain_ids": props["target_chain_ids"],
        "cdrs_present": [cdr_type],
        "source_cif": source_cif if source_cif is not None else f"{pdb_id}.cif",
    }

    return {
        "sample_id": f"{sample_id}_rabd_{cdr_type.lower()}",
        "task": "AB_CDR_REDESIGN_SFT_V2",
        "input_text": input_text,
        "target_text": target_text,
        "meta": meta,
    }


# ============================================================
#  cdr_info generation
# ============================================================
# [ref] docs/claude/implements/rabd-h3-eval-pipeline-plan.md#ADR-5


def build_cdr_info(sample_id: str, props: dict, cplx: Complex, cdr_type: str) -> dict:
    """Build cdr_info dict from RabD metadata and Complex structure.

    The indices are 0-based positions in the heavy chain where mark == cdr_digit.
    These correspond directly to label_seq_id - 1 in the CIF (sequential from 1).
    """
    protein_id = sample_id.upper()
    mark = props["heavy_chain_mark"]
    cdr_digit = cdr_type[-1]
    cdr_indices = [i for i, c in enumerate(mark) if c == cdr_digit]

    # Build chain_mapping from the Complex's molecule order (identity mapping for RabD)
    sequence_order = []
    protenix_to_original = {}
    original_to_protenix = {}
    for mol in cplx:
        chain_id = mol.id
        seq_len = len(mol.blocks)
        protenix_to_original[chain_id] = chain_id
        original_to_protenix[chain_id] = chain_id
        sequence_order.append({
            "protenix_chain": chain_id,
            "original_chain": chain_id,
            "seq_len": seq_len,
            "entity_type": "proteinChain",
        })

    cdr_key = f"cdr{cdr_digit}"
    return {
        "entry_name": protein_id,
        "heavy_chain_id": props["heavy_chain_id"],
        "light_chain_id": props["light_chain_id"],
        "chain_mapping": {
            "protenix_to_original": protenix_to_original,
            "original_to_protenix": original_to_protenix,
            "sequence_order": sequence_order,
        },
        "cdr_info": {
            "H_chain": {
                "cdr_indices": cdr_indices,
                cdr_key: {
                    "indices": cdr_indices,
                },
            }
        },
    }


# ============================================================
#  Main
# ============================================================


def build_cdr_eval_entry(sample_id: str, props: dict) -> str:
    """Build cdr_eval-format entry name: {pdb}_{heavy}_{light}_{antigen}.

    RabD sample_id format is "{pdb}_{antigen}_{heavy}_{light}" but cdr_eval
    expects "{pdb}_{heavy}_{light}_{antigen}".
    """
    pdb_id = sample_id.split("_", 1)[0]
    heavy = props["heavy_chain_id"]
    light = props["light_chain_id"]
    antigen = "".join(props["target_chain_ids"])
    return f"{pdb_id}_{heavy}_{light}_{antigen}"


def select_best_sample(predictions_dir: str, cdr_eval_entry: str, seed: int) -> Optional[int]:
    """Select the best sample index by ranking_score from Protenix confidence JSONs.

    Args:
        predictions_dir: Base predictions directory (contains entry subdirs)
        cdr_eval_entry: Entry name matching Protenix output dir (e.g. "4cmh_B_C_A")
        seed: Protenix seed number

    Returns:
        Best sample index, or None if no confidence files found.
    """
    seed_dir = os.path.join(predictions_dir, cdr_eval_entry, f"seed_{seed}")
    # Protenix puts outputs in a "predictions" subdirectory
    pred_subdir = os.path.join(seed_dir, "predictions")
    if os.path.isdir(pred_subdir):
        seed_dir = pred_subdir
    if not os.path.isdir(seed_dir):
        return None

    best_idx = None
    best_score = float("-inf")
    for fname in os.listdir(seed_dir):
        if not fname.endswith(".json") or "summary_confidence" not in fname:
            continue
        # Pattern: {entry}_summary_confidence_sample_{N}.json
        try:
            sample_num = int(fname.rsplit("_sample_", 1)[1].replace(".json", ""))
        except (IndexError, ValueError):
            continue
        conf_path = os.path.join(seed_dir, fname)
        with open(conf_path, "r") as f:
            conf_data = json.load(f)
        score = conf_data.get("ranking_score", 0)
        if score > best_score:
            best_score = score
            best_idx = sample_num

    return best_idx


def run_prepare(args):
    """Prepare mode: generate CIFs, cdr_info, split JSON, and preliminary JSONL."""
    cif_dir = os.path.join(args.output_dir, "cifs")
    cdr_info_dir = os.path.join(args.output_dir, "cdr_info")
    split_json_path = os.path.join(args.output_dir, "rabd_cdr_eval_split.json")
    jsonl_path = os.path.join(args.output_dir, f"rabd_{args.cdr_type.lower()}_test.jsonl")

    os.makedirs(cif_dir, exist_ok=True)
    os.makedirs(cdr_info_dir, exist_ok=True)

    entries = read_test_index(args.data_dir)
    print(f"[INFO] Read {len(entries)} test entries from {args.data_dir}")

    jsonl_lines = []
    split_entries = []

    for i, (sample_id, byte_start, byte_end, props) in enumerate(entries):
        cdr_digit = args.cdr_type[-1]
        if cdr_digit not in props["heavy_chain_mark"]:
            print(f"[WARN] Skipping {sample_id}: no {args.cdr_type} in heavy_chain_mark")
            continue

        pdb_id = sample_id.split("_", 1)[0]

        # CIF export (named by pdb_id for cdr_masking compatibility)
        cplx = load_complex(args.data_dir, byte_start, byte_end)
        cif_path = os.path.join(cif_dir, f"{pdb_id}.cif")
        complex_to_protenix_cif(cplx, cif_path, entry_id=pdb_id)

        # cdr_info
        protein_id = sample_id.upper()
        cdr_info = build_cdr_info(sample_id, props, cplx, args.cdr_type)
        cdr_info_path = os.path.join(cdr_info_dir, f"{protein_id}_cdr_info.json")
        with open(cdr_info_path, "w") as f:
            json.dump(cdr_info, f, indent=2)

        # cdr_eval split entry: {pdb}_{heavy}_{light}_{antigen}
        cdr_eval_entry = build_cdr_eval_entry(sample_id, props)
        split_entries.append(cdr_eval_entry)

        # Preliminary JSONL with source_cif = GT CIF
        jsonl_entry = build_jsonl_entry(sample_id, props, args.cdr_type, source_cif=f"{pdb_id}.cif")
        jsonl_lines.append(json.dumps(jsonl_entry))

        print(f"  [{i+1}/{len(entries)}] {sample_id} -> {pdb_id}.cif | {cdr_eval_entry}")

    # Write split JSON
    with open(split_json_path, "w") as f:
        json.dump(split_entries, f, indent=2)

    # Write JSONL
    with open(jsonl_path, "w") as f:
        f.write("\n".join(jsonl_lines) + "\n")

    print(f"\n[DONE] Prepare output:")
    print(f"  CIFs:      {cif_dir}/ ({len(os.listdir(cif_dir))} files)")
    print(f"  cdr_info:  {cdr_info_dir}/ ({len(os.listdir(cdr_info_dir))} files)")
    print(f"  Split JSON:{split_json_path} ({len(split_entries)} entries)")
    print(f"  JSONL:     {jsonl_path} ({len(jsonl_lines)} entries)")


def run_finalize(args):
    """Finalize mode: select best Protenix samples and generate updated JSONL."""
    if args.predictions_dir is None:
        raise ValueError("--predictions_dir is required for --finalize mode")

    jsonl_path = os.path.join(args.output_dir, f"rabd_{args.cdr_type.lower()}_test.jsonl")

    entries = read_test_index(args.data_dir)
    print(f"[INFO] Read {len(entries)} test entries from {args.data_dir}")

    jsonl_lines = []
    for i, (sample_id, byte_start, byte_end, props) in enumerate(entries):
        cdr_digit = args.cdr_type[-1]
        if cdr_digit not in props["heavy_chain_mark"]:
            print(f"[WARN] Skipping {sample_id}: no {args.cdr_type} in heavy_chain_mark")
            continue

        cdr_eval_entry = build_cdr_eval_entry(sample_id, props)

        # Select best sample
        if args.sample_idx is not None:
            best_idx = args.sample_idx
        else:
            best_idx = select_best_sample(args.predictions_dir, cdr_eval_entry, args.seed)
            if best_idx is None:
                print(f"[WARN] No predictions found for {cdr_eval_entry}, defaulting to sample_0")
                best_idx = 0

        source_cif = f"{cdr_eval_entry}_sample_{best_idx}.cif"

        # Verify remap CIF exists if remap_dir provided
        if args.remap_dir is not None:
            remap_path = os.path.join(args.remap_dir, source_cif)
            if not os.path.isfile(remap_path):
                print(f"[WARN] Remap CIF not found: {remap_path}")

        jsonl_entry = build_jsonl_entry(sample_id, props, args.cdr_type, source_cif=source_cif)
        jsonl_lines.append(json.dumps(jsonl_entry))

        print(f"  [{i+1}/{len(entries)}] {cdr_eval_entry} -> {source_cif}")

    with open(jsonl_path, "w") as f:
        f.write("\n".join(jsonl_lines) + "\n")

    print(f"\n[DONE] Finalize output:")
    print(f"  JSONL: {jsonl_path} ({len(jsonl_lines)} entries)")


def main():
    parser = argparse.ArgumentParser(
        description="Convert RabD MMAP data to CIF + JSONL + cdr_info for cdr_eval pipeline"
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        required=True,
        help="Path to RabD processed directory (containing data.bin, test_index.txt)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="data/rabd_eval/",
        help="Output directory for CIFs, JSONL, and cdr_info",
    )
    parser.add_argument(
        "--cdr_type",
        type=str,
        default="H3",
        help="CDR type to evaluate (default: H3)",
    )
    parser.add_argument(
        "--finalize",
        action="store_true",
        help="Finalize mode: read protenix predictions, select best sample, generate JSONL",
    )
    parser.add_argument(
        "--predictions_dir",
        type=str,
        default=None,
        help="Protenix predict output directory (required for --finalize)",
    )
    parser.add_argument(
        "--remap_dir",
        type=str,
        default=None,
        help="Remapped CIF directory from remap_cif_chains.py (for verification only)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=101,
        help="Protenix seed to use (default: 101)",
    )
    parser.add_argument(
        "--sample_idx",
        type=int,
        default=None,
        help="Force specific sample index (default: auto-select best by ranking_score)",
    )
    args = parser.parse_args()

    if args.finalize:
        run_finalize(args)
    else:
        run_prepare(args)


if __name__ == "__main__":
    main()
