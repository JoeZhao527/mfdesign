"""
Convert PDB files to CIF format using BioPython.

Usage:
    python prepare_protenix_reference.py input_dir output_dir

Example:
    python prepare_protenix_reference.py test_entry_pdb_files output_cif_files
"""

import os
import sys
from pathlib import Path
from Bio.PDB import PDBParser, MMCIFIO


def pdb_to_cif(input_dir: str, output_dir: str):
    """Convert all PDB files in input_dir to CIF format in output_dir."""
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    parser = PDBParser(QUIET=True)
    io = MMCIFIO()

    pdb_files = list(input_path.glob("*_reference.pdb"))
    print(f"Found {len(pdb_files)} PDB files to convert")

    for pdb_file in pdb_files:
        # Remove _reference suffix: 7vgs_D_C_A_reference.pdb -> 7vgs_D_C_A.cif
        base_name = pdb_file.stem.replace("_reference", "")
        output_file = output_path / f"{base_name}.cif"

        structure = parser.get_structure(base_name, pdb_file)
        io.set_structure(structure)
        io.save(str(output_file))

        print(f"Converted: {pdb_file.name} -> {output_file.name}")

    print(f"\nDone! Converted {len(pdb_files)} files to {output_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python prepare_protenix_reference.py input_dir output_dir")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    pdb_to_cif(input_dir, output_dir)
