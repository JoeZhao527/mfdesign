"""
Prepare AntiBMPNN input from predicted structures and YAML sequence files.

This script:
1. Reads YAML files containing sequences with 'X' marks (positions to design)
2. Copies corresponding PDB files to an input folder
3. Generates the three required JSONL files for AntiBMPNN:
   - parsed_pdbs.jsonl: parsed PDB structure info
   - assigned_pdbs.jsonl: chain assignment (which chains to design)
   - fixed_pdbs.jsonl: fixed positions (positions NOT to design)
"""

import argparse
import os
import json
import shutil
import glob
import numpy as np
import yaml
from pathlib import Path


# Amino acid mappings
ALPHA_1 = list("ARNDCQEGHILKMFPSTWYV-")
ALPHA_3 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'GAP']
AA_3_N = {a: n for n, a in enumerate(ALPHA_3)}
AA_N_1 = {n: a for n, a in enumerate(ALPHA_1)}


def N_to_AA(x):
    """Convert numerical representation to amino acid sequence."""
    x = np.array(x)
    if x.ndim == 1:
        x = x[None]
    return ["".join([AA_N_1.get(a, "-") for a in y]) for y in x]


def parse_PDB_biounits(pdb_file, atoms=['N', 'CA', 'C', 'O'], chain=None):
    """
    Parse PDB file and extract coordinates and sequence for a specific chain.
    
    Returns:
        (coords, sequence) or ('no_chain', 'no_chain') if chain not found
    """
    xyz, seq, min_resn, max_resn = {}, {}, 1e6, -1e6
    
    with open(pdb_file, "rb") as f:
        for line in f:
            line = line.decode("utf-8", "ignore").rstrip()
            
            # Handle MSE (selenomethionine) as MET
            if line[:6] == "HETATM" and line[17:17+3] == "MSE":
                line = line.replace("HETATM", "ATOM  ")
                line = line.replace("MSE", "MET")
            
            if line[:4] == "ATOM":
                ch = line[21:22]
                if ch == chain or chain is None:
                    atom = line[12:12+4].strip()
                    resi = line[17:17+3]
                    resn = line[22:22+5].strip()
                    x, y, z = [float(line[i:(i+8)]) for i in [30, 38, 46]]
                    
                    if resn[-1].isalpha():
                        resa, resn = resn[-1], int(resn[:-1]) - 1
                    else:
                        resa, resn = "", int(resn) - 1
                    
                    if resn < min_resn:
                        min_resn = resn
                    if resn > max_resn:
                        max_resn = resn
                    if resn not in xyz:
                        xyz[resn] = {}
                    if resa not in xyz[resn]:
                        xyz[resn][resa] = {}
                    if resn not in seq:
                        seq[resn] = {}
                    if resa not in seq[resn]:
                        seq[resn][resa] = resi
                    
                    if atom not in xyz[resn][resa]:
                        xyz[resn][resa][atom] = np.array([x, y, z])
    
    # Convert to numpy arrays
    seq_, xyz_ = [], []
    if min_resn > max_resn:
        return 'no_chain', 'no_chain'
    
    for resn in range(int(min_resn), int(max_resn) + 1):
        if resn in seq:
            for k in sorted(seq[resn]):
                seq_.append(AA_3_N.get(seq[resn][k], 20))
        else:
            seq_.append(20)
        if resn in xyz:
            for k in sorted(xyz[resn]):
                for atom in atoms:
                    if atom in xyz[resn][k]:
                        xyz_.append(xyz[resn][k][atom])
                    else:
                        xyz_.append(np.full(3, np.nan))
        else:
            for atom in atoms:
                xyz_.append(np.full(3, np.nan))
    
    return np.array(xyz_).reshape(-1, len(atoms), 3), N_to_AA(np.array(seq_))


def parse_yaml(yaml_path):
    """
    Parse YAML file to extract chain info and X positions.
    
    Returns:
        chain_info: {chain_id: {'sequence': str, 'x_positions': list of 1-indexed positions}}
        chain_order: list of chain IDs in YAML order
    """
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    
    chain_info = {}
    chain_order = []  # Keep track of chain order in YAML
    for seq_entry in data.get('sequences', []):
        if 'protein' in seq_entry:
            protein = seq_entry['protein']
            chain_id = protein['id']
            sequence = protein['sequence']
            
            # Find positions of 'X' (1-indexed for AntiBMPNN)
            x_positions = [i + 1 for i, aa in enumerate(sequence) if aa == 'X']
            
            chain_info[chain_id] = {
                'sequence': sequence,
                'x_positions': x_positions,
                'has_x': len(x_positions) > 0
            }
            chain_order.append(chain_id)
    
    return chain_info, chain_order


def parse_pdb_for_antibmpnn(pdb_path):
    """
    Parse a PDB file into the format required by AntiBMPNN.
    
    Returns:
        dict: parsed PDB dictionary
    """
    chain_alphabet = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") + [str(i) for i in range(300)]
    
    my_dict = {}
    concat_seq = ''
    s = 0
    
    for letter in chain_alphabet:
        sidechain_atoms = ['N', 'CA', 'C', 'O']
        xyz, seq = parse_PDB_biounits(pdb_path, atoms=sidechain_atoms, chain=letter)
        
        if type(xyz) != str:
            concat_seq += seq[0]
            my_dict['seq_chain_' + letter] = seq[0]
            coords_dict_chain = {
                'N_chain_' + letter: xyz[:, 0, :].tolist(),
                'CA_chain_' + letter: xyz[:, 1, :].tolist(),
                'C_chain_' + letter: xyz[:, 2, :].tolist(),
                'O_chain_' + letter: xyz[:, 3, :].tolist()
            }
            my_dict['coords_chain_' + letter] = coords_dict_chain
            s += 1
    
    # Extract name from path
    pdb_name = Path(pdb_path).stem
    my_dict['name'] = pdb_name
    my_dict['num_of_chains'] = s
    my_dict['seq'] = concat_seq
    
    return my_dict


def main():
    parser = argparse.ArgumentParser(
        description='Prepare AntiBMPNN input from predicted structures and YAML files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--yaml_dir', type=str, required=True,
                        help='Directory containing YAML files with sequence info')
    parser.add_argument('--pdb_dir', type=str, required=True,
                        help='Directory containing predicted PDB files (e.g., boltz_results_xxx/predictions/)')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for prepared input files')
    parser.add_argument('--pdb_suffix', type=str, default='_model_0.pdb',
                        help='Suffix for PDB files (default: _model_0.pdb for non-relaxed)')
    
    args = parser.parse_args()
    
    # Create output directories
    output_dir = Path(args.output_dir)
    pdb_input_dir = output_dir / 'pdbs'
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_input_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all YAML files
    yaml_files = glob.glob(os.path.join(args.yaml_dir, '*.yaml'))
    print(f"Found {len(yaml_files)} YAML files")
    
    # Process each YAML file
    parsed_pdbs = []
    assigned_chains = {}
    fixed_positions = {}
    
    processed_count = 0
    skipped_count = 0
    
    for yaml_path in yaml_files:
        yaml_name = Path(yaml_path).stem  # e.g., "7vgs_D_C_A"
        
        # Find corresponding PDB file
        pdb_subdir = os.path.join(args.pdb_dir, yaml_name)
        pdb_path = os.path.join(pdb_subdir, f"{yaml_name}{args.pdb_suffix}")
        
        if not os.path.exists(pdb_path):
            print(f"Warning: PDB file not found: {pdb_path}")
            skipped_count += 1
            continue
        
        # Parse YAML to get chain info and X positions
        chain_info, yaml_chain_order = parse_yaml(yaml_path)
        
        # Find YAML chains that need to be designed (have X in sequence)
        yaml_chains_to_design = [ch for ch in yaml_chain_order if chain_info[ch]['has_x']]
        
        if not yaml_chains_to_design:
            print(f"Warning: No chains with 'X' positions in {yaml_name}, skipping")
            skipped_count += 1
            continue
        
        # Copy PDB to input directory
        dest_pdb_path = pdb_input_dir / f"{yaml_name}.pdb"
        shutil.copy2(pdb_path, dest_pdb_path)
        
        # Parse PDB for AntiBMPNN
        parsed_dict = parse_pdb_for_antibmpnn(str(dest_pdb_path))
        # Use yaml_name as the name to match the copied PDB
        parsed_dict['name'] = yaml_name
        parsed_pdbs.append(parsed_dict)
        
        # Get all chains from parsed PDB (in order)
        pdb_chains = sorted([k.split('_')[-1] for k in parsed_dict.keys() if k.startswith('seq_chain_')])
        
        # Build mapping: YAML chain ID -> PDB chain ID (based on order)
        # YAML order: [E, D, A] -> PDB order: [A, B, C] => E->A, D->B, A->C
        yaml_to_pdb = {}
        for i, yaml_ch in enumerate(yaml_chain_order):
            if i < len(pdb_chains):
                yaml_to_pdb[yaml_ch] = pdb_chains[i]
        
        # Map chains to design from YAML IDs to PDB IDs
        chains_to_design = [yaml_to_pdb[ch] for ch in yaml_chains_to_design if ch in yaml_to_pdb]
        fixed_chains = [ch for ch in pdb_chains if ch not in chains_to_design]
        
        # Assigned chains format: {"name": [[designed], [fixed]]}
        assigned_chains[yaml_name] = [chains_to_design, fixed_chains]
        
        # Fixed positions format: {"name": {"chain": [fixed_positions_list]}}
        fixed_pos_dict = {}
        for pdb_chain in pdb_chains:
            seq_key = f'seq_chain_{pdb_chain}'
            if seq_key in parsed_dict:
                seq_length = len(parsed_dict[seq_key])
                all_positions = list(range(1, seq_length + 1))  # 1-indexed
                
                # Find corresponding YAML chain
                yaml_chain = None
                for y_ch, p_ch in yaml_to_pdb.items():
                    if p_ch == pdb_chain:
                        yaml_chain = y_ch
                        break
                
                if pdb_chain in chains_to_design and yaml_chain and yaml_chain in chain_info:
                    # For design chains: fix all positions EXCEPT X positions
                    x_positions = chain_info[yaml_chain]['x_positions']
                    fixed_pos = [p for p in all_positions if p not in x_positions]
                    fixed_pos_dict[pdb_chain] = fixed_pos
                else:
                    # For non-design chains: fix all positions
                    fixed_pos_dict[pdb_chain] = all_positions
        
        fixed_positions[yaml_name] = fixed_pos_dict
        processed_count += 1
        
        # Print summary for this entry
        design_summary = []
        for yaml_ch in yaml_chains_to_design:
            pdb_ch = yaml_to_pdb.get(yaml_ch)
            if yaml_ch in chain_info and pdb_ch:
                x_pos = chain_info[yaml_ch]['x_positions']
                if x_pos:
                    design_summary.append(f"{yaml_ch}->{pdb_ch}:[{min(x_pos)}-{max(x_pos)}]({len(x_pos)} pos)")
        print(f"  {yaml_name}: design {', '.join(design_summary)}")
    
    # Write output JSONL files
    parsed_path = output_dir / 'parsed_pdbs.jsonl'
    assigned_path = output_dir / 'assigned_pdbs.jsonl'
    fixed_path = output_dir / 'fixed_pdbs.jsonl'
    
    with open(parsed_path, 'w') as f:
        for entry in parsed_pdbs:
            f.write(json.dumps(entry) + '\n')
    
    with open(assigned_path, 'w') as f:
        f.write(json.dumps(assigned_chains) + '\n')
    
    with open(fixed_path, 'w') as f:
        f.write(json.dumps(fixed_positions) + '\n')
    
    print(f"\n{'='*50}")
    print(f"Processed: {processed_count} entries")
    print(f"Skipped: {skipped_count} entries")
    print(f"\nOutput files created in {output_dir}:")
    print(f"  - parsed_pdbs.jsonl")
    print(f"  - assigned_pdbs.jsonl")
    print(f"  - fixed_pdbs.jsonl")
    print(f"  - pdbs/ (directory with copied PDB files)")


if __name__ == "__main__":
    main()

