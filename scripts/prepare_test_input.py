"""
输入一个像examples/yaml包含input yaml的目录，输出三个目录：

output_dir/inpaint: 包含把sequence替换为ground truth的yaml，抗原spec_mask全为1
output_dir/fold: 包含把sequence替换为ground truth，抗体spec_mask全为1，抗原spec_mask全为1
output_dir/ag_inpaint: 抗原spec_mask全为0，其他跟输入yaml保持一致

Usage:
    python scripts/prepare_test_input.py --input_dir examples/yaml --output_dir output
"""

import argparse
from pathlib import Path
from typing import Optional

import yaml


def parse_yaml_filename(filename: str) -> tuple[str, str, Optional[str], str]:
    """Parse YAML filename to extract pdb_id, vh, vl, ag.
    
    Filename format: <pdb_id>_<vh>_<vl>_<ag>.yaml or <pdb_id>_<vh>__<ag>.yaml
    Returns: (pdb_id, vh, vl or None, ag)
    """
    name = filename.replace('.yaml', '').replace('.yml', '')
    parts = name.split('_')
    
    pdb_id = parts[0]
    vh = parts[1]
    vl = parts[2] if parts[2] else None  # Empty string means no light chain
    ag = parts[3] if len(parts) > 3 else parts[2]  # Handle edge cases
    
    # Re-join remaining parts for antigen (in case of multiple chains like AB, CD)
    if len(parts) > 4:
        ag = '_'.join(parts[3:])
    
    return pdb_id, vh, vl, ag


def get_chain_ids(vh: str, vl: Optional[str], ag: str) -> tuple[set, set]:
    """Get antibody and antigen chain IDs.
    
    Returns: (antibody_chain_ids, antigen_chain_ids)
    """
    antibody_ids = {vh}
    if vl:
        antibody_ids.add(vl)
    
    # Antigen can have multiple chains (e.g., "AB" or "A")
    antigen_ids = set(ag) if ag else set()
    
    return antibody_ids, antigen_ids


def process_yaml(input_path: Path, output_dir: Path, mode: str):
    """Process a single yaml file.
    
    Args:
        input_path: Path to input yaml file
        output_dir: Output directory
        mode: 'inpaint', 'fold', or 'ag_inpaint'
    """
    with open(input_path, 'r') as f:
        data = yaml.safe_load(f)
    
    # Parse filename to get chain IDs
    pdb_id, vh, vl, ag = parse_yaml_filename(input_path.name)
    antibody_ids, antigen_ids = get_chain_ids(vh, vl, ag)
    
    # Process each sequence
    for seq_item in data.get('sequences', []):
        # Handle different sequence types (protein, dna, rna, etc.)
        for seq_type, seq_data in seq_item.items():
            if isinstance(seq_data, dict):
                chain_id = seq_data.get('id', '')
                is_antibody = chain_id in antibody_ids
                is_antigen = chain_id in antigen_ids
                
                if mode == 'inpaint':
                    # Antibody: use ground_truth sequence, keep spec_mask
                    # Antigen: spec_mask all 1s (don't inpaint)
                    if is_antibody and 'ground_truth' in seq_data:
                        seq_data['sequence'] = seq_data['ground_truth']
                    if is_antigen and 'spec_mask' in seq_data:
                        seq_len = len(seq_data.get('sequence', ''))
                        seq_data['spec_mask'] = '1' * seq_len
                        
                elif mode == 'fold':
                    # Antibody: use ground_truth sequence, spec_mask all 1s
                    # Antigen: spec_mask all 1s
                    if is_antibody:
                        if 'ground_truth' in seq_data:
                            seq_data['sequence'] = seq_data['ground_truth']
                        if 'spec_mask' in seq_data:
                            seq_len = len(seq_data.get('sequence', ''))
                            seq_data['spec_mask'] = '1' * seq_len
                    if is_antigen and 'spec_mask' in seq_data:
                        seq_len = len(seq_data.get('sequence', ''))
                        seq_data['spec_mask'] = '1' * seq_len
                        
                elif mode == 'ag_inpaint':
                    # Antigen: spec_mask all 0s (inpaint antigen)
                    # Antibody: keep original (no changes)
                    if is_antigen and 'spec_mask' in seq_data:
                        seq_len = len(seq_data.get('sequence', ''))
                        seq_data['spec_mask'] = '0' * seq_len
    
    # Create output directory if needed
    mode_output_dir = output_dir / mode
    mode_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Write output yaml
    output_path = mode_output_dir / input_path.name
    with open(output_path, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
    
    print(f"  [{mode}] -> {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Prepare test input YAML files for structure prediction')
    parser.add_argument('--input_dir', type=str, required=True,
                        help='Input directory containing YAML files')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for processed YAML files')
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    
    # Find all yaml files
    yaml_files = list(input_dir.glob('*.yaml')) + list(input_dir.glob('*.yml'))
    
    if not yaml_files:
        print(f"No YAML files found in {input_dir}")
        return
    
    print(f"Found {len(yaml_files)} YAML files in {input_dir}")
    print(f"Output directory: {output_dir}")
    print()
    
    for yaml_file in yaml_files:
        print(f"Processing: {yaml_file.name}")
        pdb_id, vh, vl, ag = parse_yaml_filename(yaml_file.name)
        print(f"  Parsed: pdb={pdb_id}, vh={vh}, vl={vl}, ag={ag}")
        process_yaml(yaml_file, output_dir, 'inpaint')
        process_yaml(yaml_file, output_dir, 'fold')
        process_yaml(yaml_file, output_dir, 'ag_inpaint')
        print()
    
    print("Done!")
    print(f"  Inpaint YAMLs: {output_dir / 'inpaint'}")
    print(f"  Fold YAMLs: {output_dir / 'fold'}")
    print(f"  AG Inpaint YAMLs: {output_dir / 'ag_inpaint'}")


if __name__ == '__main__':
    main()
