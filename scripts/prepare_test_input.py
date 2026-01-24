"""
输入一个像examples/yaml包含input yaml的目录，输出两个目录：

output_dir/inpaint: 包含把sequence替换为ground truth的yaml
output_dir/fold: 包含把sequence替换为ground truth，并把spec mask全部换为0的的yaml

Usage:
    python scripts/prepare_test_input.py --input_dir examples/yaml --output_dir output
"""

import argparse
import os
from pathlib import Path

import yaml


def process_yaml(input_path: Path, output_dir: Path, mode: str):
    """Process a single yaml file.
    
    Args:
        input_path: Path to input yaml file
        output_dir: Output directory
        mode: 'inpaint' or 'fold'
    """
    with open(input_path, 'r') as f:
        data = yaml.safe_load(f)
    
    # Process each sequence
    for seq_item in data.get('sequences', []):
        # Handle different sequence types (protein, dna, rna, etc.)
        for seq_type, seq_data in seq_item.items():
            if isinstance(seq_data, dict):
                # Replace sequence with ground_truth if available
                if 'ground_truth' in seq_data:
                    seq_data['sequence'] = seq_data['ground_truth']
                
                # For fold mode, set spec_mask to all zeros
                if mode == 'fold' and 'spec_mask' in seq_data:
                    seq_len = len(seq_data['sequence'])
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
        process_yaml(yaml_file, output_dir, 'inpaint')
        process_yaml(yaml_file, output_dir, 'fold')
        print()
    
    print("Done!")
    print(f"  Inpaint YAMLs: {output_dir / 'inpaint'}")
    print(f"  Fold YAMLs: {output_dir / 'fold'}")


if __name__ == '__main__':
    main()
