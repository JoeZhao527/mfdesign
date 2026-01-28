"""
Prepare Boltz YAML and Protenix input from CDR predictions.

This script:
1. Loads YAML files from boltz_data directory
2. Embeds CDR predictions into sequences based on spec_mask
3. Creates predicted YAML files with embedded CDR sequences
4. Creates YAML files with spec_mask all set to 1
5. Builds Protenix input directories with input.json and cdr_info.json
"""

import json
import os
import re
from pathlib import Path
from typing import Optional


def parse_yaml_filename(filename: str) -> tuple[str, str, Optional[str], str]:
    """Parse YAML filename to extract pdb_id, vh, vl, ag.
    
    Filename format: <pdb_id>_<vh>_<vl>_<ag>.yaml or <pdb_id>_<vh>__<ag>.yaml
    Returns: (pdb_id, vh, vl or None, ag)
    """
    name = filename.replace('.yaml', '')
    parts = name.split('_')
    
    pdb_id = parts[0]
    vh = parts[1]
    vl = parts[2] if parts[2] else None  # Empty string means no light chain
    ag = parts[3] if len(parts) > 3 else parts[2]  # Handle edge cases
    
    # Re-join remaining parts for antigen (in case of multiple chains like AB, CD)
    if len(parts) > 4:
        ag = '_'.join(parts[3:])
    
    return pdb_id, vh, vl, ag


def load_yaml(filepath: Path) -> dict:
    """Load YAML file and return parsed content."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Parse simple YAML format
    result = {'version': 1, 'sequences': []}
    
    lines = content.strip().split('\n')
    current_protein = None
    
    for line in lines:
        if line.strip().startswith('version:'):
            result['version'] = int(line.split(':')[1].strip())
        elif line.strip() == '- protein:':
            current_protein = {}
            result['sequences'].append({'protein': current_protein})
        elif current_protein is not None:
            if line.strip().startswith('id:'):
                current_protein['id'] = line.split(':')[1].strip()
            elif line.strip().startswith('sequence:'):
                current_protein['sequence'] = line.split(':', 1)[1].strip()
            elif line.strip().startswith('spec_mask:'):
                mask = line.split(':', 1)[1].strip()
                # Remove quotes if present
                current_protein['spec_mask'] = mask.strip("'\"")
            elif line.strip().startswith('ground_truth:'):
                current_protein['ground_truth'] = line.split(':', 1)[1].strip()
    
    return result


def save_yaml(data: dict, filepath: Path):
    """Save data to YAML format."""
    lines = [f"version: {data['version']}", "sequences:"]
    
    for seq_item in data['sequences']:
        protein = seq_item['protein']
        lines.append("  - protein:")
        lines.append(f"      id: {protein['id']}")
        lines.append(f"      sequence: {protein['sequence']}")
        lines.append(f"      spec_mask: '{protein['spec_mask']}'")
        if 'ground_truth' in protein:
            lines.append(f"      ground_truth: {protein['ground_truth']}")
    
    with open(filepath, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def find_cdr_regions(spec_mask: str) -> list[tuple[int, int]]:
    """Find start and end indices of CDR regions (where spec_mask is 1)."""
    regions = []
    start = None
    
    for i, c in enumerate(spec_mask):
        if c == '1' and start is None:
            start = i
        elif c == '0' and start is not None:
            regions.append((start, i))
            start = None
    
    if start is not None:
        regions.append((start, len(spec_mask)))
    
    return regions


def embed_cdr_prediction(sequence: str, spec_mask: str, cdr_predictions: list[str]) -> str:
    """Embed CDR predictions into sequence based on spec_mask.
    
    - If predicted CDR is too long, cut from the right
    - If predicted CDR is too short, pad with 'X' on the right
    """
    result = list(sequence)
    regions = find_cdr_regions(spec_mask)
    
    for i, (start, end) in enumerate(regions):
        if i >= len(cdr_predictions) or cdr_predictions[i] is None:
            continue
        
        cdr = cdr_predictions[i]
        region_len = end - start
        
        # Handle length mismatch
        if len(cdr) > region_len:
            # Cut from the right if too long
            cdr = cdr[:region_len]
        elif len(cdr) < region_len:
            # Pad with 'X' on the right if too short
            cdr = cdr + 'X' * (region_len - len(cdr))
        
        # Replace the region
        for j, c in enumerate(cdr):
            result[start + j] = c
    
    return ''.join(result)


def create_protenix_input(yaml_data: dict, yaml_name: str, vh: str, vl: Optional[str], ag: str,
                          output_dir: Path):
    """Create Protenix input files: protenix_json/<name>.json and cdr_info/<name>_cdr_info.json.
    
    Note: The YAML files contain only antibody variable domains. The antigen chain(s)
    are parsed from the filename but their sequences are not in the YAML.
    We include all chains from the YAML in alphabetical order for protenix mapping.
    """
    protenix_json_dir = output_dir / "protenix_json"
    cdr_info_dir = output_dir / "cdr_info"
    protenix_json_dir.mkdir(parents=True, exist_ok=True)
    cdr_info_dir.mkdir(parents=True, exist_ok=True)
    
    # Build chain mapping (protenix uses A, B, C... in order)
    chain_ids = []
    for seq_item in yaml_data['sequences']:
        chain_ids.append(seq_item['protein']['id'])
    
    # Sort to get alphabetical order for original chains
    sorted_chain_ids = sorted(chain_ids)
    
    # Build protenix input.json
    sequences = []
    sequence_order = []
    
    for i, orig_chain in enumerate(sorted_chain_ids):
        protenix_chain = chr(ord('A') + i)
        
        # Find the protein data for this chain
        protein_data = None
        for seq_item in yaml_data['sequences']:
            if seq_item['protein']['id'] == orig_chain:
                protein_data = seq_item['protein']
                break
        
        seq_len = len(protein_data['sequence'])
        
        chain_data = {
            "proteinChain": {
                "sequence": protein_data['sequence'],
                "count": 1,
                "label_entity_id": str(i + 1),
                "label_asym_id": [protenix_chain],
                "auth_asym_id": [orig_chain]
            }
        }
        sequences.append(chain_data)
        
        sequence_order.append({
            "protenix_chain": protenix_chain,
            "original_chain": orig_chain,
            "seq_len": seq_len,
            "entity_type": "proteinChain"
        })
    
    input_json = [{
        "sequences": sequences,
        "name": yaml_name
    }]
    
    # Build cdr_info.json
    protenix_to_original = {}
    original_to_protenix = {}
    
    for i, orig_chain in enumerate(sorted_chain_ids):
        protenix_chain = chr(ord('A') + i)
        protenix_to_original[protenix_chain] = orig_chain
        original_to_protenix[orig_chain] = protenix_chain
    
    cdr_info = {
        "entry_name": yaml_name,
        "heavy_chain_id": vh,
        "light_chain_id": vl,
        "chain_mapping": {
            "protenix_to_original": protenix_to_original,
            "original_to_protenix": original_to_protenix,
            "sequence_order": sequence_order
        },
        "cdr_info": {}
    }
    
    # Extract CDR indices from spec_mask for heavy and light chains
    for seq_item in yaml_data['sequences']:
        protein = seq_item['protein']
        chain_id = protein['id']
        spec_mask = protein['spec_mask']
        
        if chain_id == vh:
            # Heavy chain
            regions = find_cdr_regions(spec_mask)
            all_indices = []
            cdr_data = {"variable_domain_start_res_id": 1}
            
            for i, (start, end) in enumerate(regions):
                indices = list(range(start, end))
                all_indices.extend(indices)
                cdr_data[f"cdr{i+1}"] = {"indices": indices}
            
            cdr_data["cdr_indices"] = sorted(all_indices)
            cdr_info["cdr_info"]["H_chain"] = cdr_data
        
        elif vl and chain_id == vl:
            # Light chain
            regions = find_cdr_regions(spec_mask)
            all_indices = []
            cdr_data = {"variable_domain_start_res_id": 1}
            
            for i, (start, end) in enumerate(regions):
                indices = list(range(start, end))
                all_indices.extend(indices)
                cdr_data[f"cdr{i+1}"] = {"indices": indices}
            
            cdr_data["cdr_indices"] = sorted(all_indices)
            cdr_info["cdr_info"]["L_chain"] = cdr_data
    
    # Save files
    with open(protenix_json_dir / f'{yaml_name}.json', 'w') as f:
        json.dump(input_json, f, indent=2)
    
    with open(cdr_info_dir / f'{yaml_name}_cdr_info.json', 'w') as f:
        json.dump(cdr_info, f, indent=2)


def process_predictions(boltz_data_dir: Path, predictions_dir: Path, output_base_dir: Path):
    """Process all predictions and generate output files."""
    
    # Create protenix outputs using raw input yaml (no predictions)
    protenix_original_mask_dir = output_base_dir / "protenix_original_mask" / "cdr_masking"
    protenix_original_gt_dir = output_base_dir / "protenix_original_gt" / "cdr_masking"
    protenix_original_mask_dir.mkdir(parents=True, exist_ok=True)
    protenix_original_gt_dir.mkdir(parents=True, exist_ok=True)
    
    for yaml_file in boltz_data_dir.glob('*.yaml'):
        yaml_name = yaml_file.stem
        pdb_id, vh, vl, ag = parse_yaml_filename(yaml_file.name)
        yaml_data = load_yaml(yaml_file)
        
        # Create ground truth yaml (replace sequence with ground_truth)
        gt_yaml = {'version': yaml_data['version'], 'sequences': []}
        for seq_item in yaml_data['sequences']:
            protein = seq_item['protein']
            gt_yaml['sequences'].append({
                'protein': {
                    'id': protein['id'],
                    'sequence': protein.get('ground_truth', protein['sequence']),
                    'spec_mask': protein['spec_mask'],
                    'ground_truth': protein.get('ground_truth', protein['sequence'])
                }
            })
        
        create_protenix_input(yaml_data, yaml_name, vh, vl, ag, protenix_original_mask_dir)
        create_protenix_input(gt_yaml, yaml_name, vh, vl, ag, protenix_original_gt_dir)
        
    # Load all prediction files
    all_predictions = {}
    for pred_file in predictions_dir.glob('*.json'):
        pred_name = pred_file.stem  # e.g., 'no_playback_0', 'playback_1'
        with open(pred_file, 'r') as f:
            all_predictions[pred_name] = json.load(f)
    
    # Process each prediction set
    for pred_name, predictions in all_predictions.items():
        print(f"Processing {pred_name}...")
        
        # Create output directories
        boltz_output_dir = output_base_dir / f"boltz_predicted_{pred_name}"
        boltz_all_mask_dir = output_base_dir / f"boltz_all_mask_{pred_name}"
        protenix_output_dir = output_base_dir / f"protenix_{pred_name}" / "cdr_masking"

        boltz_output_dir.mkdir(parents=True, exist_ok=True)
        boltz_all_mask_dir.mkdir(parents=True, exist_ok=True)
        protenix_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Process each YAML file
        for yaml_file in boltz_data_dir.glob('*.yaml'):
            yaml_name = yaml_file.stem
            pdb_id, vh, vl, ag = parse_yaml_filename(yaml_file.name)
            
            # Check if we have predictions for this pdb_id
            if pdb_id not in predictions:
                print(f"  Skipping {yaml_name}: no prediction for {pdb_id}")
                continue
            
            pred = predictions[pdb_id]
            
            # Load original YAML
            yaml_data = load_yaml(yaml_file)
            
            # Create copies for modification
            predicted_yaml = {'version': yaml_data['version'], 'sequences': []}
            all_mask_yaml = {'version': yaml_data['version'], 'sequences': []}
            
            # Parse antigen chain IDs (can be multiple chars like "AB", "CD")
            ag_chain_ids = set(ag) if ag else set()
            
            for seq_item in yaml_data['sequences']:
                protein = seq_item['protein']
                chain_id = protein['id']
                sequence = protein['sequence']
                spec_mask = protein['spec_mask']
                ground_truth = protein.get('ground_truth', sequence)
                
                # Determine if this is an antibody chain or antigen chain
                is_antigen = chain_id in ag_chain_ids
                
                if is_antigen:
                    # Antigen chain: keep sequence as-is
                    new_sequence = sequence
                    # For predicted: all 0s (no masking)
                    predicted_mask = '0' * len(spec_mask)
                    # For all_mask: all 1s
                    all_ones_mask = '1' * len(spec_mask)
                else:
                    # Antibody chain: embed CDR predictions
                    if chain_id == vh:
                        cdr_preds = [pred.get('HCDR1'), pred.get('HCDR2'), pred.get('HCDR3')]
                    elif vl and chain_id == vl:
                        cdr_preds = [pred.get('LCDR1'), pred.get('LCDR2'), pred.get('LCDR3')]
                    else:
                        cdr_preds = []
                    
                    new_sequence = embed_cdr_prediction(sequence, spec_mask, cdr_preds)
                    predicted_mask = spec_mask
                    all_ones_mask = '1' * len(spec_mask)

                # Add to predicted YAML
                predicted_yaml['sequences'].append({
                    'protein': {
                        'id': chain_id,
                        'sequence': new_sequence,
                        'spec_mask': predicted_mask,
                        'ground_truth': ground_truth
                    }
                })
                
                # Add to all-mask YAML
                all_mask_yaml['sequences'].append({
                    'protein': {
                        'id': chain_id,
                        'sequence': new_sequence,
                        'spec_mask': all_ones_mask,
                        'ground_truth': ground_truth
                    }
                })
            
            # Save predicted YAML
            save_yaml(predicted_yaml, boltz_output_dir / f"{yaml_name}.yaml")
            
            # Save all-mask YAML
            save_yaml(all_mask_yaml, boltz_all_mask_dir / f"{yaml_name}.yaml")
            
            # Create Protenix input
            create_protenix_input(predicted_yaml, yaml_name, vh, vl, ag, protenix_output_dir)
        
        print(f"  Completed {pred_name}")


def main():
    # Define paths
    script_dir = Path(__file__).parent
    boltz_data_dir = script_dir / 'boltz_data'
    predictions_dir = script_dir.parent / 'final_result_cdr'
    output_base_dir = script_dir / 'output'
    
    # Process all predictions
    process_predictions(boltz_data_dir, predictions_dir, output_base_dir)
    
    print("Done!")


if __name__ == '__main__':
    main()
