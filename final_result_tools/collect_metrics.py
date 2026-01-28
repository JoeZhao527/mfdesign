import glob
import json
import numpy as np
import pandas as pd

def load_best_metric_by_entry(file_pattern_base, mode='lddt'):
    """Load metric data grouped by entry and return the best (max) metric for each entry.
    
    Args:
        file_pattern_base: Glob pattern for files
        mode: 'lddt' or 'rmsd'
    
    Returns:
        Tuple of (best_data, per_cdr_data)
        - best_data: Dictionary with metrics as keys and lists of best values as values
        - per_cdr_data: Dictionary with per-CDR metrics
    """
    all_files = glob.glob(file_pattern_base)
    
    entries_data = {}
    
    for file_path in all_files:
        filename = file_path.split('/')[-1]
        
        # Determine file type and extract entry name
        # Support both regular and backbone files (_lddt.json, _lddt_backbone.json, etc.)
        entry = None
        if mode == 'lddt':
            if '_lddt_backbone.json' in filename:
                entry = filename.replace('_lddt_backbone.json', '').rsplit('_sample_', 1)[0]
            elif '_lddt.json' in filename:
                entry = filename.replace('_lddt.json', '').rsplit('_sample_', 1)[0]
        elif mode == 'rmsd':
            if '_rmsd_backbone.json' in filename:
                entry = filename.replace('_rmsd_backbone.json', '').rsplit('_sample_', 1)[0]
            elif '_rmsd.json' in filename:
                entry = filename.replace('_rmsd.json', '').rsplit('_sample_', 1)[0]
        
        if entry is None:
            continue
        
        if entry not in entries_data:
            entries_data[entry] = {
                'total': [], 'cdr': [], 'framework': [],
                'hcdr1': [], 'hcdr2': [], 'hcdr3': [],
                'lcdr1': [], 'lcdr2': [], 'lcdr3': []
            }
        
        with open(file_path, 'r') as f:
            json_data = json.load(f)
            if json_data.get('success') and json_data.get(mode):
                metric_data = json_data[mode]
                entries_data[entry]['total'].append(metric_data['total'])
                entries_data[entry]['cdr'].append(metric_data['cdr'])
                entries_data[entry]['framework'].append(metric_data['framework'])
                
                # Extract per-CDR metrics from heavy and light chains
                heavy_chain = json_data.get('heavy_chain', {})
                light_chain = json_data.get('light_chain', {})

                entries_data[entry]['hcdr1'].append(heavy_chain.get('cdr1', None))
                entries_data[entry]['hcdr2'].append(heavy_chain.get('cdr2', None))
                entries_data[entry]['hcdr3'].append(heavy_chain.get('cdr3', None))
                entries_data[entry]['lcdr1'].append(light_chain.get('cdr1', None))
                entries_data[entry]['lcdr2'].append(light_chain.get('cdr2', None))
                entries_data[entry]['lcdr3'].append(light_chain.get('cdr3', None))
    
    # Compute best metric for each entry
    # For lDDT: higher is better (use max), for RMSD: lower is better (use min)
    best_data = {
        'total': [], 'cdr': [], 'framework': [],
        'hcdr1': [], 'hcdr2': [], 'hcdr3': [],
        'lcdr1': [], 'lcdr2': [], 'lcdr3': []
    }
    
    use_max = (mode == 'lddt')
    
    for entry, data in entries_data.items():
        if len(data['total']) > 0:
            if use_max:
                # Higher is better for lDDT
                best_data['total'].append(np.max(data['total']))
                best_data['cdr'].append(np.max(data['cdr']))
                best_data['framework'].append(np.max(data['framework']))
                for cdr_key in ['hcdr1', 'hcdr2', 'hcdr3', 'lcdr1', 'lcdr2', 'lcdr3']:
                    valid_values = [v for v in data[cdr_key] if v is not None]
                    if valid_values:
                        best_data[cdr_key].append(np.max(valid_values))
            else:  # rmsd
                # Lower is better for RMSD
                best_data['total'].append(np.min(data['total']))
                best_data['cdr'].append(np.min(data['cdr']))
                best_data['framework'].append(np.min(data['framework']))
                for cdr_key in ['hcdr1', 'hcdr2', 'hcdr3', 'lcdr1', 'lcdr2', 'lcdr3']:
                    valid_values = [v for v in data[cdr_key] if v is not None]
                    if valid_values:
                        best_data[cdr_key].append(np.min(valid_values))
    
    return best_data

def plot_metric_comparison(mode='lddt', dirs=None):
    """Plot comparison of metrics across multiple directories.
    
    Args:
        mode: 'lddt' or 'rmsd'
        dirs: List of directory names, or None to use defaults
    """
    if dirs is None:
        dirs = [f'test_{mode}', f'test_gt_{mode}']
    
    # Support both regular and backbone files
    file_suffixes_backbone = f'*sample_*_{mode}_backbone.json'
    
    # Load best metrics from all directories
    all_best_data = {}
    for dir_name in dirs:
        pattern_backbone = f'{dir_name}/{file_suffixes_backbone}'
        best_data = load_best_metric_by_entry(pattern_backbone, mode=mode)
        all_best_data[dir_name] = best_data
        print(f"Loaded best {mode.upper()} for {len(best_data['total'])} entries from {dir_name}")
    
    # Prepare data for violin plots
    metrics = ['total', 'cdr', 'framework']
    metric_labels = [f'Total {mode.upper()}', f'CDR {mode.upper()}', f'Framework {mode.upper()}']
    
    # Print summary statistics
    print(f"\nSummary Statistics ({mode.upper()}):")
    print(f"\nOverall Metrics:")
    for metric_name in metrics:
        print(f"\n{metric_name.capitalize()} {mode.upper()}:")
        for dir_name in dirs:
            best_data = all_best_data[dir_name]
            print(f"  {dir_name}: mean={np.mean(best_data[metric_name]):.4f}, std={np.std(best_data[metric_name]):.4f}")
    
    # Print per-CDR metrics
    print(f"\n\nHeavy Chain CDR Metrics:")
    for cdr_idx in range(1, 4):
        cdr_key = f'hcdr{cdr_idx}'
        # Check if any directory has data for this CDR
        has_data = any(all_best_data[dir_name][cdr_key] for dir_name in dirs)
        if has_data:
            print(f"\nHCDR{cdr_idx}:")
            for dir_name in dirs:
                best_data = all_best_data[dir_name]
                if best_data[cdr_key]:
                    print(f"  {dir_name}: mean={np.mean(best_data[cdr_key]):.4f}, std={np.std(best_data[cdr_key]):.4f}, count={len(best_data[cdr_key])}")
                else:
                    print(f"  {dir_name}: No data available")
    
    print(f"\n\nLight Chain CDR Metrics:")
    for cdr_idx in range(1, 4):
        cdr_key = f'lcdr{cdr_idx}'
        # Check if any directory has data for this CDR
        has_data = any(all_best_data[dir_name][cdr_key] for dir_name in dirs)
        if has_data:
            print(f"\nLCDR{cdr_idx}:")
            for dir_name in dirs:
                best_data = all_best_data[dir_name]
                if best_data[cdr_key]:
                    print(f"  {dir_name}: mean={np.mean(best_data[cdr_key]):.4f}, std={np.std(best_data[cdr_key]):.4f}, count={len(best_data[cdr_key])}")
                else:
                    print(f"  {dir_name}: No data available")
        else:
            print(f"\nLCDR{cdr_idx}: No data available")

if __name__ == '__main__':
    import sys
    
    # Default to lddt mode, can be overridden with command line argument
    mode = 'lddt'
    if len(sys.argv) > 1:
        mode = sys.argv[1].lower()
    
    if mode not in ['lddt', 'rmsd']:
        print(f"Error: mode must be 'lddt' or 'rmsd', got '{mode}'")
        sys.exit(1)
    
    # Get directories from command line arguments (after mode)
    # If not provided, use defaults
    if len(sys.argv) > 2:
        dirs = sys.argv[2:]
    else:
        # Default directories
        dirs = [f'results/playback_0_{mode}', f'results/playback_1_{mode}', f'results/no_playback_0_{mode}', f'results/no_playback_1_{mode}']
        # dirs = [f'test_{mode}', f'test_gt_{mode}']
        # dirs = [f'test_{mode}',  f'train_mask_{mode}', f'test_gt_{mode}', f'train_gt_{mode}']
        # [f'test_random_30_{mode}', f'key_res_{mode}', f'replace_30_{mode}', f'all_noise_{mode}']
    
    # Plot the comparison
    plot_metric_comparison(mode=mode, dirs=dirs)
