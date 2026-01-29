import os
import argparse
import functools
import multiprocessing as mp
import json
import pandas as pd
import numpy as np
import re
import csv
from tqdm import tqdm
from abx.metric import (boltz_eval_metric, 
                        boltz_make_coords, 
                        boltz_cdr_numbering, 
                        InterfaceEnergy, 
                        cdr_numbering)
import glob

from ruamel.yaml import YAML
import warnings
warnings.filterwarnings("ignore")

def read_yaml_file(file_path, heavy_id, light_id):
    """
    Reads and loads a YAML file using ruamel.yaml.

    Parameters:
    - file_path (str): Path to the YAML file.

    Returns:
    - dict: Parsed content of the YAML file.
    """
    yaml = YAML()
    with open(file_path, 'r') as file:
        data = yaml.load(file)
    seq_data = data['sequences']
    heavy_spec_mask = ''
    light_spec_mask = ''
    for pro in seq_data:
        pro_data = pro['protein']
        if pro_data['id'] == heavy_id:
            heavy_spec_mask = pro_data['spec_mask']
            # heavy_seq = pro_data['ground_truth']
        elif pro_data['id'] == light_id:
            light_spec_mask = pro_data['spec_mask']
            # light_seq = pro_data['ground_truth']
    return heavy_spec_mask, light_spec_mask # heavy_seq, light_seq


def parse_list(data_dir):
    input_fname_pattern = '\.pdb$'
    relax_fname_pattern = '\_relaxed.pdb$'
    visited = set()
    for parent, _, files in os.walk(data_dir):
        for fname in files:
            # print(f"fname: {fname}")
            fpath = os.path.join(parent, fname)
            if not re.search(input_fname_pattern, fname):
                continue
            if re.search(relax_fname_pattern, fname):
                continue
            if os.path.getsize(fpath) == 0:
                continue
            if fpath in visited:
                continue
            visited.add(fpath)

            yield fpath


def iggm_extract_base_name(file_name):
    # Remove the extension
    base_name = os.path.splitext(file_name)[0]
    # Find the portion before "_All"
    base_name = base_name.split('_All')[0]
    # Remove "_NA" if present
    base_name = base_name.replace('NA', '')
    return base_name


def get_iggm_files_list(data_dir):
    pdb_files_path_list = []
    pdb_files_name_list = []    
    for pdb_file in os.listdir(data_dir):
        if pdb_file.endswith('.pdb'):
            pdb_file_path = os.path.join(data_dir, pdb_file)
            pdb_files_path_list.append(pdb_file_path)
            
            base_name = iggm_extract_base_name(pdb_file)
            pdb_files_name_list.append(base_name)
            
    return pdb_files_path_list, pdb_files_name_list


def get_rfab_files_list(data_dir, suffix="_relaxed.pdb"):
    pdb_files_path_list = []        
    pdb_files_name_list = []
    pdb_fname_list =[f for f in os.listdir(data_dir) if f.endswith('.pdb')]
    for fname in pdb_fname_list:
        if fname.endswith(suffix):
            fname_list = fname.split('_')
            pdb_name = '_'.join([fname_list[0], fname_list[1], fname_list[2], fname_list[3]])
            fpath = os.path.join(data_dir, fname)
            pdb_files_path_list.append(fpath)
            pdb_files_name_list.append(pdb_name)
    return pdb_files_path_list, pdb_files_name_list


def get_diffab_files_list(data_dir, suffix="_reference.pdb"):
    pdb_files_path_list = []
    pdb_files_name_list = []
    suffix_length = 21
    # For diffab
    for gen_dir in os.listdir(data_dir):
        if 'log' not in gen_dir and '.csv' not in gen_dir and '.json' not in gen_dir:
            pdb_name = gen_dir[:-suffix_length]
            pdb_dir = os.path.join(args.data_dir, f'{gen_dir}/MultipleCDRs/')
            # pdb_dir = os.path.join(args.data_dir, f'{gen_dir}/H_CDR3/')
            count_relax_pdb = 0
            # Only consider about the pdb.
            list_pdb_dir = [pdb for pdb in os.listdir(pdb_dir) if 
                            pdb.endswith('_filter_relaxed.pdb')]
            # list_pdb_dir = [pdb for pdb in os.listdir(pdb_dir) if 
            #                 not pdb.endswith('filter.pdb') and not pdb.endswith('relaxed.pdb')]
            # list_pdb_dir = [pdb for pdb in os.listdir(pdb_dir) if pdb.endswith('relaxed.pdb')]                
            # print(list_pdb_dir)
            for pdb in sorted(list_pdb_dir):
                # pdb = '0001_filter_relaxed.pdb'
                pdb_file_path = os.path.join(pdb_dir, pdb)
                pdb_files_path_list.append(pdb_file_path)
                pdb_files_name_list.append(pdb_name)
                count_relax_pdb += 1
                if count_relax_pdb >= 20:
                    break
                
    return pdb_files_path_list, pdb_files_name_list 


def get_pdb_files(directory):
    return [f for f in os.listdir(directory) if (not f.endswith('json') 
                                                 and not f.endswith('log')
                                                 and not f.endswith('csv'))]


def get_dymean_files_list(data_dir,):
    pdb_files_path_list = []
    pdb_files_name_list = []
    pdb_name_list = get_pdb_files(data_dir)
    for pdb_name in tqdm(pdb_name_list):
        pdb_dir = os.path.join(data_dir, pdb_name)
        # relaxed_pdb_fpath = os.path.join(pdb_dir, f'{pdb_name}_filter_relaxed.pdb')
        relaxed_pdb_fpath = os.path.join(pdb_dir, f'{pdb_name}_openmm_relaxed.pdb')
        
        pdb_files_path_list.append(relaxed_pdb_fpath)
        pdb_files_name_list.append(pdb_name)
    return pdb_files_path_list, pdb_files_name_list


def get_ref_files_list(dir_path, suffix="_reference.pdb"):
    ref_files_path_list = []
    pdb_files_name_list = []
    for filename in os.listdir(dir_path):
        # print(filename)
        if filename.endswith(suffix):
            fpath = os.path.join(dir_path, filename)
            ref_files_path_list.append(fpath)

            file_name = filename[: -len(suffix)]
            pdb_files_name_list.append(file_name)
            
    return ref_files_path_list, pdb_files_name_list


def get_dymean_ref_files_list(dir_path, suffix="_reference.pdb"):
    ref_files_path_list = []
    pdb_files_name_list = []
    pdb_name_list = get_pdb_files(dir_path)
    for filename in pdb_name_list:
        if not filename.endswith('.csv'):
            pdb_files_name_list.append(filename)
            ref_file_path = os.path.join(dir_path, filename + "/" + f'{filename}_original.pdb')
            ref_files_path_list.append(ref_file_path)
            
    return ref_files_path_list, pdb_files_name_list

def get_boltz_files_list(data_dir, suffix="_relaxed.pdb"):
    pdb_files_path_list = []
    pdb_files_name_list = []
    pdb_name_list =[pdb_dir for pdb_dir in os.listdir(data_dir) if (not pdb_dir.endswith('.log')
                                                                    and not pdb_dir.endswith('.csv'))]
    for pdb_dir_name in pdb_name_list:
        pdb_dir = os.path.join(data_dir, pdb_dir_name)
        count_relax_pdb = 0
        for pdb_fname in os.listdir(pdb_dir):
            if pdb_fname.endswith(suffix):
                # When suffix is ".pdb", exclude "_relaxed.pdb" files
                if suffix == ".pdb" and pdb_fname.endswith("_relaxed.pdb"):
                    continue
                pdb_fpath = os.path.join(pdb_dir, pdb_fname)
                pdb_files_path_list.append(pdb_fpath)
                pdb_files_name_list.append(pdb_dir_name)
                count_relax_pdb += 1
                if count_relax_pdb >= 20:
                    break
    # pdb_fname_list =[f for f in os.listdir(data_dir) if not f.endswith('.log') and not f.endswith('.csv')]
    # for fname in pdb_fname_list:
    #     if fname.endswith(suffix):
    #         fname_list = fname.split('_')
    #         pdb_name = '_'.join([fname_list[0], fname_list[1], fname_list[2], fname_list[3]])
    #         fpath = os.path.join(data_dir, fname)
    #         pdb_files_path_list.append(fpath)
    #         pdb_files_name_list.append(pdb_name)
    
    return pdb_files_path_list, pdb_files_name_list


def parse_antibmpnn_fasta(fasta_path):
    """
    Parse AntiBMPNN output FASTA file.
    Returns: (pdb_name, gt_heavy, gt_light, list of (pred_heavy, pred_light, seq_recovery))
    """
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    
    pdb_name = None
    gt_heavy, gt_light = None, None
    predictions = []
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('>'):
            # Parse header
            if pdb_name is None:
                # First header is ground truth
                # Format: >8r9y_H_L_A, score=2.4788, ...
                pdb_name = line.split(',')[0][1:]  # Remove '>' and get name
                i += 1
                seq_line = lines[i].strip()
                parts = seq_line.split('/')
                gt_heavy = parts[0] if len(parts) > 0 else ''
                gt_light = parts[1] if len(parts) > 1 else ''
            else:
                # Subsequent headers are predictions
                # Format: >T=0.1, sample=1, score=0.4146, global_score=1.0845, seq_recovery=0.3846
                seq_recovery = 0.0
                for part in line.split(','):
                    part = part.strip()
                    if 'seq_recovery=' in part:
                        seq_recovery = float(part.split('=')[1])
                        break
                i += 1
                seq_line = lines[i].strip()
                parts = seq_line.split('/')
                pred_heavy = parts[0] if len(parts) > 0 else ''
                pred_light = parts[1] if len(parts) > 1 else ''
                predictions.append((pred_heavy, pred_light, seq_recovery))
        i += 1
    
    return pdb_name, gt_heavy, gt_light, predictions


def get_antibmpnn_files_list(data_dir):
    """Get list of AntiBMPNN FASTA output files."""
    fasta_files = glob.glob(os.path.join(data_dir, '*.fa'))
    fasta_path_list = []
    pdb_name_list = []
    
    for fpath in fasta_files:
        # Extract pdb_name from filename: "8r9y_H_L_A|antibmpnn_000|N-0.0|T-0.1.fa"
        fname = os.path.basename(fpath)
        pdb_name = fname.split('|')[0]
        fasta_path_list.append(fpath)
        pdb_name_list.append(pdb_name)
    
    return fasta_path_list, pdb_name_list


def calc_seq_aar(gt_seq, pred_seq, cdr_def, scheme='chothia'):
    """
    Calculate AAR from sequences without coordinates.
    Similar to calc_ab_metrics but only computes AAR.
    """
    from collections import OrderedDict
    
    ret = OrderedDict()
    ret.update({'full_len': len(gt_seq)})
    
    _schema = {'cdr1': 1, 'cdr2': 3, 'cdr3': 5}
    cdr_idx = {v: 'heavy_' + k for k, v in _schema.items()}
    cdr_idx.update({v + 7: 'light_' + k for k, v in _schema.items()})
    
    for k, v in cdr_idx.items():
        indices = (cdr_def == k)
        gt_s = ''.join([char for char, keep in zip(gt_seq, indices) if keep])
        pred_s = ''.join([char for char, keep in zip(pred_seq, indices) if keep])
        
        if len(gt_s) > 0 and len(pred_s) > 0:
            AAR = np.mean([a == b for a, b in zip(gt_s, pred_s)])
            ret.update({v + '_AAR': AAR})
            
            # For CDR3, also compute Loop AAR
            if k == 5:  # heavy_cdr3
                if scheme == 'imgt' and len(gt_s) > 6:
                    AAR_loop = np.mean([a == b for a, b in zip(gt_s[4:-2], pred_s[4:-2])])
                elif scheme == 'chothia' and len(gt_s) > 4:
                    AAR_loop = np.mean([a == b for a, b in zip(gt_s[2:-2], pred_s[2:-2])])
                else:
                    AAR_loop = AAR
                ret.update({v + '_Loop_AAR': AAR_loop})
            elif k == 12:  # light_cdr3
                if scheme == 'imgt' and len(gt_s) > 6:
                    AAR_loop = np.mean([a == b for a, b in zip(gt_s[4:-2], pred_s[4:-2])])
                elif scheme == 'chothia' and len(gt_s) > 4:
                    AAR_loop = np.mean([a == b for a, b in zip(gt_s[2:-2], pred_s[2:-2])])
                else:
                    AAR_loop = AAR
                ret.update({v + '_Loop_AAR': AAR_loop})
    
    return ret


def antibmpnn_eval_aar(fasta_path, pdb_name, reference_data, args):
    """
    Evaluate AAR for AntiBMPNN output.
    """
    _, gt_heavy, gt_light, predictions = parse_antibmpnn_fasta(fasta_path)
    
    if pdb_name not in reference_data:
        print(f"Warning: {pdb_name} not in reference_data")
        return None
    
    ref_data = reference_data[pdb_name]
    cdr_def = ref_data['cdr_def']
    
    results = []
    for pred_heavy, pred_light, seq_recovery in predictions:
        gt_seq = gt_heavy + gt_light
        pred_seq = pred_heavy + pred_light
        
        if len(gt_seq) != len(pred_seq) or len(gt_seq) != len(cdr_def):
            continue
        
        metrics = calc_seq_aar(gt_seq, pred_seq, cdr_def, scheme=args.scheme)
        metrics.update({
            'code': pdb_name,
            'file_path': fasta_path,
            'mpnn_seq_recovery': seq_recovery
        })
        results.append(metrics)
    
    return results

def main(args):
    reference_data = {}
    results = []
    
    if args.model == 'diffab':
        predict_pdb_fpath_list, pred_pdb_name_list = get_diffab_files_list(args.data_dir)
        ref_pdb_files_list, ref_pdb_name_list = get_ref_files_list(args.ref_dir)
    elif args.model == 'dymean':
        predict_pdb_fpath_list, pred_pdb_name_list = get_dymean_files_list(args.data_dir)
        ref_pdb_files_list, ref_pdb_name_list = get_dymean_ref_files_list(args.dymean_ref_dir)
    elif args.model == 'boltz':
        predict_pdb_fpath_list, pred_pdb_name_list = get_boltz_files_list(args.data_dir, suffix=args.suffix)
        ref_pdb_files_list, ref_pdb_name_list = get_ref_files_list(args.ref_dir)
    elif args.model == 'iggm':
        predict_pdb_fpath_list, pred_pdb_name_list = get_iggm_files_list(args.data_dir)
        ref_pdb_files_list, ref_pdb_name_list = get_ref_files_list(args.ref_dir)
    elif args.model == 'rfab':
        predict_pdb_fpath_list, pred_pdb_name_list = get_rfab_files_list(args.data_dir)
        ref_pdb_files_list, ref_pdb_name_list = get_ref_files_list(args.ref_dir)
    elif args.model == 'antibmpnn':
        # AntiBMPNN: read FASTA files, no need for ref PDB files
        predict_pdb_fpath_list, pred_pdb_name_list = get_antibmpnn_files_list(args.data_dir)
        ref_pdb_files_list, ref_pdb_name_list = [], []
    print(f"[DEBUG] Found {len(predict_pdb_fpath_list)} prediction files")
    print(f"[DEBUG] Found {len(ref_pdb_files_list)} reference files")
    if len(predict_pdb_fpath_list) > 0:
        print(f"[DEBUG] First pred file: {predict_pdb_fpath_list[0]}")
        print(f"[DEBUG] First pred name: {pred_pdb_name_list[0]}")
    
    # # Filtering for antibody or nanobody
    with open(args.test_json_fpath, 'r') as f:
        test_pdb_name_list = json.load(f)  
    
    print(f"[DEBUG] test_json has {len(test_pdb_name_list)} entries")
    if len(test_pdb_name_list) > 0:
        print(f"[DEBUG] First test entry: {test_pdb_name_list[0]}")
        
    test_pdb_name_set = set(test_pdb_name_list)
    filtered_predict = [(f, n) for f, n in zip(predict_pdb_fpath_list, pred_pdb_name_list) if n in test_pdb_name_set]
    print(f"[DEBUG] After filtering: {len(filtered_predict)} files match")
    predict_pdb_fpath_list, pred_pdb_name_list = zip(*filtered_predict) if filtered_predict else ([], [])

    # Build reference_data
    if args.model == 'antibmpnn':
        # For antibmpnn, build reference_data from YAML files only (ground truth is in FASTA)
        for fasta_path, pdb_name in zip(predict_pdb_fpath_list, pred_pdb_name_list):
            if pdb_name in reference_data:
                continue
            code, heavy_chain_id, light_chain_id, antigen_chain = pdb_name.split('_')
            yaml_fpath = os.path.join(args.test_yaml_dir, f'{pdb_name}.yaml')
            if not os.path.exists(yaml_fpath):
                print(f"Warning: YAML file not found: {yaml_fpath}")
                continue
            # Parse FASTA to get ground truth sequences for length
            _, gt_heavy, gt_light, _ = parse_antibmpnn_fasta(fasta_path)
            heavy_cdr_mask, light_cdr_mask = read_yaml_file(yaml_fpath, heavy_chain_id, light_chain_id)
            cdr_def, _, _ = boltz_cdr_numbering(gt_heavy, gt_light, heavy_cdr_mask, light_cdr_mask)
            data = {
                'cdr_def': cdr_def,
                'heavy_str': gt_heavy,
                'light_str': gt_light,
            }
            reference_data[pdb_name] = data
    else:
        # Original logic for other models
        for ref_pdb_fpath, pdb_name in zip(ref_pdb_files_list, ref_pdb_name_list):
            code, heavy_chain_id, light_chain_id, antigen_chain = pdb_name.split('_')
            antibody_id_list = [heavy_chain_id, light_chain_id]
            ref_ab_ca, ref_ab_str_seq, ref_heavy_str, ref_light_str = boltz_make_coords(ref_pdb_fpath, antibody_id_list)
            # Using that our cdr definition.
            if (args.model == 'diffab' 
                or args.model == 'boltz' 
                or args.model == 'iggm' 
                or args.model == 'rfab'):
                yaml_fpath = os.path.join(args.test_yaml_dir, f'{pdb_name}.yaml')
                heavy_cdr_mask, light_cdr_mask = read_yaml_file(yaml_fpath, heavy_chain_id, light_chain_id)
                cdr_def, _, _ = boltz_cdr_numbering(ref_heavy_str, ref_light_str, heavy_cdr_mask, light_cdr_mask)
            elif args.model == 'dymean':
                cdr_def = cdr_numbering(ref_heavy_str, ref_light_str)
                    
            data = {
                'cdr_def': cdr_def,
                'coords': ref_ab_ca,
                'str_seq': ref_ab_str_seq,
                'heavy_str': ref_heavy_str,
                'light_str': ref_light_str,
            }
            reference_data[f'{pdb_name}'] = data
    # Run evaluation
    if args.model == 'antibmpnn':
        # For antibmpnn, each FASTA file contains multiple samples
        results = []
        for fasta_path, pdb_name in tqdm(zip(predict_pdb_fpath_list, pred_pdb_name_list), 
                                          total=len(predict_pdb_fpath_list), desc="Evaluating"):
            file_results = antibmpnn_eval_aar(fasta_path, pdb_name, reference_data, args)
            if file_results:
                results.extend(file_results)
    else:
        func = functools.partial(boltz_eval_metric, args=args, reference_data=reference_data)
        with mp.Pool(args.cpus) as p:
            results = p.starmap(func, ((pdb_file, pdb_name) for pdb_file, pdb_name,
                                       in zip(predict_pdb_fpath_list, pred_pdb_name_list)))

    # Average Results for each Metric
    print(f"[DEBUG] Raw results count: {len(results)}")
    none_count = sum(1 for r in results if r is None)
    print(f"[DEBUG] None results: {none_count}")
    results = [r for r in results if r is not None]
    print(f"[DEBUG] Valid results: {len(results)}")
    if len(results) > 0:
        print(f"[DEBUG] First result keys: {results[0].keys()}")
    df = pd.DataFrame(results)
    print(df)
    
    grouped = df.groupby('code').mean(numeric_only=True)

    
    column_means = grouped.mean(numeric_only=True)
    column_sems = grouped.sem(numeric_only=True)  # Standard Error of the Mean (SEM)

    if args.model == 'antibmpnn':
        # Only AAR metrics for antibmpnn
        filtered_means = column_means.filter(like='AAR')
        filtered_sems = column_sems.filter(like='AAR')
    else:
        filtered_means = pd.concat([
            column_means.filter(like='RMSD'),
            column_means.filter(like='AAR')
        ])
        filtered_sems = pd.concat([
            column_sems.filter(like='RMSD'),
            column_sems.filter(like='AAR')
        ])
    print(f"---------------------")
    print(f"Average Results for each Metric (mean ± SEM, grouped by code)")
    print(f"---------------------")
    for metric in filtered_means.index:
        mean = filtered_means[metric]
        sem = filtered_sems[metric]
        print(f"{metric}: {mean:.4f} ± {sem:.4f}")

    csv_file_path = os.path.join(args.data_dir, 'results.csv')
    with open(csv_file_path, mode='w', newline='') as file:
        fieldnames = results[0].keys() if results else []
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cpus', type=int, default=1)
    parser.add_argument('-e', '--energy', type=bool, default=False)
    parser.add_argument('-v', '--verbose', type=bool, default=False)
    parser.add_argument('-test_json_fpath', type=str,
                        # default='/mnt/nas-new/home/yangnianzu/icml/baseline/AbX/regular_list.json'
                        default='/mnt/nas-new/home/yangnianzu/icml/baseline/AbX/nano_list.json'
                        )
    parser.add_argument('--data_dir', type=str, 
                        default='/mnt/nas-new/home/yangnianzu/icml/baseline/uniform_result'
                        )
    parser.add_argument('--ref_dir', type=str,
                        default='/mnt/nas-new/home/yangnianzu/icml/data/all_structures/test_entry_pdb_files'
                        )
    parser.add_argument('--model', type=str, 
                        default='boltz',   # diffab, dymean, boltz, iggm, rfab, antibmpnn
                        )
    parser.add_argument('--test_yaml_dir', type=str,
                        default='/mnt/nas-new/home/yangnianzu/icml/boltz/src/boltz/data/data_split/test_yml_dir/ab',
                        # default='/mnt/nas-new/home/yangnianzu/icml/boltz/src/boltz/data/data_split/rabd_yml_dir/ab',
                        )
    parser.add_argument('--dymean_ref_dir', type=str,
                        # default='/mnt/nas-new/home/yangnianzu/icml/baseline/dyMEAN/imgt_results',
                        default='/mnt/nas-new/home/yangnianzu/icml/baseline/dyMEAN/out_init_bolt_template_ckpt356_rabd_results'
                        )
    parser.add_argument('--scheme', type=str,
                        default='chothia'
                        )
    parser.add_argument('--suffix', type=str,
                        default='_relaxed.pdb',
                        help='PDB file suffix to match (e.g., "_relaxed.pdb" or ".pdb")'
                        )
    args = parser.parse_args()
    
    main(args)
    # copy_reference_pdb_files()


