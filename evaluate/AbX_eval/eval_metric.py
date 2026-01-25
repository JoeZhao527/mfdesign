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

from ruamel.yaml import YAML

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

    # print(ref_pdb_files_list)
    # print(ref_pdb_name_list)
    for ref_pdb_fpath, pdb_name in zip(ref_pdb_files_list, ref_pdb_name_list):
        # if 'relaxed' not in ref_pdb:
        # pdb_name = ((ref_pdb.split('/')[-1]).split('.pdb'))[0]
        # print(pdb_name)
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
            try:
                cdr_def, _, _ = boltz_cdr_numbering(ref_heavy_str, ref_light_str, heavy_cdr_mask, light_cdr_mask)
                # print(cdr_def)
            except AssertionError:
                print(f'The pdb {pdb_name} have problem')
        elif args.model == 'dymean':
            try:
                cdr_def = cdr_numbering(ref_heavy_str, ref_light_str)
            except AssertionError:
                print(f'The pdb {pdb_name} have problem')
                
        data = {
            'cdr_def': cdr_def,
            'coords': ref_ab_ca,
            'str_seq': ref_ab_str_seq,
            'heavy_str': ref_heavy_str,
            'light_str': ref_light_str,
        }
        reference_data[f'{pdb_name}'] = data
    # print(f"ref: {reference_data}")
    func = functools.partial(boltz_eval_metric, args=args, reference_data=reference_data)
    with mp.Pool(args.cpus) as p:
        results = p.starmap(func, ((pdb_file, pdb_name) for pdb_file, pdb_name,
                                   in zip(predict_pdb_fpath_list, pred_pdb_name_list)))
    # results = []
    # # print(pred_pdb_name_list)
    # print(len(predict_pdb_fpath_list))
    # for pdb_file, pdb_name in tqdm(zip(predict_pdb_fpath_list, pred_pdb_name_list), total=len(predict_pdb_fpath_list)):
    #     # print(pdb_name)
    #     # print(pdb_file)
    #     # try:
    #     result = boltz_eval_metric(pdb_file, pdb_name, args=args, reference_data=reference_data)
    #     results.append(result)
    #     # except:
    #     #     continue

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
                        default='boltz',   # diffab, dymean, boltz, iggm, rfab
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


