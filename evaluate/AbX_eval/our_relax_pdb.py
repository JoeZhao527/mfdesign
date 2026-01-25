from abx.single_relax import Rosetta_relax
import os
import argparse
import functools
import multiprocessing as mp
import logging
import itertools
import json
import pandas as pd
import traceback
import pickle
import numpy as np
import pdb
import re
from tqdm import tqdm
import traceback

from ruamel.yaml import YAML
from abx.metric import boltz_eval_metric, boltz_make_coords, boltz_cdr_numbering


def parse_list(data_dir):
    input_fname_pattern = '\.pdb$'
    visited = set()
    for parent, _, files in os.walk(data_dir):
        for fname in files:
            # print(f"fname: {fname}")
            fpath = os.path.join(parent, fname)
            if not re.search(input_fname_pattern, fname):
                continue
            if os.path.getsize(fpath) == 0:
                continue
            if fpath in visited:
                continue
            visited.add(fpath)

            yield fpath


def main(args):
    # fpath = parse_list(args.data_dir)

    func = functools.partial(Rosetta_relax, args=args, ref_data=ref_data, logger=logger, type=args.type)
    # print(f"fpath: {fpath}")

    with mp.Pool(args.cpus) as p:
        p.starmap(func, ((pdb_file, pdb_name, original_pdb_fpath) for pdb_file, pdb_name, original_pdb_fpath in 
                         zip(task_pdb_list, pdb_name_list, original_pdb_list)))
    
    # for pdb_file, pdb_name, original_pdb_fpath in zip(task_pdb_list, pdb_name_list, original_pdb_list):
    #     print(f'Finish {pdb_name}')
    #     print(pdb_file)
    #     print(original_pdb_fpath)
    #     Rosetta_relax(pdb_file, pdb_name, original_pdb_fpath, args=args, ref_data=ref_data, type=args.type, logger=logger)
    # for pdb_file, pdb_name in zip(task_pdb_list[:3], pdb_name_list[:3]):
    #     try:
    #         print(pdb_file)
    #         print(f'True {pdb_name}')
    #         Rosetta_relax(pdb_file, pdb_name, args=args, ref_data=ref_data, type='boltz', logger=logger)
    #     except:
    #         traceback.print_exc()
    #         print(f'Error pdb is {pdb_file}')
            

                    
def get_pdb_files(directory):
    return [f for f in os.listdir(directory) if not f.endswith('json') and not f.endswith('log')]


def get_diffab_pdb_names(input_dir, suffix_length=21):
    list_dir = os.listdir(input_dir)
    pdb_name_list = [adir[:-suffix_length] for adir in list_dir]
    return pdb_name_list


def get_boltz_pdb_names(input_dir):
    pdb_name_list = []
    pdb_fpath_list = []
    
    # Devide to run.
    with open(args.split_json, 'r') as f:
        split_list = json.load(f)
    print(args.split_json)
    print(len(split_list))
    
    list_dir = [dir_name for dir_name in os.listdir(input_dir)
                if not dir_name.endswith('log')]
    for pdb_name in list_dir:
        if pdb_name in split_list:
            print(pdb_name)
            pdb_dir = os.path.join(input_dir, pdb_name)
            pdb_list_dir = os.listdir(pdb_dir)
            for num_pdb in pdb_list_dir:
                if (not num_pdb.endswith('relaxed.pdb')
                and not num_pdb.endswith('.json')
                and not num_pdb.endswith('.npz')
                and not num_pdb.endswith('.seq')
                ):
                    print(num_pdb)
                    pdb_fpath = os.path.join(pdb_dir, num_pdb)
                    pdb_name_list.append(pdb_name)
                    pdb_fpath_list.append(pdb_fpath)

    # pdb_fname_list =[f for f in os.listdir(input_dir) if not f.endswith('relaxed.pdb')
    #                  and not f.endswith('.log') and not f.endswith('.csv')]
    # for fname in pdb_fname_list:
    #     # print(fname)
    #     fname_list = fname.split('_')
    #     pdb_name = '_'.join([fname_list[0], fname_list[1], fname_list[2], fname_list[3]])
    #     if pdb_name in split_list:
    #         fpath = os.path.join(input_dir, fname)
    #         pdb_fpath_list.append(fpath)
    #         pdb_name_list.append(pdb_name)
    #         print(fpath)
    #         print(pdb_name)
    return pdb_fpath_list, pdb_name_list


def get_rfab_pdb_names(input_dir):
    pdb_name_list = []
    pdb_fpath_list = []
    # Devide to run.
    with open(args.split_json, 'r') as f:
        split_list = json.load(f)
    print(args.split_json)
    print(len(split_list))
    pdb_fname_list =[f for f in os.listdir(input_dir) if not f.endswith('relaxed.pdb') 
                 and not f.endswith('.log') and not f.endswith('.csv')]
    for fname in pdb_fname_list:
        # print(fname)
        fname_list = fname.split('_')
        pdb_name = '_'.join([fname_list[0], fname_list[1], fname_list[2], fname_list[3]])
        if pdb_name in split_list:
            fpath = os.path.join(input_dir, fname)
            # print('check')
            # print(pdb_name)
            # print(fpath)
            pdb_fpath_list.append(fpath)
            pdb_name_list.append(pdb_name)
            # print(fpath)
            # print(pdb_name)
    return pdb_fpath_list, pdb_name_list


def get_iggm_pdb_names(input_dir):
    pdb_name_list = []  
    pdb_fpath_list = []
    # Devide to run.
    with open(args.split_json, 'r') as f:
        split_list = json.load(f)
    print(args.split_json)
    print(len(split_list))
    pdb_fname_list =[f for f in os.listdir(input_dir) if not f.endswith('relaxed.pdb')
            and f.endswith('.pdb')]
    for fname in pdb_fname_list:
        # print(fname)
        if "_NA_" in fname:  # Nanobody case
            new_fname = fname.replace("_NA_", "__")
        else:
            new_fname = fname
        new_fname_list = new_fname.split('_')
        pdb_name = '_'.join([new_fname_list[0], 
                            new_fname_list[1], 
                            new_fname_list[2], 
                            new_fname_list[3]])
        if pdb_name in split_list:
            fpath = os.path.join(input_dir, fname)
            # print('check')
            # print(pdb_name)
            # print(fpath)
            pdb_fpath_list.append(fpath)
            pdb_name_list.append(pdb_name)
            # print(fpath)
            # print(pdb_name)
    return pdb_fpath_list, pdb_name_list        


def filter_valid_atoms(input_pdb, output_pdb):
    if not os.path.exists(output_pdb):
        with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
            for line in infile:
                if line.startswith("ATOM"):
                    try:
                        # Extract the x, y, z coordinates
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                    except:
                        continue
                    # Check if the coordinates are not (0.000, 0.000, 0.000)
                    if (x, y, z) != (0.0, 0.0, 0.0):
                        outfile.write(line)
                        

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


def get_ref_files_list(dir_path, suffix="_reference.pdb"):
    ref_files_path_list = []
    pdb_files_name_list = []
    for filename in os.listdir(dir_path):
        if filename.endswith(suffix):
            fpath = os.path.join(dir_path, filename)
            ref_files_path_list.append(fpath)

            file_name = filename[: -len(suffix)]
            pdb_files_name_list.append(file_name)
            
    return ref_files_path_list, pdb_files_name_list
    
def get_ref_info():
    
    reference_data = {}
    ref_pdb_files_list, ref_pdb_name_list = get_ref_files_list(args.ref_dir)
    
    for ref_pdb_fpath, pdb_name in tqdm(
        zip(ref_pdb_files_list, ref_pdb_name_list), 
        total=len(ref_pdb_files_list)
        ):
        # if 'relaxed' not in ref_pdb:
        # pdb_name = ((ref_pdb.split('/')[-1]).split('.pdb'))[0]
        # print(pdb_name)
        code, heavy_chain_id, light_chain_id, antigen_chain = pdb_name.split('_')
        antibody_id_list = [heavy_chain_id, light_chain_id]
        ref_ab_ca, ref_ab_str_seq, ref_heavy_str, ref_light_str = boltz_make_coords(ref_pdb_fpath, antibody_id_list)
        # Using that our cdr definition.
        yaml_fpath = os.path.join(args.test_yaml_dir, f'{pdb_name}.yaml')
        heavy_cdr_mask, light_cdr_mask = read_yaml_file(yaml_fpath, heavy_chain_id, light_chain_id)
        try:
            _, heavy_cdr, light_cdr = boltz_cdr_numbering(ref_heavy_str, ref_light_str, heavy_cdr_mask, light_cdr_mask)
        except AssertionError:
            print(f'The pdb {pdb_name} have problem')
        data = {
            'heavy_cdr': heavy_cdr,
            'light_cdr': light_cdr,
            'coords': ref_ab_ca,
            'str_seq': ref_ab_str_seq
        }
        reference_data[f'{pdb_name}'] = data

    return reference_data
 
def boltz_relax():
    pass
 
def iggm_relax():
    pass

def dymean_relax():
    pass

def diffab_relax():
    pass

                        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, 
                        # default='/mnt/nas-new/home/yangnianzu/icml/baseline/diffab/our_single_cdr_h3_results_best_ckpt/codesign_single_boltz',
                        # default='/mnt/nas-new/home/yangnianzu/icml/baseline/rabd_c_result'
                        default='/mnt/nas-new/home/yangnianzu/icml/baseline/uniform_result'
                        )
    parser.add_argument('--output_dir', type=str, 
                        default='./tmp'
                        )
    parser.add_argument('--split_json', type=str, # split idx for multi nodes to relax.    
                        default='/mnt/nas-new/home/yangnianzu/icml/baseline/AbX/regular_list.json' # nano_list.json
                        # default='/mnt/nas-new/home/yangnianzu/icml/baseline/rabd_entry.json'
                        )
    parser.add_argument('--type', type=str,
                        default='boltz',  # like rfab
                        )
    parser.add_argument(
        "--cpus",
        type=int,
        default=mp.cpu_count(),
        help="The number of processes.",
    )
    parser.add_argument('--generate_area', type=str, choices=['cdrs', 'H3'], default='cdr3')
    parser.add_argument('-v', '--verbose', type=bool, default=True)
    
    parser.add_argument('--ref_dir', type=str,
                        default='/mnt/nas-new/home/yangnianzu/icml/data/all_structures/test_entry_pdb_files'
                        )
    parser.add_argument('--test_yaml_dir', type=str,  # for boltz test set or for rabd test set.
                        default='/mnt/nas-new/home/yangnianzu/icml/boltz/src/boltz/data/data_split/test_yml_dir/ab',
                        # default='/mnt/nas-new/home/yangnianzu/icml/boltz/src/boltz/data/data_split/rabd_yml_dir/ab', 
                        )
    parser.add_argument('--scheme', type=str, 
                        # default='chothia',
                        default='chothia'
                        )
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)
    logger = logging.getLogger(__name__)
    log_file = os.path.join(args.data_dir,'relax.log')
    handler_test = logging.FileHandler(log_file) # stdout to file
    handler_control = logging.StreamHandler()    # stdout to console

    selfdef_fmt = '%(asctime)s - %(funcName)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(selfdef_fmt)
    handler_test.setFormatter(formatter)
    handler_control.setFormatter(formatter)
    logger.setLevel('DEBUG')           #设置了这个才会把debug以上的输出到控制台
    logger.addHandler(handler_test)    #添加handler
    logger.addHandler(handler_control)

    
    # For diffab 
    # pdb_name_list = []
    # task_pdb_list = []
    # suffix_length = 21
    # with open(args.split_json, 'r') as f:
    #     split_list = json.load(f)
   
    # for gen_dir in os.listdir(args.data_dir):
    #     if 'log' not in gen_dir and '.csv' not in gen_dir:
    #         pdb_name = gen_dir[:-suffix_length]
    #         if pdb_name not in split_list:
    #             continue
    #         pdb_dir = os.path.join(args.data_dir, f'{gen_dir}/MultipleCDRs/')
    #         # pdb_dir = os.path.join(args.data_dir, f'{gen_dir}/H_CDR3/')
    #         count_relax_pdb = 0
    #         # Only consider about the pdb.
    #         list_pdb_dir = [pdb for pdb in os.listdir(pdb_dir) if 
    #                         not pdb.endswith('filter.pdb') and not pdb.endswith('relaxed.pdb')]
    #         # print(list_pdb_dir)q
    #         for pdb in sorted(list_pdb_dir):
    #             # print(pdb)
    #             if 'REF' not in pdb and count_relax_pdb <20:
    #                 count_relax_pdb += 1
    #                 input_pdb_fpath = os.path.join(pdb_dir, f'{pdb}')
    #                 filter_pdb_fpath = os.path.join(pdb_dir, f'{pdb[:4]}_filter.pdb')
    #                 filter_valid_atoms(input_pdb_fpath, filter_pdb_fpath)
    #                 # print(filter_pdb_fpath)
    #                 task_pdb_list.append(filter_pdb_fpath)
    #                 pdb_name_list.append(pdb_name)
    # original_pdb_list = task_pdb_list
    # print(len(split_list))
    # print(args.split_json)
    # print(len(original_pdb_list))
    
    # For dymean
    # task_pdb_list = []
    # pdb_name_list = get_pdb_files(args.data_dir)
    # original_pdb_list = []
    # for pdb_name in tqdm(pdb_name_list):
    #     # if pdb_name == '7vgs_D_C_A':
    #     #     continue
    #     pdb_dir = os.path.join(args.data_dir, pdb_name)
    #     original_pdb_fpath = os.path.join(pdb_dir, f'{pdb_name}_original.pdb')
    #     need_relax_pdb_fpath = os.path.join(pdb_dir, f'{pdb_name}.pdb')
    #     filter_out_pdb_fpath = os.path.join(pdb_dir, f'{pdb_name}_filter.pdb')
    #     filter_valid_atoms(need_relax_pdb_fpath, filter_out_pdb_fpath)
    #     task_pdb_list.append(filter_out_pdb_fpath)
    #     original_pdb_list.append(original_pdb_fpath)
    
    # For boltz
    task_pdb_list, pdb_name_list = get_boltz_pdb_names(args.data_dir)
    original_pdb_list = task_pdb_list
    print(len(task_pdb_list))
    print(len(pdb_name_list))


    # For rfab
    # task_pdb_list, pdb_name_list = get_rfab_pdb_names(args.data_dir)
    # original_pdb_list = task_pdb_list
    # print(len(task_pdb_list))   
    # print(len(pdb_name_list))

    # For iggm
    # task_pdb_list, pdb_name_list = get_iggm_pdb_names(args.data_dir)
    # original_pdb_list = task_pdb_list
    # print(len(task_pdb_list)) 
    # print(len(pdb_name_list))
    
    

    if args.scheme == 'imgt':
        ref_data = None
    elif args.scheme == 'chothia':
        ref_data = get_ref_info()
        
    main(args)


