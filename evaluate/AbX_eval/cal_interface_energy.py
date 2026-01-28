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
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB import PDBParser
from Bio import PDB
import random
from tqdm import tqdm
from copy import copy

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

#Core Includes
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.simple_metrics import metrics
from pyrosetta.rosetta.core.select import residue_selector as selections
from pyrosetta.rosetta.core import select
from pyrosetta.rosetta.core.select.movemap import *
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, AndResidueSelector, OrResidueSelector
from pyrosetta import create_score_function

#Protocol Includes
from pyrosetta.rosetta.protocols import minimization_packing as pack_min
from pyrosetta.rosetta.protocols import relax as rel
from pyrosetta.rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from pyrosetta.rosetta.protocols.antibody import *
from pyrosetta.rosetta.protocols.loops import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, move_map_action
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

init('-use_input_sc -input_ab_scheme AHo_Scheme -ignore_unrecognized_res \
    -ignore_zero_occupancy false -load_PDB_components true -relax:default_repeats 2 -no_fconfig')

# have_name = list(pd.read_csv('./na.idx', names=['name'], header=None)['name'])

def parse_list(name_idx, pdb_dir):
    names = list(pd.read_csv(name_idx, names=['name'], header=None)['name'])
    names = names[::-1]
    random.shuffle(names)
    
    for name in names:
        if name in have_name:
            continue
        code, heavy_chain, light_chain, antigen_chain = name.split('_')
        def _parse_chain_id(heavy_chain_id, light_chain_id):
            if heavy_chain_id.islower() and heavy_chain_id.upper() == light_chain_id:
                heavy_chain_id = heavy_chain_id.upper()
            elif light_chain_id.islower() and light_chain_id.upper() == heavy_chain_id:
                light_chain_id = light_chain_id.upper()
            return heavy_chain_id, light_chain_id
        heavy_chain_, light_chain_ = _parse_chain_id(heavy_chain, light_chain)
        pdb_file = os.path.join(pdb_dir, f'{code}_{heavy_chain_}_{light_chain_}_{antigen_chain}' ,f'{code}_{heavy_chain_}{light_chain_}{antigen_chain}_ab_ag.pdb')
        yield pdb_file 
        
        
def parse_boltz_raw_pdb_list(pdb_dir, suffix='_reference.pdb'):
    sub_path_list = [os.path.join(pdb_dir, item) for item in os.listdir(pdb_dir)]
    for sub_fpath in sub_path_list:
        pdb_name = sub_fpath.split('/')[-1]
        pdb_core_name = pdb_name.rsplit(suffix, 1)[0]
        yield sub_fpath, pdb_core_name, pdb_core_name
        
def parse_our_design_pdb_list(pdb_dir, suffix='_relaxed.pdb'):
    # sub_path_list = [os.path.join(pdb_dir, item) for item in os.listdir(pdb_dir)
    #                  if not item.endswith('.log') and not item.endswith('.csv')]
    # for sub_fpath in sub_path_list:
    #     # pdb_fpath_list = []
    #     for item in os.listdir(sub_fpath):
    #         if not item.endswith('_relaxed_relaxed.pdb'):
    #             if item.endswith(suffix):
    #                 pdb_fpath = os.path.join(sub_fpath, item)
    #                 pdb_core_name = sub_fpath.split('/')[-1]
                    
    #                 order_name = ['A', 'B', 'C', 'D', 'E', 'F']
    #                 pdb_core_name_list = pdb_core_name.split('_')
    #                 code = pdb_core_name_list[0]
    #                 sum_chain = len(''.join(pdb_core_name_list[1:]))
    #                 select_name = order_name[:sum_chain]
    #                 anti_count = len(pdb_core_name_list[1] + pdb_core_name_list[2])
    #                 if anti_count == 1:
    #                     heavy_name = select_name[0]
    #                     light_name = ''
    #                     ag_name = ''.join(select_name[1:])
    #                 elif anti_count == 2:
    #                     heavy_name = select_name[0]
    #                     light_name = select_name[1]
    #                     ag_name = ''.join(select_name[2:])
    #                 new_pdb_core_name = f'{code}_{heavy_name}_{light_name}_{ag_name}'  
    #                 # print(pdb_core_name)
    #                 # print(new_pdb_core_name)
    #                 yield pdb_fpath, new_pdb_core_name, pdb_core_name
    sub_path_list = [item for item in os.listdir(pdb_dir)
                     if not item.endswith('.log') and not item.endswith('.csv')]
    for item_pdb in sub_path_list:
        # print(item_pdb)
        # pdb_fpath_list = []
        item_pdb = os.path.join(pdb_dir, item_pdb)
        
        if not item_pdb.endswith('_relaxed_relaxed.pdb'):
            if item_pdb.endswith(suffix):
                pdb_fpath = os.path.join(pdb_dir, item_pdb)
                # print(pdb_fpath)
                item_pdb_list = item_pdb.split('_')
                pdb_core_name = f'{item_pdb_list[0]}_{item_pdb_list[1]}_{item_pdb_list[2]}_{item_pdb_list[3]}'
                # print(pdb_core_name)
                # pdb_core_name = sub_fpath.split('/')[-1]
                order_name = ['A', 'B', 'C', 'D', 'E', 'F']
                pdb_core_name_list = pdb_core_name.split('_')
                code = pdb_core_name_list[0]
                sum_chain = len(''.join(pdb_core_name_list[1:]))
                select_name = order_name[:sum_chain]
                anti_count = len(pdb_core_name_list[1] + pdb_core_name_list[2])
                if anti_count == 1:
                    heavy_name = select_name[0]
                    light_name = ''
                    ag_name = ''.join(select_name[1:])
                elif anti_count == 2:
                    heavy_name = select_name[0]
                    light_name = select_name[1]
                    ag_name = ''.join(select_name[2:])
                new_pdb_core_name = f'{code}_{heavy_name}_{light_name}_{ag_name}'  
                # print(new_pdb_core_name)
                # print(pdb_core_name)
                # print(new_pdb_core_name)
                yield pdb_fpath, new_pdb_core_name, pdb_core_name


def parse_diffab_design_pdb_list(pdb_dir, suffix_length=21):
        gen_dir_list = [
            gen_dir for gen_dir in os.listdir(pdb_dir)
            if (
                os.path.isdir(os.path.join(pdb_dir, gen_dir))  # Only directories
                and 'log' not in gen_dir
                and '.csv' not in gen_dir
            )
                        
        ]
        for gen_dir in gen_dir_list:
            pdb_name = gen_dir[:-suffix_length]
            sub_pdb_dir = os.path.join(pdb_dir, f'{gen_dir}/MultipleCDRs/')
            count_relax_pdb = 0
            # Only consider about the pdb.
            list_pdb_dir = [pdb for pdb in os.listdir(sub_pdb_dir) 
                            if pdb.endswith('_filter_relaxed.pdb')]
            # print(list_pdb_dir)q
            for pdb in sorted(list_pdb_dir):  
                if count_relax_pdb < 20:
                    count_relax_pdb += 1
                    input_pdb_fpath = os.path.join(sub_pdb_dir, f'{pdb}')

                    # print(input_pdb_fpath)
                    # print(pdb_name)
                    yield input_pdb_fpath, pdb_name, pdb_name
                    
def parse_dymean_design_pdb_list(pdb_dir,):
    pdb_name_list = [f for f in os.listdir(pdb_dir) if 
                     not f.endswith('json') and not f.endswith('log')
                     and not f.endswith('csv')]
    for pdb_name in tqdm(pdb_name_list):
        # if pdb_name == '7vgs_D_C_A':
        #     continue
        sub_pdb_dir = os.path.join(pdb_dir, pdb_name)
        relaxed_pdb_fpath = os.path.join(sub_pdb_dir, f'{pdb_name}_openmm_relaxed.pdb')

        # print(relaxed_pdb_fpath)
        # print(pdb_name)
        yield relaxed_pdb_fpath, pdb_name, pdb_name


def parse_iggm_design_pdb_list(pdb_dir,):
    list_pdb_dir = [pdb for pdb in os.listdir(pdb_dir) if pdb.endswith('relaxed.pdb')]
    for pdb in sorted(list_pdb_dir):
        input_pdb_fpath = os.path.join(pdb_dir, f'{pdb}')
        # print(input_pdb_fpath)
        # print(pdb_name)
        pdb_name = pdb.split('_All')[0]
        if "_NA_" in pdb_name:  # Nanobody case
            pdb_name = pdb_name.replace("_NA_", "__")
        old_pdb_name = copy(pdb_name)

        pdb_name_list = pdb_name.split('_')
        antigen_ids = pdb_name_list[3]
        if len(antigen_ids) > 1:
            antigen_id = antigen_ids[0]
        else:
            antigen_id = antigen_ids
        pdb_name = f'{pdb_name_list[0]}_{pdb_name_list[1]}_{pdb_name_list[2]}_{antigen_id}'
        yield input_pdb_fpath, pdb_name, old_pdb_name

def parse_rfab_design_pdb_list(pdb_dir,):
    list_pdb_dir = [pdb for pdb in os.listdir(pdb_dir) if pdb.endswith('relaxed.pdb')]
    for pdb in sorted(list_pdb_dir):
        input_pdb_fpath = os.path.join(pdb_dir, f'{pdb}')
        # print(input_pdb_fpath)
        # print(pdb_name)
        old_pdb_name = '_'.join(pdb.split('_')[:4])

        pdb_list = pdb.split('_')[:4]
        
        heavy_id = 'H'
        if pdb_list[1] != '':
            light_id = 'L'
        else:
            light_id = ''
        
        antigen_id = 'T'
        pdb_name = f'{pdb_list[0]}_{heavy_id}_{light_id}_{antigen_id}'
        yield input_pdb_fpath, pdb_name, old_pdb_name


def pyrosetta_interface_energy(pdb_path, interface):
    pose = pyrosetta.pose_from_pdb(pdb_path)
    mover = InterfaceAnalyzerMover()
    mover.set_interface(interface)
    mover.set_scorefunction(pyrosetta.create_score_function('ref2015'))
    mover.apply(pose)
    print(pose.scores)
    return pose.scores['dG_separated']


def InterfaceEnergy(origin_pdb_file, pdb_name):
    logger.info(f"Calculate Rosetta Interface Energy for {origin_pdb_file}")
    # code, heavy_chain_id, light_chain_id, antigen_chain_ids = (origin_pdb_file.split('/')[-2]).split('_')
    code, heavy_chain_id, light_chain_id, antigen_chain = pdb_name.split('_')
    parser = PDBParser()
    struc = parser.get_structure(code, origin_pdb_file)[0]
    antibody_chains = [heavy_chain_id, light_chain_id]
    # antigen_chains = set()
    # for chain in struc:
    #     if chain.id in antibody_chains:
    #         continue
    #     if chain.id in antigen_chain_ids:
    #         antigen_chains.add(chain.id)
    antigen_chains = ''.join(antigen_chain)
    antibody_chains = ''.join(antibody_chains)
    interface = f"{antibody_chains}_{antigen_chains}"

    print(f'pdb_name: {pdb_name}')
    print(f'interface: {interface}')
    dG_ref = pyrosetta_interface_energy(origin_pdb_file, interface)
    # dG_ref = 0
    return dG_ref

def process(pdb_file, pdb_name, old_pdb_name):
    # try:
    dG_ref = InterfaceEnergy(pdb_file, pdb_name)
    logger.info(f"{pdb_file.split('/')[-2]}@dG_wild: {dG_ref}")
    # dG_ref = 0
    return old_pdb_name, dG_ref
    # except:
    #     return old_pdb_name, 0

def main(args):
   
    # TODO: parser_file_function if different.
    # 1. boltz: parse_boltz_raw_pdb_list
    # 2. our_design: parse_our_design_pdb_list
    # 3. diffab_design: parse_diffab_design_pdb_list
    # 4. dymean_design:  parse_dymean_design_pdb_list
    # 5. iggm_design: parse_iggm_design_pdb_list
    # 6. rfab_design: parse_rfab_design_pdb_list
    parse_function = parse_our_design_pdb_list
    # for pdb_file, pdb_name, old_pdb_name in parse_function(args.pdb_dir):
        # print(pdb_file)
        # print(pdb_name)

    with mp.Pool(args.cpus) as p:
        results = p.starmap(process, [(pdb_file, pdb_name, old_pdb_name) for pdb_file, pdb_name, old_pdb_name in 
                            parse_function(args.pdb_dir)])

    # 保存到 CSV 文件
    df = pd.DataFrame(results, columns=["PDB_Name", "dG"])
    csv_path = os.path.join(args.output_dir, args.csv_name)  
    df.to_csv(csv_path, index=False)
    print(f"Results saved to {csv_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--pdb_dir', type=str, 
                        required=False,
                        default='/mnt/nas-new/home/yangnianzu/icml/baseline/uniform_result'
                        )
    parser.add_argument('-n', '--name_idx', type=str, 
                        required=False,
                        default=None,
                        )
    parser.add_argument('-o', '--output_dir', type=str, 
                        required=False,
                        default='./imp_results/',
                        )
    parser.add_argument('--csv_name', type=str,
                        default='mfdesign_uniform_imp_energy_validation.csv',
                        )
    parser.add_argument('-c', '--cpus', type=int, default=mp.cpu_count())
    parser.add_argument('-v', '--verbose', type=bool, default=False)
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)
    logger = logging.getLogger(__name__)
    log_file = os.path.join(args.output_dir,'diffab_new_native_energy.log')
    handler_test = logging.FileHandler(log_file) # stdout to file
    handler_control = logging.StreamHandler()    # stdout to console

    selfdef_fmt = '%(asctime)s - %(funcName)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(selfdef_fmt)
    handler_test.setFormatter(formatter)
    handler_control.setFormatter(formatter)
    logger.setLevel('DEBUG')           #设置了这个才会把debug以上的输出到控制台
    logger.addHandler(handler_test)    #添加handler
    logger.addHandler(handler_control)
    
    main(args)

