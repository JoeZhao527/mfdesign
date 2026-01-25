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
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import PDB
from tqdm import tqdm
import traceback

from multiprocessing import Pool
from joblib import Parallel, delayed
import multiprocessing
import joblib

import sys
sys.path.insert(0, '/mnt/nas-new/home/yangnianzu/icml/baseline/AbX')

from abx.common import residue_constants
from abx.data.mmcif_parsing import parse as mmcif_parse
from abx.preprocess.numbering import renumber_ab_seq, get_ab_regions

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
import pyrosetta.distributed

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
three_to_one = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'
}

init('-use_input_sc -input_ab_scheme AHo_Scheme -ignore_unrecognized_res \
    -ignore_zero_occupancy false -load_PDB_components true -relax:default_repeats 2 -no_fconfig')


def find_all_indices(input_list, value_to_find):
    return [index for index, value in enumerate(input_list) if value == value_to_find]


def get_seqres_from_pdb(structure, chain_id):
    """Extract SEQRES sequence for a specific chain from a PDB file."""
    # parser = PDBParser(QUIET=True)
    # structure = parser.get_structure('structure', pdb_filename)
    model = structure[0]  # Assuming only one model in the PDB
    chain = model[chain_id]
    seq = ""
    for res in chain.get_unpacked_list():
        resname = res.get_resname().strip()
        if resname in three_to_one:
            seq += three_to_one[resname]

    return seq


def cdr_domain(pred_antibody, code, orginal_antibody, heavy_chain_id, light_chain_id, antigen_chain_ids, 
               heavy_cdr_def, light_cdr_def, scheme, logger,
               ):
    def _make_domain(feature, chain_id, cdr_def, scheme='chothia'):
        allow = ['H'] if chain_id == 'H' else ['K', 'L']
        # print(scheme)
        if cdr_def is None:
            anarci_res = renumber_ab_seq(feature['str_seq'], allow=allow, scheme=scheme)
            print(f"str_seq: {feature['str_seq']}, anarci_res: {anarci_res}")
            domain_numbering, domain_start, domain_end = map(anarci_res.get, ['domain_numbering', 'start', 'end'])
            assert domain_numbering is not None, print(chain_id)
            
            cdr_def = get_ab_regions(domain_numbering, chain_id=chain_id)
            
        CDR_indices = dict()
        if allow == ['H']:
            CDR_H1_indice = find_all_indices(cdr_def, 1)
            CDR_H2_indice = find_all_indices(cdr_def, 3)
            CDR_H3_indice = find_all_indices(cdr_def, 5)
            CDR_indices.update(
                CDR_H1 = [min(CDR_H1_indice)+1, max(CDR_H1_indice)+1],
                CDR_H2 = [min(CDR_H2_indice)+1, max(CDR_H2_indice)+1],
                CDR_H3 = [min(CDR_H3_indice)+1, max(CDR_H3_indice)+1],
            )
        if allow == ['K', 'L']:
            CDR_L1_indice = find_all_indices(cdr_def, 8)
            CDR_L2_indice = find_all_indices(cdr_def, 10)
            CDR_L3_indice = find_all_indices(cdr_def, 12)
            CDR_indices.update(
                CDR_L1 = [min(CDR_L1_indice)+1, max(CDR_L1_indice)+1],
                CDR_L2 = [min(CDR_L2_indice)+1, max(CDR_L2_indice)+1],
                CDR_L3 = [min(CDR_L3_indice)+1, max(CDR_L3_indice)+1],
            )
            
        return CDR_indices
    
    try:
        parser = PDBParser()
        # struc = parser.get_structure(code, origin_antibody)
        pred_struc = parser.get_structure(code, pred_antibody)
        orginal_struc = parser.get_structure(code, orginal_antibody)
    except:
        raise ValueError('PDB_parse: %s', pred_antibody)
    

    antigen_chain_ids = list(antigen_chain_ids)
    try:
        if heavy_chain_id != '':
            heavy_str_seq = get_seqres_from_pdb(pred_struc, heavy_chain_id)
            heavy_data = dict(
                str_seq = heavy_str_seq,
            )
            gt_heavy_str_seq = get_seqres_from_pdb(orginal_struc, heavy_chain_id)
            gt_heavy_data = dict(
                str_seq = gt_heavy_str_seq,
            )
        else:
            heavy_data = None
            gt_heavy_data = None
            
        if light_chain_id != '':
            light_str_seq = get_seqres_from_pdb(pred_struc, light_chain_id)
            light_data = dict(
                str_seq = light_str_seq,
            )
            gt_light_str_seq = get_seqres_from_pdb(orginal_struc, light_chain_id)
            gt_light_data = dict( 
                str_seq = gt_light_str_seq,
            )
        else:
            light_data = None
            gt_light_data = None
    except:
        print(f'heavy{heavy_chain_id}, light{light_chain_id}')
        
    try:
        if heavy_data and heavy_cdr_def is not None:
            assert len(heavy_str_seq) == len(heavy_cdr_def)
        if light_data and light_cdr_def is not None:
            assert len(light_str_seq) == len(light_cdr_def)
    except:
        if heavy_data:
            heavy_cdr_def = None
            
        if light_data:
            light_cdr_def = None
        
        # if len(heavy_str_seq) + 1 == len(heavy_cdr_def):
        #     heavy_cdr_def = heavy_cdr_def[1:]   # This consider for the diffab situation.
        logger.info("May exist two antibody that the length is not equal to cdr def")
        logger.info(f"{code}")
        # pdb
        
    CDR_indice = dict()
    if heavy_data:
        heavy_cdr_indice = _make_domain(gt_heavy_data, 'H', heavy_cdr_def, scheme=scheme)
        CDR_indice.update(
            heavy_cdr_indice,
        )
    if light_data:
        light_cdr_indice = _make_domain(gt_light_data, 'L', light_cdr_def, scheme=scheme)
        CDR_indice.update(
            light_cdr_indice,
        )

    return CDR_indice

def apply_threshold_to_cdr_dict(cdr_dict, threshold, chianed_id):
    """
    Adjust the indices in the CDR dictionary by applying a threshold.

    Args:
        cdr_dict (dict): The original CDR dictionary with residue ranges.
        threshold (int): The number of residues to expand the range by.
        pose (Pose): The PyRosetta pose object to ensure valid indices.

    Returns:
        dict: The updated CDR dictionary with adjusted indices.
    """
    adjusted_cdr_dict = {}

    print(cdr_dict)
    for cdr_name, indices in cdr_dict.items():
        start_idx = indices[0]
        end_idx = indices[1]

        if 'L' in cdr_name:
            start_idx = start_idx + threshold
            end_idx = end_idx + threshold
            adjusted_cdr_dict[cdr_name] = [start_idx, end_idx]
        else:   
            adjusted_cdr_dict[cdr_name] = [start_idx, end_idx]

    return adjusted_cdr_dict

def Rosetta_relax(pdb_file, pdb_name, original_pdb_file, args, ref_data, type=None, logger=None):
    _, heavy_chain_id, light_chain_id, antigen_ids = pdb_name.split('_')
    print(f'Heavy_{heavy_chain_id}')
    print(f'Light_{light_chain_id}')
    print('---------------')
    output_dir, output = os.path.dirname(pdb_file), os.path.basename(pdb_file).split('.')[0]
    # print(pdb_file)
    # print(output_dir)
    # print(output)
    output_file = os.path.join(output_dir, f'{output}_relaxed.pdb')
    if os.path.exists(output_file):
        print('exists')
        return None

    if args.scheme == 'chothia':
        h_cdr_def = ref_data[pdb_name]['heavy_cdr']
        print(h_cdr_def)
        l_cdr_def = ref_data[pdb_name]['light_cdr']
        print(l_cdr_def)
    elif args.scheme == 'imgt':
        h_cdr_def = None
        l_cdr_def = None
    if type == 'boltz':
        # Here need to reorder the chain id.
        order_chain = ['A', 'B']
        if heavy_chain_id != '':
            heavy_chain_id = order_chain[0]
        if light_chain_id != '':
            light_chain_id = order_chain[1]
    elif type == 'rfab':
        order_chain = ['H', 'L', 'T']
        if heavy_chain_id != '':
            heavy_chain_id = order_chain[0]
        if light_chain_id != '':
            light_chain_id = order_chain[1]
        if antigen_ids != '':
            antigen_ids = order_chain[2]
    
    try:
        cdr_dict = cdr_domain(pdb_file, pdb_name, original_pdb_file, heavy_chain_id, light_chain_id, list(antigen_ids),
                              h_cdr_def, l_cdr_def, args.scheme, logger,
                              )
    except:
        return None
    # cdr_dict = cdr_domain(pdb_file, pdb_name, original_pdb_file, heavy_chain_id, light_chain_id, list(antigen_ids),
    #                       h_cdr_def, l_cdr_def, args.scheme, logger,
    #                       )
    # print(h_cdr_def)
    if type == 'rfab':
        cdr_dict = apply_threshold_to_cdr_dict(cdr_dict, len(h_cdr_def), light_chain_id)

    logging.info(f'Rosetta processing {pdb_file} for Relax')
    pose = pose_from_pdb(pdb_file)
    scorefxn = create_score_function('ref2015')
    pack_mover = PackRotamersMover()
    
    # pdb_name = pdb_file.split('/')[-1]
    # if '@' in pdb_name:
    #     code, heavy_chain_id, light_chain_id,antigen_chain_ids = (pdb_name.split('@')[0]).split('_')
    #     output = '.'.join(pdb_name.split('.')[:-1])
    # else:
    #     code, heavy_chain_id, light_chain_id,antigen_chain_ids = (pdb_name.split('.')[0]).split('_')
    #     output = pdb_name.split('.')[0]
    

    
    original_pose = pose.clone()
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking()) 
    tf.push_back(operation.PreventRepacking()) 
    flexible_dict = dict()
    count = 0

    
    if args.generate_area == 'H3':
        cdr_dict = dict(
            CDR_H3 = cdr_dict['CDR_H3'],
        )

    for cdr_name, indice in cdr_dict.items():
        if count < 3:
            flexible_dict.update(
                {cdr_name: [(heavy_chain_id, indice[0]), (heavy_chain_id, indice[1])]}
            )
        else:
            flexible_dict.update(
                {cdr_name: [(light_chain_id, indice[0]), (light_chain_id, indice[1])]}
            ) 
        count += 1       
    gen_selector = selections.ResidueIndexSelector('1')
    # print(f"generate_area: {flexible_dict}")
    print(flexible_dict)
    for cdr_name, indice in flexible_dict.items():
        flexible_residue_first = indice[0]
        flexible_residue_last = indice[1]
        gen_selector1 = selections.ResidueIndexSelector()
        gen_selector1.set_index_range(
            pose.pdb_info().pdb2pose(*flexible_residue_first), 
            pose.pdb_info().pdb2pose(*flexible_residue_last), 
        )
        gen_selector = OrResidueSelector(gen_selector, gen_selector1)
    nbr_selector = selections.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(gen_selector)
    nbr_selector.set_include_focus_in_subset(True)
    subset_selector = nbr_selector
    prevent_repacking_rlt = operation.PreventRepackingRLT()
    prevent_subset_repacking = operation.OperateOnResidueSubset(
        prevent_repacking_rlt, 
        subset_selector,
        flip_subset=True,
    )
    tf.push_back(prevent_subset_repacking)
    print(pose)
    packer_task = tf.create_task_and_apply_taskoperations(pose)

    movemap = MoveMapFactory()
    movemap.add_bb_action(move_map_action.mm_enable, gen_selector)  # 允许第i个残基的骨架移动
    movemap.add_chi_action(move_map_action.mm_enable, subset_selector) # 允许第i个残基的侧链移动
    mm = movemap.create_movemap_from_pose(pose)
    fastrelax = FastRelax()
    fastrelax.set_scorefxn(scorefxn)
    fastrelax.set_movemap(mm) #使用默认的Movemap()
    fastrelax.set_task_factory(tf)
    fastrelax.apply(pose)
    
    pose.dump_pdb(f'{output_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--generate_area', type=str, default='all')
    parser.add_argument('--data_dir', type=str, 
                        default='/mnt/nas-new/home/yangnianzu/icml/baseline/dyMEAN/results/',
                        )
    parser.add_argument('--pdb', type=str, 
                        default='/mnt/nas-new/home/yangnianzu/icml/baseline/dyMEAN/results/8zc1_E_C_B_filter.pdb'
                        )
    parser.add_argument('--pdb-name', type=str, 
                        default='8zc1_E_C_B.pdb'
                        )
    args = parser.parse_args()
    
    Rosetta_relax(args.pdb, args.pdb_name, args=args)
    
    # relax_pdb_fpath = '/mnt/nas-new/home/yangnianzu/icml/baseline/diffab/diffab/tools/runner/results/codesign_multicdrs_boltz/7vgs_D_C_A_2025_01_21__21_52_35/MultipleCDRs/0000.pdb'
    # pdb_name_list = get_pdb_files(args.data_dir)
    # task_pdb_list = []
    # for pdb_name in tqdm(pdb_name_list):
    #     need_relax_pdb_fpath = os.path.join(args.data_dir, pdb_name)
    #     filter_out_pdb_fpath = os.path.join(args.data_dir, f'{pdb_name[:-4]}_filter.pdb')
    #     filter_valid_atoms(need_relax_pdb_fpath, filter_out_pdb_fpath)
    #     task_pdb_list.append(filter_out_pdb_fpath)
    

