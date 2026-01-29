#!/usr/bin/env python
"""
计算每个CDR区的AAR (Amino Acid Recovery)。
CDR区域由spec_mask中的'1'标识。
支持多sample模式：对每个case的多个sample分别计算AAR，然后取平均。
"""

import argparse
import json
import re
from pathlib import Path
from collections import defaultdict
import yaml


def find_cdr_regions(spec_mask):
    """从spec_mask中找出连续的'1'区域，返回[(start, end), ...]"""
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


def calc_aar(sequence, ground_truth, region):
    """计算指定区域的AAR"""
    start, end = region
    seq_region = sequence[start:end]
    gt_region = ground_truth[start:end]
    
    if len(seq_region) == 0:
        return 0.0
    
    matches = sum(1 for s, g in zip(seq_region, gt_region) if s == g)
    return matches / len(seq_region)


def calc_loop_aar(sequence, ground_truth, region, scheme='chothia'):
    """
    计算CDR-H3中间loop区域的AAR。
    IMGT scheme: 排除前4个和后2个残基 [4:-2]
    Chothia scheme: 排除前2个和后2个残基 [2:-2]
    """
    start, end = region
    seq_region = sequence[start:end]
    gt_region = ground_truth[start:end]
    
    # 根据scheme选择loop区域
    if scheme == 'imgt':
        seq_loop = seq_region[4:-2] if len(seq_region) > 6 else seq_region
        gt_loop = gt_region[4:-2] if len(gt_region) > 6 else gt_region
    else:  # chothia
        seq_loop = seq_region[2:-2] if len(seq_region) > 4 else seq_region
        gt_loop = gt_region[2:-2] if len(gt_region) > 4 else gt_region
    
    if len(seq_loop) == 0:
        return 0.0
    
    matches = sum(1 for s, g in zip(seq_loop, gt_loop) if s == g)
    return matches / len(seq_loop)


def parse_yaml_name(yaml_name):
    """解析yaml文件名，获取抗体链ID（去掉_sampleN后缀）"""
    name = yaml_name.replace('.yaml', '')
    # 去掉_sampleN后缀
    name = re.sub(r'_sample\d+$', '', name)
    parts = name.split('_')
    chain_parts = parts[1:]
    antigen_part = chain_parts[-1]
    antibody_parts = chain_parts[:-1]
    antibody_chains = [c for c in antibody_parts if c]
    return antibody_chains


def get_base_case_name(yaml_name):
    """从yaml文件名获取base case name（去掉_sampleN后缀）"""
    name = yaml_name.replace('.yaml', '')
    return re.sub(r'_sample\d+$', '', name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--yaml_dir', required=True, help='包含yaml文件的目录')
    parser.add_argument('--list_file', required=True, help='要评估的case列表(json)')
    args = parser.parse_args()
    
    yaml_dir = Path(args.yaml_dir)
    
    # 读取列表
    with open(args.list_file, 'r') as f:
        case_list = json.load(f)
    case_set = set(case_list)
    
    # CDR名称：重链H1,H2,H3；轻链L1,L2,L3
    heavy_cdr_names = ['H1', 'H2', 'H3']
    light_cdr_names = ['L1', 'L2', 'L3']
    
    # 收集每个case的所有sample结果: {case_name: [{cdr: aar, ...}, ...]}
    case_samples = defaultdict(list)
    
    # 遍历yaml_dir中的所有yaml文件
    for yaml_path in yaml_dir.glob('*.yaml'):
        base_case_name = get_base_case_name(yaml_path.name)
        
        # 只处理在list中的case
        if base_case_name not in case_set:
            continue
        
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        
        # 获取抗体链ID
        antibody_chains = parse_yaml_name(yaml_path.name)
        
        # 按顺序找到抗体链的蛋白质数据
        chain_data = {}
        for protein in data['sequences']:
            chain_id = protein['protein']['id']
            if chain_id in antibody_chains:
                chain_data[chain_id] = protein['protein']
        
        result = {}
        
        # 处理重链（第一条抗体链）
        if len(antibody_chains) >= 1:
            heavy_id = antibody_chains[0]
            if heavy_id in chain_data:
                heavy = chain_data[heavy_id]
                cdr_regions = find_cdr_regions(heavy['spec_mask'])
                for i, region in enumerate(cdr_regions):
                    cdr_name = heavy_cdr_names[i] if i < len(heavy_cdr_names) else f'H{i+1}'
                    aar = calc_aar(heavy['sequence'], heavy['ground_truth'], region)
                    result[cdr_name] = aar
                    # 对H3计算Loop-AAR
                    if cdr_name == 'H3':
                        loop_aar = calc_loop_aar(heavy['sequence'], heavy['ground_truth'], region)
                        result['H3_Loop'] = loop_aar
        
        # 处理轻链（第二条抗体链，如果有）
        if len(antibody_chains) >= 2:
            light_id = antibody_chains[1]
            if light_id in chain_data:
                light = chain_data[light_id]
                cdr_regions = find_cdr_regions(light['spec_mask'])
                for i, region in enumerate(cdr_regions):
                    cdr_name = light_cdr_names[i] if i < len(light_cdr_names) else f'L{i+1}'
                    aar = calc_aar(light['sequence'], light['ground_truth'], region)
                    result[cdr_name] = aar
        
        if result:
            case_samples[base_case_name].append(result)
    
    # 对每个case计算sample平均值
    all_results = []
    for case_name in case_list:
        if case_name not in case_samples:
            continue
        
        samples = case_samples[case_name]
        if not samples:
            continue
        
        # 计算每个CDR的平均值
        avg_result = {'name': case_name}
        all_cdrs = set()
        for s in samples:
            all_cdrs.update(s.keys())
        
        for cdr in all_cdrs:
            values = [s[cdr] for s in samples if cdr in s]
            if values:
                avg_result[cdr] = sum(values) / len(values)
        
        all_results.append(avg_result)
    
    # 打印结果
    if not all_results:
        print("No results found.")
        return
    
    # 获取所有CDR列
    all_cdrs = set()
    for r in all_results:
        all_cdrs.update(k for k in r.keys() if k != 'name')
    
    def sort_key(x):
        # H1, H2, H3, H3_Loop, L1, L2, L3
        if '_' in x:  # H3_Loop
            base, suffix = x.split('_', 1)
            return (base[0], int(base[1:]), suffix)
        return (x[0], int(x[1:]), '')
    
    cdr_cols = sorted(all_cdrs, key=sort_key)
    
    # 打印表头
    header = ['name'] + cdr_cols
    print('\t'.join(header))
    
    # 打印平均值
    print('-' * 80)
    avg_row = ['Average']
    for cdr in cdr_cols:
        values = [r[cdr] for r in all_results if cdr in r]
        if values:
            avg_row.append(f'{sum(values)/len(values):.4f}')
        else:
            avg_row.append('-')
    print('\t'.join(avg_row))
    
    # 打印统计信息
    print(f'\nTotal cases: {len(all_results)}')
    print(f'Samples per case: {len(case_samples[all_results[0]["name"]]) if all_results else 0}')


if __name__ == '__main__':
    main()
