#!/usr/bin/env python
"""
从fasta文件中选取score最高的序列，更新对应的yaml文件。
- 替换抗体链的sequence
- 将抗原链的spec_mask全改成0
"""

import argparse
import re
from pathlib import Path
import yaml


def parse_fasta(fasta_path):
    """解析fasta文件，返回score最高的sample序列"""
    samples = []
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        header = lines[i].strip()
        if header.startswith('>'):
            seq = lines[i + 1].strip() if i + 1 < len(lines) else ''
            # 只考虑sample行，不要原始序列
            if 'sample=' in header:
                score_match = re.search(r'score=([0-9.]+)', header)
                if score_match:
                    score = float(score_match.group(1))
                    samples.append({'header': header, 'sequence': seq, 'score': score})
            i += 2
        else:
            i += 1
    
    if not samples:
        return None
    
    # 选score最高的
    best = max(samples, key=lambda x: x['score'])
    return best


def parse_yaml_name(yaml_name):
    """
    解析yaml文件名，获取抗体链和抗原链ID。
    例如: 8r9y_B_C_A.yaml -> antibody=['B','C'], antigen=['A']
         9bu8_C_D_AB.yaml -> antibody=['C','D'], antigen=['A','B']
         8zes_E__D.yaml -> antibody=['E'], antigen=['D'] (双下划线表示只有一条抗体链)
    """
    name = yaml_name.replace('.yaml', '')
    parts = name.split('_')
    # 第一部分是PDB ID，剩余是链ID
    chain_parts = parts[1:]
    
    # 最后一个部分是抗原链（可能是多个字母如AB）
    antigen_part = chain_parts[-1]
    antibody_parts = chain_parts[:-1]
    
    # 过滤空字符串（处理双下划线情况如 8zes_E__D）
    antibody_chains = [c for c in antibody_parts if c]
    antigen_chains = list(antigen_part)
    
    return antibody_chains, antigen_chains


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_dir', required=True, help='包含fasta文件的目录')
    parser.add_argument('--yaml_dir', required=True, help='包含yaml文件的目录')
    parser.add_argument('--output_dir', required=True, help='输出目录')
    args = parser.parse_args()
    
    fasta_dir = Path(args.fasta_dir)
    yaml_dir = Path(args.yaml_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 遍历所有fasta文件
    for fasta_file in fasta_dir.glob('*.fa'):
        # 从文件名解析base name: 8r9y_H_L_A|xxx.fa -> 8r9y_H_L_A
        base_name = fasta_file.stem.split('|')[0]
        
        # 找对应的yaml文件
        yaml_path = yaml_dir / f'{base_name}.yaml'
        if not yaml_path.exists():
            print(f'Warning: {yaml_path} not found, skipping {fasta_file.name}')
            continue
        
        # 解析fasta获取最佳序列
        best_sample = parse_fasta(fasta_file)
        if not best_sample:
            print(f'Warning: No samples found in {fasta_file.name}, skipping')
            continue
        
        # 读取yaml
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)
        
        # 从yaml文件名获取链信息
        antibody_chains, antigen_chains = parse_yaml_name(yaml_path.name)
        
        # 分割序列（用/分隔）
        seqs = best_sample['sequence'].split('/')
        assert len(seqs) == len(antibody_chains), \
            f'{fasta_file.name}: 序列数量({len(seqs)})与抗体链数量({len(antibody_chains)})不匹配'
        
        # 更新yaml
        skip_this = False
        for protein in data['sequences']:
            chain_id = protein['protein']['id']
            
            if chain_id in antibody_chains:
                # 替换抗体链序列
                idx = antibody_chains.index(chain_id)
                new_seq = seqs[idx]
                old_seq = protein['protein']['sequence']
                if len(new_seq) != len(old_seq):
                    print(f'Warning: {fasta_file.name}: 链{chain_id}序列长度不匹配({len(new_seq)} vs {len(old_seq)}), skipping')
                    skip_this = True
                    break
                protein['protein']['sequence'] = new_seq
            
            elif chain_id in antigen_chains:
                # 抗原链的spec_mask全改成0
                old_mask = protein['protein']['spec_mask']
                protein['protein']['spec_mask'] = '0' * len(old_mask)
        
        if skip_this:
            continue
        
        # 写入输出
        output_path = output_dir / yaml_path.name
        with open(output_path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        
        # print(f'Processed: {fasta_file.name} (score={best_sample["score"]:.4f}) -> {output_path.name}')


if __name__ == '__main__':
    main()

