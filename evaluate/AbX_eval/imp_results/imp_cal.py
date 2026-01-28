import pandas as pd
import json
import math as ma
import numpy as np

import sys

ground_truth_result_csv = 'evaluate/AbX_eval/imp_results/test_ground_binding.csv'
sample_imp_result_csv = sys.argv[1]

ground_truth_result = pd.read_csv(ground_truth_result_csv)
sample_imp_result = pd.read_csv(sample_imp_result_csv)

def count_success(group, ref_value):
    # print(group['dG'])
    # print(ref_value)
    # print(group["dG"] < ref_value)
    return (group["dG"] < ref_value).sum()

def get_ref_dG(df, pdb_name):
    result = df[df["PDB_Name"] == pdb_name]["dG"]
    # print('ref_result', result)
    return result.values[0] if not result.empty else None 

anti_json_fpath = 'evaluate/AbX_eval/regular_list.json'
with open(anti_json_fpath, 'r') as f:
    anti_json = json.load(f)

nano_json_fpath = 'evaluate/AbX_eval/nano_list.json'
with open(nano_json_fpath, 'r') as f:
    nano_json = json.load(f)

def find_anti_idx(df: object, idx_json: list) -> bool:
    idx_list = []
    for f in df['PDB_Name']:
        base_name = f
        if base_name in idx_json:
            idx_list.append(True)
        else:
            idx_list.append(False)
    return idx_list

def split_get_imp_result(df: object, anti_json) -> tuple:
    anti_idx = find_anti_idx(df, anti_json)
    nano_idx = ~np.array(anti_idx)
    anti_idx = np.array(anti_idx)
    anti_df = df[anti_idx]
    nano_df = df[nano_idx]    
    return anti_df, nano_df

sample_anti_df, sample_nano_df = split_get_imp_result(sample_imp_result, anti_json)


results = []
acc_succ_count = 0
cal_num = 0
for pdb_name, group in sample_anti_df.groupby("PDB_Name"):
    # print(pdb_name)
    if pdb_name in list(ground_truth_result['PDB_Name']):  # Ensure we have a reference value
        cal_num += 1
        success_count = count_success(group, get_ref_dG(ground_truth_result, pdb_name))
        results.append({"PDB_Name": pdb_name, "success_imp": success_count})
        acc_succ_count += success_count
N = sample_anti_df['PDB_Name'].nunique()  
p = acc_succ_count / len(sample_anti_df)
se = np.sqrt(p * (1 - p) / N)
print(f"IMP (anti): {p:.4f} ± {se:.4f}")

results = []
acc_succ_count = 0
cal_num = 0
for pdb_name, group in sample_nano_df.groupby("PDB_Name"):
    # print(pdb_name)
    if pdb_name in list(ground_truth_result['PDB_Name']):  # Ensure we have a reference value
        cal_num += 1
        success_count = count_success(group, get_ref_dG(ground_truth_result, pdb_name))
        results.append({"PDB_Name": pdb_name, "success_imp": success_count})
        acc_succ_count += success_count
N_nano = sample_nano_df['PDB_Name'].nunique()
p_nano = acc_succ_count / len(sample_nano_df)  
se_nano = np.sqrt(p_nano * (1 - p_nano) / N_nano)
print(f"IMP (nano): {p_nano:.4f} ± {se_nano:.4f}")

# Calculate IMP (all)
acc_succ_count_all = 0
for pdb_name, group in sample_imp_result.groupby("PDB_Name"):
    if pdb_name in list(ground_truth_result['PDB_Name']):
        acc_succ_count_all += count_success(group, get_ref_dG(ground_truth_result, pdb_name))
N_all = sample_imp_result['PDB_Name'].nunique()
p_all = acc_succ_count_all / len(sample_imp_result)
se_all = np.sqrt(p_all * (1 - p_all) / N_all)
print(f"IMP (all): {p_all:.4f} ± {se_all:.4f}")


