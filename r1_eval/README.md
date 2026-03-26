# R1 Evaluation Pipeline

评估 3 个设计模型 × 2 种结构类型（设计直出结构 + Protenix 折叠结构）的 CDR RMSD/AAR 指标。

## 目录结构

```
r1_eval/
├── convert_protenix_cif_to_pdb.py   # Step 1: CIF → PDB 转换
├── run_all_eval.sh                  # Step 2: 运行 6 组评估
├── collect_results.py               # Step 3: 汇总结果
└── README.md
```

## 3 个设计模型

| 模型 | 说明 |
|---|---|
| `gen_trial_0` | 生成试验 0 |
| `gen_trial_1` | 生成试验 1 |
| `mfdesign_test_bagel_protein_boltz_base` | MFDesign bagel boltz base |

## 数据路径

### 输入

| 数据 | 路径 |
|---|---|
| 设计直出结构 (PDB) | `/hai/scratch/fangwu97/protein_data_3stages_v3/r1_exp_data/r1_design_structures/{model}/predictions/{ID}/{ID}_model_0.pdb` |
| Protenix 折叠结构 (CIF) | `/hai/scratch/fangwu97/protein_data_3stages_v3/r1_protenix/{model}/prediction/{ID}/seed_101/predictions/{ID}_sample_0.cif` |
| 参考结构 | `data/test_entry_pdb_files/{ID}_reference.pdb` |
| CDR mask | `data/test_yaml_dir/ab/{ID}.yaml` |
| 目标列表 | `data/test_entry.json` |

### 输出

| 数据 | 路径 |
|---|---|
| 转换后 Protenix PDB | `/hai/scratch/fangwu97/protein_data_3stages_v3/r1_protenix_pdb/{model}/predictions/{ID}/{ID}_sample_0.pdb` |
| 设计评估结果 | `/hai/scratch/fangwu97/protein_data_3stages_v3/r1_exp_data/r1_design_structures/{model}/predictions/results.csv` |
| Protenix 评估结果 | `/hai/scratch/fangwu97/protein_data_3stages_v3/r1_protenix_pdb/{model}/predictions/results.csv` |
| 汇总对比表 | `r1_eval/r1_comparison.csv`, `r1_eval/r1_comparison_full.csv` |

## 运行方式

```bash
# 在项目根目录 (/hai/scratch/fangwu97/mfdesign) 下运行

# Step 1: 将 Protenix CIF 转为 PDB（自动验证 chain 映射）
python r1_eval/convert_protenix_cif_to_pdb.py

# Step 2: 运行 6 组评估
bash r1_eval/run_all_eval.sh

# Step 3: 汇总结果
python r1_eval/collect_results.py
```

## 各脚本说明

### convert_protenix_cif_to_pdb.py

将 Protenix 输出的 CIF 文件转为 PDB 格式，同时将嵌套目录结构（`{ID}/seed_101/predictions/`）展平为 `eval_metric.py` 期望的一层结构（`{ID}/`）。

- 仅转换 `sample_0`
- 运行前自动校验 chain 映射（CIF chain A = heavy chain）
- 40 进程并行，支持断点续转（已存在的 PDB 自动跳过）
- `--skip-verify` 可跳过校验

### run_all_eval.sh

调用 `evaluate/AbX_eval/eval_metric.py` 依次运行 6 组评估：

| Run | 类型 | suffix |
|---|---|---|
| 1-3 | 设计直出 (unrelaxed) | `_model_0.pdb` |
| 4-6 | Protenix 折叠 (raw) | `_sample_0.pdb` |

所有 run 使用 `--model boltz`（chain A=heavy, B=light）、`test_entry.json` 全集、40 核并行。

### collect_results.py

读取 6 个 `results.csv`，按 code 分组计算 mean +/- SEM，输出对比表：

```
                         | heavy_cdr3_RMSD | heavy_cdr3_AAR | ... | full_RMSD
gen_trial_0 (design)     |                 |                |     |
gen_trial_0 (protenix)   |                 |                |     |
gen_trial_1 (design)     |                 |                |     |
gen_trial_1 (protenix)   |                 |                |     |
mfdesign_bagel (design)  |                 |                |     |
mfdesign_bagel (protenix)|                 |                |     |
```

## 评估指标

由 `eval_metric.py` 计算（基于 framework 对齐后的 CA 坐标）：

- **RMSD**: 各 CDR 区域的均方根偏差
- **AAR**: 氨基酸恢复率（设计序列 vs 天然序列）
- **Loop_RMSD / Loop_AAR**: CDR3 loop 核心区域（去掉两端 flanking 残基）
- **full_RMSD**: 全抗体 RMSD
