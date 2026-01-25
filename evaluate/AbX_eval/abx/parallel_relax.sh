#!/bin/bash

# Input and output directories
DATA_DIR="/mnt/nas-new/home/yangnianzu/icml/baseline/dyMEAN/results/"
PYTHON_SCRIPT="/mnt/nas-new/home/yangnianzu/icml/baseline/AbX/abx/single_relax.py" # 替换为你的 Python 脚本路径
MAX_JOBS=24 # 最大并行任务数

# 获取所有后缀为 _filter.pdb 的文件
filter_pdb_files=($(find "$DATA_DIR" -maxdepth 1 -name "*_filter.pdb"))

# 定义函数，控制并行任务数
run_task() {
    local pdb_file=$1
    local pdb_name=$(basename "$pdb_file" "_filter.pdb") # 提取原始名称

    # 调用 Python 脚本
    python "$PYTHON_SCRIPT" --pdb "$pdb_file" --pdb-name "$pdb_name" &
}

# 遍历所有 filter 文件并并行处理
for filter_pdb_file in "${filter_pdb_files[@]}"; do
    # 限制并行任务数
    while [[ $(jobs -r | wc -l) -ge $MAX_JOBS ]]; do
        sleep 1
    done

    # 调用任务
    run_task "$filter_pdb_file"
done

# 等待所有后台任务完成
wait

echo "All tasks completed."

