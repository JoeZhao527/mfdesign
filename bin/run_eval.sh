#!/bin/bash

# Fold results
FOLD_DIRS=(
    "boltz_results_boltz_all_mask_no_playback_0"
    "boltz_results_boltz_all_mask_no_playback_1"
    "boltz_results_boltz_all_mask_playback_0"
    "boltz_results_boltz_all_mask_playback_1"
)

# Inpaint results
INPAINT_DIRS=(
    "boltz_results_boltz_predicted_no_playback_0"
    "boltz_results_boltz_predicted_no_playback_1"
    "boltz_results_boltz_predicted_playback_0"
    "boltz_results_boltz_predicted_playback_1"
)

# Evaluate fold results
for dir in "${FOLD_DIRS[@]}"; do
    echo "Evaluating fold: $dir"
    python evaluate/AbX_eval/eval_metric.py \
        --data_dir relaxed_out/fold/$dir/predictions \
        --ref_dir ./data/test_entry_pdb_files \
        --test_yaml_dir data/test_yaml_dir/ab \
        -test_json_fpath data/test_entry.json \
        --model boltz \
        --suffix _model_0_relaxed.pdb \
        -c 40
done

# Evaluate inpaint results
for dir in "${INPAINT_DIRS[@]}"; do
    echo "Evaluating inpaint: $dir"
    python evaluate/AbX_eval/eval_metric.py \
        --data_dir relaxed_out/inpaint/$dir/predictions \
        --ref_dir ./data/test_entry_pdb_files \
        --test_yaml_dir data/test_yaml_dir/ab \
        -test_json_fpath data/test_entry.json \
        --model boltz \
        --suffix _model_0_relaxed.pdb \
        -c 40
done


python evaluate/AbX_eval/eval_metric.py \
        --data_dir 0128_out/inpaint/boltz_results_boltz_predicted_trail_1/predictions \
        --ref_dir ./data/test_entry_pdb_files \
        --test_yaml_dir data/test_yaml_dir/ab \
        -test_json_fpath data/test_entry.json \
        --model boltz \
        --suffix _model_0.pdb \
        -c 40


python evaluate/AbX_eval/eval_metric.py \
    --data_dir 0128_out/inpaint/boltz_results_boltz_predicted_trail_1/predictions \
    --ref_dir ./data/test_entry_pdb_files \
    --test_yaml_dir data/test_yaml_dir/ab \
    -test_json_fpath evaluate/AbX_eval/regular_list.json \
    --model boltz \
    --suffix _model_0.pdb \
    -c 10

python evaluate/AbX_eval/eval_metric.py \
    --data_dir 0128_out/inpaint/boltz_results_boltz_predicted_trail_1/predictions \
    --ref_dir ./data/test_entry_pdb_files \
    --test_yaml_dir data/test_yaml_dir/ab \
    -test_json_fpath evaluate/AbX_eval/nano_list.json \
    --model boltz \
    --suffix _model_0.pdb \
    -c 10