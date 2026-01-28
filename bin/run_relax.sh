#!/bin/bash

python evaluate/AbX_eval/our_relax_pdb.py \
        --data_dir native_output/fold/boltz_results_fold/predictions/ \
        --output_dir native_output/fold/boltz_results_fold/predictions/ \
        --split_json data/test_entry.json \
        --type boltz \
        --ref_dir ./data/test_entry_pdb_files \
        --test_yaml_dir data/test_yaml_dir/ab \
        --scheme chothia \
        --cpus 20 

python evaluate/AbX_eval/our_relax_pdb.py \
        --data_dir native_output/inpaint/boltz_results_inpaint/predictions/ \
        --output_dir native_output/inpaint/boltz_results_inpaint/predictions/ \
        --split_json data/test_entry.json \
        --type boltz \
        --ref_dir ./data/test_entry_pdb_files \
        --test_yaml_dir data/test_yaml_dir/ab \
        --scheme chothia \
        --cpus 20 

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

# Relax fold results
for dir in "${FOLD_DIRS[@]}"; do
    echo "Relaxing fold: $dir"
    python evaluate/AbX_eval/our_relax_pdb.py \
        --data_dir output/fold/$dir/predictions \
        --output_dir output/fold/$dir/relaxed \
        --split_json data/test_entry.json \
        --type boltz \
        --ref_dir ./data/test_entry_pdb_files \
        --test_yaml_dir data/test_yaml_dir/ab \
        --scheme chothia \
        --cpus 20
done

# Relax inpaint results
for dir in "${INPAINT_DIRS[@]}"; do
    echo "Relaxing inpaint: $dir"
    python evaluate/AbX_eval/our_relax_pdb.py \
        --data_dir output/inpaint/$dir/predictions \
        --output_dir output/inpaint/$dir/relaxed \
        --split_json data/test_entry.json \
        --type boltz \
        --ref_dir ./data/test_entry_pdb_files \
        --test_yaml_dir data/test_yaml_dir/ab \
        --scheme chothia \
        --cpus 20
done
