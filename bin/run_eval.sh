python evaluate/AbX_eval/eval_metric.py \
    --data_dir boltz_results_fold/predictions \
    --ref_dir ./data/test_entry_pdb_files \
    --test_yaml_dir data/test_yaml_dir/ab \
    -test_json_fpath data/test_entry.json \
    --model boltz \
    --suffix _model_0.pdb \
    -c 8 