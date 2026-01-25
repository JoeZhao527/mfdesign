python evaluate/AbX_eval/eval_metric.py \
    --data_dir boltz_results_fold/predictions \
    --ref_dir ./raw_data/chothia \
    --test_yaml_dir data/test_yaml_dir/ab \
    -test_json_fpath data/test_entry.json \
    --model boltz \
    -c 8