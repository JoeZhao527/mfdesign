python evaluate/AbX_eval/cal_interface_energy.py \
    --pdb_dir relaxed_out/trial0/predictions \
    --output_dir binding_res/raw_gen_trial_0 \
    --cpus 50 &

python evaluate/AbX_eval/cal_interface_energy.py \
    --pdb_dir relaxed_out/trial1/predictions \
    --output_dir binding_res/raw_gen_trial_1 \
    --cpus 50 &

wait

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir relaxed_out/fold/boltz_results_boltz_all_mask_playback_0/predictions \
#     --output_dir binding_res/pb_fold_0 \
#     --cpus 30

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir relaxed_out/fold/boltz_results_boltz_all_mask_playback_1/predictions \
#     --output_dir binding_res/pb_fold_1 \
#     --cpus 30

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir relaxed_out/fold/boltz_results_boltz_all_mask_no_playback_0/predictions \
#     --output_dir binding_res/no_pb_fold_0 \
#     --cpus 30

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir relaxed_out/fold/boltz_results_boltz_all_mask_no_playback_1/predictions \
#     --output_dir binding_res/no_pb_fold_1 \
#     --cpus 30

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir relaxed_out/inpaint/boltz_results_boltz_predicted_playback_0/predictions \
#     --output_dir binding_res/pb_inpaint_0 \
#     --cpus 30

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir relaxed_out/inpaint/boltz_results_boltz_predicted_playback_1/predictions \
#     --output_dir binding_res/pb_inpaint_1 \
#     --cpus 30

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir relaxed_out/inpaint/boltz_results_boltz_predicted_no_playback_0/predictions \
#     --output_dir binding_res/no_pb_inpaint_0 \
#     --cpus 30

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir relaxed_out/inpaint/boltz_results_boltz_predicted_no_playback_1/predictions \
#     --output_dir binding_res/no_pb_inpaint_1 \
#     --cpus 30

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir native_output/inpaint/boltz_results_inpaint/predictions \
#     --output_dir binding_res/boltz_inpaint \
#     --cpus 30

# python evaluate/AbX_eval/cal_interface_energy.py \
#     --pdb_dir native_output/fold/boltz_results_fold/predictions \
#     --output_dir binding_res/boltz_fold \
#     --cpus 30