python evaluate/AbX_eval/cal_interface_energy.py \
    --pdb_dir relaxed_out/inpaint/boltz_results_boltz_predicted_playback_0/predictions \
    --output_dir binding_res/pb_inpaint_0 \
    --cpus 30

python evaluate/AbX_eval/cal_interface_energy.py \
    --pdb_dir relaxed_out/inpaint/boltz_results_boltz_predicted_playback_1/predictions \
    --output_dir binding_res/pb_inpaint_1 \
    --cpus 30

python evaluate/AbX_eval/cal_interface_energy.py \
    --pdb_dir relaxed_out/inpaint/boltz_results_boltz_predicted_no_playback_0/predictions \
    --output_dir binding_res/no_pb_inpaint_0 \
    --cpus 30

python evaluate/AbX_eval/cal_interface_energy.py \
    --pdb_dir relaxed_out/inpaint/boltz_results_boltz_predicted_no_playback_1/predictions \
    --output_dir binding_res/no_pb_inpaint_1 \
    --cpus 30