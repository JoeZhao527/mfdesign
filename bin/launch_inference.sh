final_result_tools/output/boltz_all_mask_trail_0
final_result_tools/output/boltz_all_mask_trail_1
final_result_tools/output/boltz_predicted_trail_0
final_result_tools/output/boltz_predicted_trail_1


sbatch -p yejin --output=logs/trail_0_inpaint.log bin/run_inpaint.sh final_result_tools/output/boltz_predicted_trail_0
sbatch -p yejin --output=logs/trail_1_inpaint.log bin/run_inpaint.sh final_result_tools/output/boltz_predicted_trail_1



sbatch -p yejin --output=logs/trail_0_fold.log bin/run_fold.sh final_result_tools/output/boltz_all_mask_trail_0




sbatch -p yejin --output=logs/back_0_fold.log bin/run_fold.sh data/boltz_all_mask_playback_0
sbatch -p yejin --output=logs/back_1_fold.log bin/run_fold.sh data/boltz_all_mask_playback_1
sbatch -p yejin --output=logs/no_back_0_fold.log bin/run_fold.sh data/boltz_all_mask_no_playback_0
sbatch -p yejin --output=logs/no_back_1_fold.log bin/run_fold.sh data/boltz_all_mask_no_playback_1

sbatch -p yejin-lo --output=logs/back_0_inpaint.log bin/run_inpaint.sh data/boltz_predicted_playback_0
sbatch -p yejin-lo --output=logs/back_1_inpaint.log bin/run_inpaint.sh data/boltz_predicted_playback_1
sbatch -p yejin-lo --output=logs/no_back_0_inpaint.log bin/run_inpaint.sh data/boltz_predicted_no_playback_0
sbatch -p yejin-lo --output=logs/no_back_1_inpaint.log bin/run_inpaint.sh data/boltz_predicted_no_playback_1
