sbatch -p yejin-lo run_antibmpnn_batch.sh --pdb_dir ../mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.0 --pdb_suffix _model_0.pdb --sampling_temp 0.2 --output_dir ./v2_antibmpnn_mfdesign_output/bn0.0_raw_t0.2

sbatch -p yejin run_antibmpnn_batch.sh --relaxed --pdb_dir ../mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.0 --pdb_suffix _model_0.pdb --sampling_temp 0.2 --output_dir ./v2_antibmpnn_mfdesign_output/relaxed_bn0.0_raw_t0.2

sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir ../mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.0 --pdb_suffix _model_0.pdb --sampling_temp 0.1 --output_dir ./v2_antibmpnn_mfdesign_output/relaxed_bn0.0_raw_t0.1

sbatch -p yejin run_antibmpnn_batch.sh --relaxed --pdb_dir ../mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.0 --pdb_suffix _model_0.pdb --sampling_temp 0.1 --output_dir ./v2_antibmpnn_mfdesign_output/relaxed_bn0.0_raw_t0.1


sbatch -p yejin-lo run_antibmpnn_batch.sh --pdb_dir ../gen_trial_1/predictions --backbone_noise 0.0 --pdb_suffix _model_0.pdb --sampling_temp 0.2 --output_dir ./v3_antibmpnn_output/bn0.0_raw_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --relaxed --pdb_dir ../gen_trial_1/predictions --backbone_noise 0.0 --pdb_suffix _model_0.pdb --sampling_temp 0.2 --output_dir ./v3_antibmpnn_output/relaxed_bn0.0_raw_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir ../gen_trial_1/predictions --backbone_noise 0.0 --pdb_suffix _model_0.pdb --sampling_temp 0.1 --output_dir ./v3_antibmpnn_output/relaxed_bn0.0_raw_t0.1
sbatch -p yejin run_antibmpnn_batch.sh --relaxed --pdb_dir ../gen_trial_1/predictions --backbone_noise 0.0 --pdb_suffix _model_0.pdb --sampling_temp 0.1 --output_dir ./v3_antibmpnn_output/relaxed_bn0.0_raw_t0.1
