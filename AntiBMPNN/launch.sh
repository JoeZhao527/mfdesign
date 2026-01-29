sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.1 --pdb_suffix _model_0.pdb --sampling_temp 0.2 --output_dir ./antibmpnn_output/bn0.1_raw_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.1 --pdb_suffix _model_0.pdb --sampling_temp 0.4 --output_dir ./antibmpnn_output/bn0.1_raw_t0.4
sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.1 --relaxed --sampling_temp 0.2 --output_dir ./antibmpnn_output/bn0.1_relaxed_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.1 --relaxed --sampling_temp 0.4 --output_dir ./antibmpnn_output/bn0.1_relaxed_t0.4
sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.5 --pdb_suffix _model_0.pdb --sampling_temp 0.2 --output_dir ./antibmpnn_output/bn0.5_raw_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.5 --pdb_suffix _model_0.pdb --sampling_temp 0.4 --output_dir ./antibmpnn_output/bn0.5_raw_t0.4
sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.5 --relaxed --sampling_temp 0.2 --output_dir ./antibmpnn_output/bn0.5_relaxed_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.5 --relaxed --sampling_temp 0.4 --output_dir ./antibmpnn_output/bn0.5_relaxed_t0.4


sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.1 --pdb_suffix _model_0.pdb --sampling_temp 0.2 --output_dir ./antibmpnn_output/bn0.1_raw_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --backbone_noise 0.1 --pdb_suffix _model_0.pdb --sampling_temp 0.4 --output_dir ./antibmpnn_output/bn0.1_raw_t0.4

sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.1 --pdb_suffix _model_0.pdb --sampling_temp 0.2 --output_dir ./antibmpnn_mfdesign_output/bn0.1_raw_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.1 --pdb_suffix _model_0.pdb --sampling_temp 0.4 --output_dir ./antibmpnn_mfdesign_output/bn0.1_raw_t0.4
sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.1 --relaxed --sampling_temp 0.2 --output_dir ./antibmpnn_mfdesign_output/bn0.1_relaxed_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.1 --relaxed --sampling_temp 0.4 --output_dir ./antibmpnn_mfdesign_output/bn0.1_relaxed_t0.4
sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.5 --pdb_suffix _model_0.pdb --sampling_temp 0.2 --output_dir ./antibmpnn_mfdesign_output/bn0.5_raw_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.5 --pdb_suffix _model_0.pdb --sampling_temp 0.4 --output_dir ./antibmpnn_mfdesign_output/bn0.5_raw_t0.4
sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.5 --relaxed --sampling_temp 0.2 --output_dir ./antibmpnn_mfdesign_output/bn0.5_relaxed_t0.2
sbatch -p yejin run_antibmpnn_batch.sh --pdb_dir mfdesign_test_bagel_protein_boltz_base/predictions --backbone_noise 0.5 --relaxed --sampling_temp 0.4 --output_dir ./antibmpnn_mfdesign_output/bn0.5_relaxed_t0.4
