sbatch run_antibmpnn_batch.sh --backbone_noise 0.1 --pdb_suffix _model_0.pdb --sampling_temp 0.2
sbatch run_antibmpnn_batch.sh --backbone_noise 0.1 --pdb_suffix _model_0.pdb --sampling_temp 0.4
sbatch run_antibmpnn_batch.sh --backbone_noise 0.1 --relaxed --sampling_temp 0.2
sbatch run_antibmpnn_batch.sh --backbone_noise 0.1 --relaxed --sampling_temp 0.4
sbatch run_antibmpnn_batch.sh --backbone_noise 0.5 --pdb_suffix _model_0.pdb --sampling_temp 0.2
sbatch run_antibmpnn_batch.sh --backbone_noise 0.5 --pdb_suffix _model_0.pdb --sampling_temp 0.4
sbatch run_antibmpnn_batch.sh --backbone_noise 0.5 --relaxed --sampling_temp 0.2
sbatch run_antibmpnn_batch.sh --backbone_noise 0.5 --relaxed --sampling_temp 0.4
