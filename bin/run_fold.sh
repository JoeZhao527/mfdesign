#!/bin/bash
#SBATCH --job-name=fold
#SBATCH --output=logs/fold.log
#SBATCH --chdir=/hai/scratch/fangwu97/mfdesign
#SBATCH --cpus-per-task=30        # Number of CPU cores per task
#SBATCH --mem=628G                # Total memory per node (or job), e.g. 16G, 3200M
#SBATCH --time=1-00:00:00       # 3 days
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=yejin

CONDA=/hai/scratch/fangwu97/miniconda3

source $CONDA/bin/activate

conda activate $CONDA/envs/mfdesign

export PYTHONPATH=$CONDA/envs/mfdesign/bin/python

# Run fold
python scripts/predict.py \
    --data data/structure/fold \
    --processed_msa_dir data/msa \
    --checkpoint ./model/boltz1.ckpt \
    --only_structure_prediction

