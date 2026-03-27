#!/bin/bash
#SBATCH --job-name=msa
#SBATCH --output=logs/msa.log
#SBATCH --chdir=/hai/scratch/fangwu97/mfdesign
#SBATCH --cpus-per-task=16        # Number of CPU cores per task
#SBATCH --mem=100G                # Total memory per node (or job), e.g. 16G, 3200M
#SBATCH --time=1-00:00:00       # 3 days
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=yejin

CONDA=/hai/scratch/fangwu97/miniconda3

source $CONDA/bin/activate

conda activate $CONDA/envs/mfdesign

export PYTHONPATH=$CONDA/envs/mfdesign/bin/python

python scripts/predict.py \
    --data rabd_eval/mfdesign_input/yamls/ \
    --out_dir rabd_eval/mfdesign_input/msa/ \
    --use_msa_server \
    --msa_server_url https://api.colabfold.com \
    --msa_pairing_strategy greedy \
    --only_process_msa