#!/bin/bash
#SBATCH --job-name=inpaint
#SBATCH --output=logs/inpaint.log
#SBATCH --chdir=/hai/scratch/fangwu97/mfdesign
#SBATCH --cpus-per-task=4        # Number of CPU cores per task
#SBATCH --mem=100G                # Total memory per node (or job), e.g. 16G, 3200M
#SBATCH --time=1-00:00:00       # 3 days
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=yejin

CONDA=/hai/scratch/fangwu97/miniconda3

source $CONDA/bin/activate

conda activate $CONDA/envs/mfdesign

export PYTHONPATH=$CONDA/envs/mfdesign/bin/python

# Get data directory from command line argument
DATA_DIR=$1
OUTPUT_DIR=$2

# Run inpaint
python scripts/predict.py \
    --data $DATA_DIR \
    --out_dir $OUTPUT_DIR \
    --processed_msa_dir data/msa \
    --checkpoint ./model/boltz1.ckpt \
    --only_structure_prediction \
    --structure_inpainting

