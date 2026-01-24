#!/bin/bash
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

data_dir=$1

# data/structure/inpaint
# data/structure/fold
# Run inpaint
python scripts/predict.py --data ${data_dir} --processed_msa_dir data/msa
