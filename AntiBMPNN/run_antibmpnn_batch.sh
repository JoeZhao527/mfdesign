#!/bin/bash
#SBATCH --job-name=fold
#SBATCH --output=logs/fold.log
#SBATCH --chdir=/hai/scratch/fangwu97/mfdesign/AntiBMPNN
#SBATCH --cpus-per-task=20        # Number of CPU cores per task
#SBATCH --mem=200G                # Total memory per node (or job), e.g. 16G, 3200M
#SBATCH --time=1-00:00:00       # 3 days
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=yejin

#
# Batch script for running AntiBMPNN on predicted structures
#
# Usage:
#   bash run_antibmpnn_batch.sh [options]
#
# This script:
# 1. Calls prepare_antibmpnn_input.py to prepare input data from YAML and PDB files
# 2. Runs AntiBMPNN for inverse folding on the prepared data
#

set -e  # Exit on error

# ==================== Configuration ====================
# Paths (modify these according to your setup)
YAML_DIR="../final_result_tools/boltz_data"
PDB_DIR="../0128_out/inpaint/boltz_results_boltz_predicted_trail_0/predictions"
OUTPUT_BASE="./antibmpnn_output"
THEME=$(date +"%m%d")_"batch_design"

# AntiBMPNN parameters
MODEL_NAME="antibmpnn_000"
NUM_SEQ_PER_TARGET=20
SAMPLING_TEMP="0.2"
BATCH_SIZE=10
BACKBONE_NOISE=0

# PDB suffix (use _model_0.pdb for non-relaxed, _model_0_relaxed.pdb for relaxed)
PDB_SUFFIX="_model_0.pdb"

# ==================== Parse arguments ====================
while [[ $# -gt 0 ]]; do
    case $1 in
        --yaml_dir)
            YAML_DIR="$2"
            shift 2
            ;;
        --pdb_dir)
            PDB_DIR="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_BASE="$2"
            shift 2
            ;;
        --pdb_suffix)
            PDB_SUFFIX="$2"
            shift 2
            ;;
        --model_name)
            MODEL_NAME="$2"
            shift 2
            ;;
        --num_seq)
            NUM_SEQ_PER_TARGET="$2"
            shift 2
            ;;
        --sampling_temp)
            SAMPLING_TEMP="$2"
            shift 2
            ;;
        --batch_size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --backbone_noise)
            BACKBONE_NOISE="$2"
            shift 2
            ;;
        --relaxed)
            PDB_SUFFIX="_model_0_relaxed.pdb"
            shift
            ;;
        --help)
            echo "Usage: bash run_antibmpnn_batch.sh [options]"
            echo ""
            echo "Options:"
            echo "  --yaml_dir DIR       Directory containing YAML files (default: ../final_result_tools/boltz_data)"
            echo "  --pdb_dir DIR        Directory containing predicted PDBs (default: ../0128_out/inpaint/boltz_results_boltz_predicted_trail_0/predictions)"
            echo "  --output_dir DIR     Output directory (default: ./antibmpnn_output)"
            echo "  --pdb_suffix SUFFIX  PDB file suffix (default: _model_0.pdb)"
            echo "  --relaxed            Use relaxed PDB files (_model_0_relaxed.pdb)"
            echo "  --model_name NAME    Model name (default: antibmpnn_000)"
            echo "  --num_seq N          Number of sequences per target (default: 1000)"
            echo "  --sampling_temp T    Sampling temperature (default: 0.1)"
            echo "  --batch_size N       Batch size (default: 10)"
            echo "  --help               Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Output directory with theme
OUTPUT_DIR="${OUTPUT_BASE}/${THEME}"

# ==================== Load environment ====================
echo "=========================================="
echo "AntiBMPNN Batch Inference Script"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  YAML Dir:      ${YAML_DIR}"
echo "  PDB Dir:       ${PDB_DIR}"
echo "  Output Dir:    ${OUTPUT_DIR}"
echo "  PDB Suffix:    ${PDB_SUFFIX}"
echo "  Model:         ${MODEL_NAME}"
echo "  Num Seqs:      ${NUM_SEQ_PER_TARGET}"
echo "  Temperature:   ${SAMPLING_TEMP}"
echo "  Batch Size:    ${BATCH_SIZE}"
echo ""

# Activate conda environment
source activate mlfold 2>/dev/null || conda activate mlfold 2>/dev/null || {
    echo "Warning: Could not activate mlfold environment"
    echo "Make sure conda environment is set up correctly"
}

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# ==================== Step 1: Prepare input ====================
echo "=========================================="
echo "Step 1: Preparing AntiBMPNN input..."
echo "=========================================="

python prepare_antibmpnn_input.py \
    --yaml_dir "${YAML_DIR}" \
    --pdb_dir "${PDB_DIR}" \
    --output_dir "${OUTPUT_DIR}" \
    --pdb_suffix "${PDB_SUFFIX}"

# Define paths to generated files
PATH_PARSED_CHAINS="${OUTPUT_DIR}/parsed_pdbs.jsonl"
PATH_ASSIGNED_CHAINS="${OUTPUT_DIR}/assigned_pdbs.jsonl"
PATH_FIXED_POSITIONS="${OUTPUT_DIR}/fixed_pdbs.jsonl"

# Check if files were created
if [ ! -f "${PATH_PARSED_CHAINS}" ]; then
    echo "Error: parsed_pdbs.jsonl was not created"
    exit 1
fi

# ==================== Step 2: Run AntiBMPNN ====================
echo ""
echo "=========================================="
echo "Step 2: Running AntiBMPNN..."
echo "=========================================="

python Running_AntiBMPNN_run.py \
    --jsonl_path "${PATH_PARSED_CHAINS}" \
    --chain_id_jsonl "${PATH_ASSIGNED_CHAINS}" \
    --fixed_positions_jsonl "${PATH_FIXED_POSITIONS}" \
    --out_folder "${OUTPUT_DIR}" \
    --model_name "${MODEL_NAME}" \
    --num_seq_per_target ${NUM_SEQ_PER_TARGET} \
    --sampling_temp "${SAMPLING_TEMP}" \
    --batch_size ${BATCH_SIZE} \
    --backbone_noise ${BACKBONE_NOISE}

# ==================== Step 3: Parse results (optional) ====================
echo ""
echo "=========================================="
echo "Step 3: Design complete!"
echo "=========================================="
echo ""
echo "Output files are in: ${OUTPUT_DIR}"
echo "  - seqs/: FASTA files with designed sequences"
echo ""
echo "To parse results for a specific position range, use:"
echo "  python helper_scripts/parsing_design_result.py --dir ${OUTPUT_DIR}/seqs/ --start START --end END"
echo ""
echo "Done!"

