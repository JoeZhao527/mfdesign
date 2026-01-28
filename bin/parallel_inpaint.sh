#!/bin/bash
# Parallel inpaint launcher script
# Usage: ./bin/parallel_inpaint.sh <data-dir> <num-partitions> [partition]
# Example: ./bin/parallel_inpaint.sh final_result_tools/output/boltz_predicted_trail_1 4 yejin

set -e

DATA_DIR=$1
NUM_PARTITIONS=${2:-3}
PARTITION=${3:-yejin}

if [ -z "$DATA_DIR" ]; then
    echo "Usage: $0 <data-dir> [num-partitions] [slurm-partition]"
    echo "Example: $0 final_result_tools/output/boltz_predicted_trail_1 4 yejin"
    exit 1
fi

# Remove trailing slash if present
DATA_DIR=${DATA_DIR%/}

# Create tmp directory
TMP_DIR="${DATA_DIR}_tmp"
mkdir -p "$TMP_DIR"

# Get all yaml files
YAML_FILES=($(ls "$DATA_DIR"/*.yaml 2>/dev/null))
TOTAL_FILES=${#YAML_FILES[@]}

if [ "$TOTAL_FILES" -eq 0 ]; then
    echo "No yaml files found in $DATA_DIR"
    exit 1
fi

echo "Found $TOTAL_FILES yaml files in $DATA_DIR"
echo "Partitioning into $NUM_PARTITIONS sub-folders..."

# Calculate files per partition
FILES_PER_PARTITION=$(( (TOTAL_FILES + NUM_PARTITIONS - 1) / NUM_PARTITIONS ))

# Create sub-folders and distribute files
for ((i=0; i<NUM_PARTITIONS; i++)); do
    SUB_DIR="${TMP_DIR}/part_${i}"
    mkdir -p "$SUB_DIR"
    
    START_IDX=$((i * FILES_PER_PARTITION))
    END_IDX=$((START_IDX + FILES_PER_PARTITION))
    
    # Copy files to sub-folder
    for ((j=START_IDX; j<END_IDX && j<TOTAL_FILES; j++)); do
        cp "${YAML_FILES[$j]}" "$SUB_DIR/"
    done
    
    # Count files in this partition
    FILE_COUNT=$(ls "$SUB_DIR"/*.yaml 2>/dev/null | wc -l)
    echo "Partition $i: $FILE_COUNT files"
done

# Create logs directory if it doesn't exist
mkdir -p logs

# Extract base name for log files
BASE_NAME=$(basename "$DATA_DIR")

# Launch sbatch jobs for each partition
echo ""
echo "Launching sbatch jobs..."
for ((i=0; i<NUM_PARTITIONS; i++)); do
    SUB_DIR="${TMP_DIR}/part_${i}"
    LOG_FILE="logs/${BASE_NAME}_part_${i}_inpaint.log"
    
    OUTPUT_DIR="./0128_inpaint_out/${BASE_NAME}"
    echo "Submitting job for partition $i -> $LOG_FILE (output: $OUTPUT_DIR)"
    sbatch -p "$PARTITION" --output="$LOG_FILE" bin/run_inpaint.sh "$SUB_DIR" "$OUTPUT_DIR"
done

echo ""
echo "All jobs submitted. Check logs in logs/ directory."
echo "Monitor jobs with: squeue -u \$USER"
