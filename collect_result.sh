#!/bin/bash
# Collect results from AntiBMPNN experiments and calculate CDR AAR

set -e

YAML_DIR="final_result_tools/boltz_data"
OUTPUT_BASE="./AntiBMPNN/v3_antibmpnn_output_metrics"

# Find all fasta directories using find
for fasta_dir in $(find ./AntiBMPNN/v3_antibmpnn_output -type d -name "seqs" 2>/dev/null); do
    
    # Structure: ./antibmpnn_output/{exp_name}/{date}_batch_design/seqs
    # Extract date_theme (e.g., 0129_batch_design)
    date_theme=$(basename $(dirname "$fasta_dir"))
    # Extract experiment name (e.g., bn0.1_raw_t0.2)
    exp_name=$(basename $(dirname $(dirname "$fasta_dir")))
    # Extract base name (antibmpnn_output or antibmpnn_mfdesign_output)
    base_name=$(basename $(dirname $(dirname $(dirname "$fasta_dir"))))
    
    output_dir="${OUTPUT_BASE}/${base_name}/${exp_name}"
    
    echo "=========================================="
    echo "Processing: $fasta_dir"
    echo "Output to: $output_dir"
    echo "=========================================="
    
    # Step 1: Update YAML from FASTA
    python final_result_tools/update_yaml_from_fasta.py \
        --fasta_dir "$fasta_dir" \
        --yaml_dir "$YAML_DIR" \
        --output_dir "$output_dir"
    
    # Step 2: Calculate CDR AAR for regular antibodies
    echo ""
    echo "--- Regular Antibodies ---"
    python final_result_tools/calc_cdr_aar.py \
        --yaml_dir "$output_dir" \
        --list_file evaluate/AbX_eval/regular_list.json
    
    # Step 3: Calculate CDR AAR for nanobodies
    echo ""
    echo "--- Nanobodies ---"
    python final_result_tools/calc_cdr_aar.py \
        --yaml_dir "$output_dir" \
        --list_file evaluate/AbX_eval/nano_list.json
    
    echo ""
done

echo "Done!"
