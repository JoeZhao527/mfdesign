#!/bin/bash
#
# Run all 6 evaluation runs: 3 design models x 2 structure types.
#
# Reads (via eval_metric.py):
#   {DESIGN_BASE}/{model}/predictions/*_model_0.pdb     (design structures)
#   {PROTENIX_PDB}/{model}/predictions/*_sample_0.pdb   (protenix structures)
#   data/test_entry_pdb_files/*_reference.pdb            (references)
#   data/test_yaml_dir/ab/*.yaml                         (CDR masks)
#   data/test_entry.json                                 (target list)
#
# Writes:
#   {DESIGN_BASE}/{model}/predictions/results.csv        (3 files)
#   {PROTENIX_PDB}/{model}/predictions/results.csv       (3 files)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

FAIL_COUNT=0

DESIGN_BASE="/hai/scratch/fangwu97/protein_data_3stages_v3/r1_exp_data/r1_design_structures"
PROTENIX_PDB="/hai/scratch/fangwu97/protein_data_3stages_v3/r1_protenix_pdb"

MODELS=(
    "gen_trial_0"
    "gen_trial_1"
    "mfdesign_test_bagel_protein_boltz_base"
)

COMMON_ARGS=(
    --ref_dir ./data/test_entry_pdb_files
    --test_yaml_dir data/test_yaml_dir/ab
    -test_json_fpath data/test_entry.json
    --model boltz
    -c 40
)

echo "======================================"
echo "  R1 Evaluation: 6 runs"
echo "======================================"

# Runs 1-3: Design structures (unrelaxed)
for model in "${MODELS[@]}"; do
    data_dir="${DESIGN_BASE}/${model}/predictions"
    echo ""
    echo "--- [Design] ${model} ---"
    echo "  data_dir: ${data_dir}"
    echo "  suffix:   _model_0.pdb"
    python evaluate/AbX_eval/eval_metric.py \
        --data_dir "${data_dir}" \
        --suffix _model_0.pdb \
        "${COMMON_ARGS[@]}" || { echo "  FAILED: ${model}"; FAIL_COUNT=$((FAIL_COUNT+1)); }
done

# Runs 4-6: Protenix folded (raw)
for model in "${MODELS[@]}"; do
    data_dir="${PROTENIX_PDB}/${model}/predictions"
    echo ""
    echo "--- [Protenix] ${model} ---"
    echo "  data_dir: ${data_dir}"
    echo "  suffix:   _sample_0.pdb"
    python evaluate/AbX_eval/eval_metric.py \
        --data_dir "${data_dir}" \
        --suffix _sample_0.pdb \
        "${COMMON_ARGS[@]}" || { echo "  FAILED: ${model}"; FAIL_COUNT=$((FAIL_COUNT+1)); }
done

echo ""
echo "======================================"
echo "  All 6 runs complete. Failures: ${FAIL_COUNT}"
echo "======================================"
exit ${FAIL_COUNT}
