# rabd_prepare

Scripts for preparing RAbD (Rosetta Antibody Design) benchmark data — from raw MMAP binary through to mfdesign-ready NPZ files.

## Input Data

`rabd/` directory containing:
- `data.bin` — MMAP binary file with compressed Complex structures
- `test_index.txt` (60 entries), `train_index.txt` (9473), `valid_index.txt` (400) — tab-separated index with `sample_id`, byte offsets, and properties JSON
- `heavy_chain_mark` in properties uses digits `1/2/3` to denote CDR1/CDR2/CDR3 positions

---

## Pipeline Overview

```
rabd/ (MMAP)
  │
  ▼  [Step 1] convert_rabd_to_jsonl.py --prepare
rabd_eval/cifs/ + rabd_h3_test.jsonl + cdr_info/
  │
  ▼  [Step 2] Protenix prediction + remap (external)
rabd_eval/predictions_remap/
  │
  ▼  [Step 3] convert_rabd_to_jsonl.py --finalize
rabd_eval/rabd_h3_test.jsonl (updated source_cif)
  │
  ▼  [Step 4] prepare_rabd_for_mfdesign.py
rabd_eval/mfdesign_input/{csv, yamls/, pdbs/, summary.json}
  │
  ▼  [Step 5] MSA generation (ColabFold API)
rabd_eval/mfdesign_input/msa/
  │
  ▼  [Step 6] convert_msa.py → MSA NPZ
  ▼  [Step 7] antibody.py   → structure NPZ
rabd_eval/mfdesign_input/npz_output/
```

---

## Step 1: Extract GT Structures (convert_rabd_to_jsonl.py)

```bash
python rabd_prepare/convert_rabd_to_jsonl.py \
    --data_dir rabd/ \
    --output_dir rabd_eval/ \
    --cdr_type H3
```

Outputs:
- `cifs/{pdb_id}.cif` — ground-truth CIF files
- `cdr_info/{SAMPLE_ID}_cdr_info.json` — CDR index info per sample
- `rabd_cdr_eval_split.json` — split list for cdr_masking
- `rabd_h3_test.jsonl` — evaluation JSONL (source_cif = GT)

## Step 2: Protenix Prediction + Remap (external)

Run Protenix structure prediction and chain remapping. Results go to `rabd_eval/predictions/` and `rabd_eval/predictions_remap/`.

## Step 3: Finalize JSONL (convert_rabd_to_jsonl.py --finalize)

Selects best Protenix sample per entry by `ranking_score`:

```bash
python rabd_prepare/convert_rabd_to_jsonl.py \
    --data_dir rabd/ \
    --output_dir rabd_eval/ \
    --finalize \
    --predictions_dir rabd_eval/predictions/ \
    --remap_dir rabd_eval/predictions_remap/ \
    --seed 101
```

Use `--sample_idx N` to force a specific sample index instead of auto-selecting.

## Step 4: Prepare mfdesign Inputs (prepare_rabd_for_mfdesign.py)

Generates the processed CSV, YAML files, PDB files, and summary.json needed by the mfdesign preprocessing scripts:

```bash
python rabd_prepare/prepare_rabd_for_mfdesign.py \
    --data_dir rabd/ \
    --rabd_eval_dir rabd_eval/ \
    --output_dir rabd_eval/mfdesign_input/
```

Outputs:
- `rabd_processed.csv` — processed CSV for `antibody.py` (60 entries)
- `summary.json` — metadata JSON for `convert_msa.py`
- `yamls/{entry}.yaml` — Boltz-compatible YAML entity files
- `pdbs/{pdb_id}.pdb` — PDB structures converted from GT CIF

## Step 5: MSA Generation

Use the ColabFold public API (`https://api.colabfold.com`, no local setup needed):

```bash
python scripts/predict.py \
    --data rabd_eval/mfdesign_input/yamls/ \
    --out_dir rabd_eval/mfdesign_input/msa/ \
    --use_msa_server \
    --msa_server_url https://api.colabfold.com \
    --msa_pairing_strategy greedy \
    --only_process_msa
```

This queries MSA for each chain sequence via the API (~180 queries, ~few hours with rate limits).

## Step 6: Convert MSA to NPZ

```bash
python scripts/process/convert_msa.py \
    --input_dir rabd_eval/mfdesign_input/msa/ \
    --output_dir rabd_eval/mfdesign_input/msa_npz/ \
    --preprocessed_data_path rabd_eval/mfdesign_input/summary.json \
    --msa_filtering_threshold 0.2
```

## Step 7: Generate Structure NPZ

```bash
# Start Redis with CCD data first
redis-server --dbfilename ccd.rdb --port 7777

python scripts/process/antibody.py \
    --datadir rabd_eval/mfdesign_input/pdbs/ \
    --processed_csv_fpath rabd_eval/mfdesign_input/rabd_processed.csv \
    --outdir rabd_eval/mfdesign_input/npz_output/ \
    --cif_path rabd_eval/cifs/ \
    --yaml_path rabd_eval/mfdesign_input/yamls/ \
    --cdr_select H3 \
    --num-processes 8
```

Final output:
- `npz_output/structures/{entry}.npz` — structure NPZ files
- `npz_output/records/{entry}.json` — metadata records
- `npz_output/manifest.json` — index of all records

---

## Dependencies

- numpy, pandas, gemmi, ruamel.yaml (prepare script)
- biotite, protenix (optional, for CIF export in convert_rabd_to_jsonl.py)
- boltz, rdkit, redis (for antibody.py NPZ conversion)
