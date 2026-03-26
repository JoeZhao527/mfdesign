"""
Convert Protenix CIF outputs to PDB format for evaluation.

Reads:
  {PROTENIX_BASE}/{model}/prediction/{ID}/seed_101/predictions/{ID}_sample_0.cif
  data/test_entry_pdb_files/{ID}_reference.pdb  (for chain verification)

Writes:
  {OUTPUT_BASE}/{model}/predictions/{ID}/{ID}_sample_0.pdb
  {OUTPUT_BASE}/conversion_log.json
"""

import os
import json
import argparse
import multiprocessing as mp

from Bio.PDB import MMCIFParser, PDBIO, PDBParser

PROTENIX_BASE = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_protenix"
OUTPUT_BASE = "/hai/scratch/fangwu97/protein_data_3stages_v3/r1_protenix_pdb"
REF_DIR = "/hai/scratch/fangwu97/mfdesign/data/test_entry_pdb_files"

MODELS = [
    "gen_trial_0",
    "gen_trial_1",
    "mfdesign_test_bagel_protein_boltz_base",
]


def verify_chain_mapping(cif_path, ref_path, target_id):
    """Spot-check that CIF chain A = heavy chain by comparing residue counts."""
    cif_parser = MMCIFParser(QUIET=True)
    pdb_parser = PDBParser(QUIET=True)

    cif_struct = cif_parser.get_structure("cif", cif_path)
    ref_struct = pdb_parser.get_structure("ref", ref_path)

    parts = target_id.split("_")
    heavy_id = parts[1]
    light_id = parts[2]

    ref_model = ref_struct[0]
    ref_heavy_len = sum(
        1 for r in ref_model[heavy_id].get_residues() if r.get_id()[0] == " "
    )

    cif_model = cif_struct[0]
    cif_chains = list(cif_model.get_chains())
    chain_info = []
    for ch in cif_chains:
        n_res = sum(1 for r in ch.get_residues() if r.get_id()[0] == " ")
        chain_info.append((ch.id, n_res))

    print(f"[Chain Verification] {target_id}")
    print(f"  Reference heavy chain ({heavy_id}): {ref_heavy_len} residues")
    for ch_id, n_res in chain_info:
        print(f"  CIF chain {ch_id}: {n_res} residues")

    cif_chain_a_len = next((n for cid, n in chain_info if cid == "A"), None)
    if cif_chain_a_len is not None and cif_chain_a_len == ref_heavy_len:
        print("  -> Chain A = heavy chain (lengths match)")
        return True

    if light_id:
        ref_light_len = sum(
            1 for r in ref_model[light_id].get_residues() if r.get_id()[0] == " "
        )
        print(f"  Reference light chain ({light_id}): {ref_light_len} residues")
        if cif_chain_a_len == ref_light_len:
            print("  -> WARNING: Chain A = light chain! Mapping is reversed!")

    print("  -> Chain mapping mismatch - needs investigation")
    return False


def convert_single(args_tuple):
    """Convert a single CIF file to PDB."""
    cif_path, pdb_path = args_tuple
    try:
        os.makedirs(os.path.dirname(pdb_path), exist_ok=True)

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("struct", cif_path)

        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)

        return {"status": "success", "cif": cif_path, "pdb": pdb_path}
    except Exception as e:
        return {"status": "error", "cif": cif_path, "error": str(e)}


def main():
    parser = argparse.ArgumentParser(description="Convert Protenix CIF to PDB")
    parser.add_argument("--cpus", type=int, default=40)
    parser.add_argument("--skip-verify", action="store_true",
                        help="Skip chain mapping verification")
    parser.add_argument("--protenix-base", type=str, default=PROTENIX_BASE)
    parser.add_argument("--output-base", type=str, default=OUTPUT_BASE)
    parser.add_argument("--ref-dir", type=str, default=REF_DIR)
    args = parser.parse_args()

    protenix_base = args.protenix_base
    output_base = args.output_base
    ref_dir = args.ref_dir

    # Step 1: Verify chain mapping on one example
    if not args.skip_verify:
        verified = False
        for model in MODELS:
            pred_dir = os.path.join(protenix_base, model, "prediction")
            if not os.path.exists(pred_dir):
                continue
            for target in sorted(os.listdir(pred_dir)):
                ref_path = os.path.join(ref_dir, f"{target}_reference.pdb")
                cif_path = os.path.join(
                    pred_dir, target, "seed_101", "predictions",
                    f"{target}_sample_0.cif"
                )
                if os.path.exists(ref_path) and os.path.exists(cif_path):
                    verified = verify_chain_mapping(cif_path, ref_path, target)
                    break
            if verified:
                break

        if not verified:
            print("Chain mapping verification FAILED. Use --skip-verify to bypass.")
            return

    # Step 2: Collect conversion tasks
    tasks = []
    for model in MODELS:
        pred_dir = os.path.join(protenix_base, model, "prediction")
        if not os.path.exists(pred_dir):
            print(f"Warning: {pred_dir} not found, skipping")
            continue

        for target_id in sorted(os.listdir(pred_dir)):
            if target_id == "ERR":
                continue
            cif_path = os.path.join(
                pred_dir, target_id, "seed_101", "predictions",
                f"{target_id}_sample_0.cif"
            )
            pdb_path = os.path.join(
                output_base, model, "predictions", target_id,
                f"{target_id}_sample_0.pdb"
            )

            if not os.path.exists(cif_path):
                continue
            if os.path.exists(pdb_path):
                continue  # skip already converted

            tasks.append((cif_path, pdb_path))

    print(f"Converting {len(tasks)} CIF files to PDB...")

    # Step 3: Convert in parallel
    if tasks:
        with mp.Pool(args.cpus) as pool:
            results = list(pool.imap_unordered(convert_single, tasks))
    else:
        results = []
        print("No new files to convert.")

    # Step 4: Log results
    success = sum(1 for r in results if r["status"] == "success")
    errors = [r for r in results if r["status"] == "error"]

    print(f"Done: {success} success, {len(errors)} errors")
    for e in errors:
        print(f"  Error: {e['cif']}: {e['error']}")

    log_path = os.path.join(output_base, "conversion_log.json")
    os.makedirs(output_base, exist_ok=True)
    with open(log_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Log saved to {log_path}")


if __name__ == "__main__":
    main()
