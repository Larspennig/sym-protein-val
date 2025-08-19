# Script for running ProteinMPNN and Boltz folding model


import argparse
import os
import subprocess
import torch
from pathlib import Path
from Bio import SeqIO
import yaml
import re
import pandas as pd


def evaluate_results(folder_path):
    sample_paths = list(Path(folder_path).glob("**/confidence*.json"))
    sample_list = []
    for sample in sample_paths:
        sample_name = sample.stem
        # load json
        try:
            df = pd.read_json(sample)
        except ValueError as e:
            print(f"Error reading {sample_name}: {e}")
            continue
        # extract relevant sample. Mean is only computed because value is the same for all chains.
        sample_df = {
            "sample_name": sample_name,
            "num_chains": len(df),
            "confidence_score": df["confidence_score"].mean(),
            "ptm_score": df["ptm_score"].mean(),
            "iptm_score": df["iptm_score"].mean(),
            "complex_plddt": df["complex_plddt"].mean(),
            "complex_iplddt": df["complex_ipl_ddt"].mean(),
            "complex_ipde": df["complex_ipde"].mean()
        }
        sample_list.append(sample_df)

    sample_df = pd.DataFrame(sample_list)
    # save to csv
    output_csv = Path(folder_path) /  "sample_evaluation.csv"
    sample_df.to_csv(output_csv, index=False)

    print(f"Sample evaluation saved to {output_csv}")


def preprocess_pdbs(input_path, output_path):
    """Preprocess input PDB files for ProteinMPNN."""
    pre_output_path = os.path.join(output_path, "parsed_pdbs.jsonl")
    process = subprocess.Popen(
        [
            "python",
            "modules/ProteinMPNN/helper_scripts/parse_multiple_chains.py",
            f"--input_path={input_path}",
            f"--output_path={pre_output_path}",
        ]
    )
    _ = process.wait()
    return pre_output_path


def parse_chains_from_record(record):
    """Parse chains from a FASTA record."""
    chain_seqs = str(record.seq).split("/")
    chains = []
    for j, chain_seq in enumerate(chain_seqs):
        entry = {
            "protein": {
                "id": chr(ord("A") + j),
                "sequence": chain_seq,
                "msa": "empty",
                "cyclic": False,
            }
        }
        chains.append(entry)
    return chains


def extract_score_from_record(record):
    """Extract score from FASTA record description."""
    match = re.search(r"score=([0-9.]+)", record.description)
    return float(match.group(1)) if match else 0.0


def process_fasta_to_yaml(fasta_file, output_yaml_folder, k=4):
    """Process a single FASTA file and create YAML files."""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        print(f"No sequences found in {fasta_file}. Skipping.")
        return []

    # Extract sequences and scores
    homomers = []
    scores = []
    for record in records:
        chains = parse_chains_from_record(record)
        score = extract_score_from_record(record)
        homomers.append(chains)
        scores.append(score)

    # Keep top k sequences
    if len(homomers) > k:
        sorted_indices = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
        top_k_indices = sorted_indices[:k]
        homomers = [homomers[i] for i in top_k_indices]

    # Create YAML files
    yaml_files = []
    for i, homomer in enumerate(homomers):
        yaml_data = {"sequences": homomer}
        yaml_file_path = output_yaml_folder / f"{fasta_file.stem}_{i}.yaml"  # Fixed variable name
        with open(yaml_file_path, "w") as yaml_file:
            yaml.dump(yaml_data, yaml_file, default_flow_style=False)
        yaml_files.append(yaml_file_path)

    return yaml_files


def convert_sequences_to_yaml(mpnn_fasta_path, output_path):
    """Convert ProteinMPNN sequences to Boltz YAML format."""
    fasta_files = [f for f in Path(mpnn_fasta_path).iterdir() if f.suffix in [".fasta", ".fa"]]

    output_yaml_folder = Path(output_path) / "boltz_yaml"
    output_yaml_folder.mkdir(parents=True, exist_ok=True)

    all_yaml_files = []
    for fasta_file in fasta_files:
        yaml_files = process_fasta_to_yaml(fasta_file, output_yaml_folder)
        all_yaml_files.extend(yaml_files)

    return output_yaml_folder, all_yaml_files


def run_protein_mpnn(pre_output_path, output_path, num_seq_per_target=8):
    """Run ProteinMPNN to generate sequences."""
    pmpnn_args = [
        "python",
        "modules/ProteinMPNN/protein_mpnn_run.py",
        "--out_folder",
        f"{output_path}/protein_mpnn",
        "--jsonl_path",
        pre_output_path,
        "--num_seq_per_target",
        str(num_seq_per_target),
        "--sampling_temp",
        "0.1",
        "--seed",
        str(123),
        "--batch_size",
        "1",
    ]

    gpu_id = os.environ.get("CUDA_VISIBLE_DEVICES", "0")
    if gpu_id:
        pmpnn_args.extend(["--device", str(gpu_id)])

    num_tries = 0
    ret = -1
    while ret != 0:  # Fixed: was ret < 0
        try:
            process = subprocess.Popen(pmpnn_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            ret = process.wait()
            if ret == 0:
                break
        except Exception as e:
            num_tries += 1
            torch.cuda.empty_cache()
            if num_tries > 4:
                raise e

    return os.path.join(output_path, "protein_mpnn", "seqs")


def run_boltz_folding(yaml_folder):
    """Run Boltz folding model."""
    process = subprocess.Popen(["boltz", "predict", str(yaml_folder), "--max_parallel_samples", "2", "--out_dir", "eval_symmetric"])
    ret = process.wait()
    if ret != 0:
        raise RuntimeError(f"Boltz folding failed with return code {ret}")

def main(input_path, output_path, num_seq_per_target=8):
    """Main pipeline function."""
    # Create output directory
    Path(output_path).mkdir(parents=True, exist_ok=True)

    # Step 1: Preprocess PDB files for ProteinMPNN
    pre_output_path = preprocess_pdbs(input_path, output_path)

    # Step 2: Run ProteinMPNN
    mpnn_fasta_path = run_protein_mpnn(pre_output_path, output_path, num_seq_per_target)

    # Step 3: Convert FASTA to YAML format for Boltz
    yaml_folder, yaml_files = convert_sequences_to_yaml(mpnn_fasta_path, output_path)

    # Step 4: Run Boltz folding
    run_boltz_folding(yaml_folder)

    # Step 5: Evaluate results
    evaluate_results(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ProteinMPNN and Boltz folding model")
    parser.add_argument("--input_path", type=str, required=True, help="Path to input PDB files")
    parser.add_argument("--output_path", type=str, required=True, help="Path to output directory")
    
    args = parser.parse_args()
    
    main(args.input_path, args.output_path)


