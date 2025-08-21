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
from pmpnn import run_pmpnn_processes
from evaluator import Evaluator
import shutil


class ValidationPipeline:
    def __init__(self, folder_path, num_seq_per_target=8):
        self.folder_path = folder_path
        self.num_seq_per_target = num_seq_per_target

        self.evaluator = Evaluator(folder_path)

    def parse_chains_from_record(self, record):
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

    def extract_score_from_record(self, record):
        """Extract score from FASTA record description."""
        match = re.search(r"score=([0-9.]+)", record.description)
        return float(match.group(1)) if match else 0.0

    def process_fasta_to_yaml(self, fasta_file, output_yaml_folder, k=4):
        """Process a single FASTA file and create YAML files."""
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            print(f"No sequences found in {fasta_file}. Skipping.")
            return []

        # Extract sequences and scores
        homomers = []
        scores = []
        for record in records[1:]:  # skip record 0 as it only contains prototype sequence
            chains = self.parse_chains_from_record(record)
            score = self.extract_score_from_record(record)
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
            yaml_file_path = output_yaml_folder / f"{fasta_file.stem}_seq_{i}.yaml"  # Fixed variable name
            with open(yaml_file_path, "w") as yaml_file:
                yaml.dump(yaml_data, yaml_file, default_flow_style=False)
            yaml_files.append(yaml_file_path)

        return yaml_files

    def convert_sequences_to_yaml(self, mpnn_fasta_path, output_path):
        """Convert ProteinMPNN sequences to Boltz YAML format."""
        fasta_files = [f for f in Path(mpnn_fasta_path, "seqs").iterdir() if f.suffix in [".fasta", ".fa"]]

        output_yaml_folder = Path(output_path) / "boltz_files" / "boltz_yaml"
        output_yaml_folder.mkdir(parents=True, exist_ok=True)

        all_yaml_files = []
        for fasta_file in fasta_files:
            yaml_files = self.process_fasta_to_yaml(fasta_file, output_yaml_folder)
            all_yaml_files.extend(yaml_files)

        # Move all fastas to the processed folder
        proc_folder = Path(mpnn_fasta_path, "processed")
        proc_folder.mkdir(parents=True, exist_ok=True)

        self.evaluator.move_seqs_to_processed_folder(fasta_files, proc_folder)

        return output_yaml_folder, all_yaml_files

    def run_protein_mpnn(self, input_path, output_path, num_seq_per_target=8):
        """Run ProteinMPNN to generate sequences."""

        self.evaluator.get_backbones_needing_mpnn(input_path)

        run_pmpnn_processes(
            input_path=input_path + "pdbs/",
            output_dir=output_path,
            symmetry=True,
            seqs=num_seq_per_target,
            sampling_temp=0.2,
            use_soluble_model=False,
        )
        # gpu_id = os.environ.get("CUDA_VISIBLE_DEVICES", "0")
        # if gpu_id:
        #     pmpnn_args.extend(["--device", str(gpu_id)])

        # Save outputs to dataframe
        self.evaluator.evaluate_mpnn_outputs(input_path, output_path)

        return os.path.join(output_path)

    def run_boltz_folding(self, yaml_folder):
        """Run Boltz folding model."""
        self.evaluator.get_yamls_needing_boltz(yaml_folder)

        process = subprocess.Popen(
            [
                "boltz",
                "predict",
                str(yaml_folder),
                "--max_parallel_samples",
                "2",
                "--out_dir",
                f"{self.folder_path}/boltz_files/boltz_out",
            ],
        )
        ret = process.wait()
        if ret != 0:
            raise RuntimeError(f"Boltz folding failed with return code {ret}")

        self.evaluator.evaluate_boltz_outputs(yaml_folder)

    def run_pipeline(self):
        """Main pipeline function."""

        # Create output directory
        Path(self.folder_path).mkdir(parents=True, exist_ok=True)

        # Step 1: Run ProteinMPNN
        mpnn_fasta_path = self.run_protein_mpnn(
            os.path.join(self.folder_path, "backbone_pdbs/"),
            os.path.join(self.folder_path, "pmpnn"),
            self.num_seq_per_target,
        )

        # Eventually build in Foldseek clustering here and give vector with num_folds per backbone to Boltz

        # Step 2: Convert FASTA to YAML format for Boltz
        yaml_folder, yaml_files = self.convert_sequences_to_yaml(mpnn_fasta_path, self.folder_path)

        # Step 3: Run Boltz folding
        self.run_boltz_folding(yaml_folder)

        print("Validation pipeline completed successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ProteinMPNN and Boltz folding model")
    parser.add_argument("--folder_path", type=str, required=True, help="Path to input PDB files")

    args = parser.parse_args()

    pipeline = ValidationPipeline(args.folder_path, num_seq_per_target=8)
    pipeline.run_pipeline()

