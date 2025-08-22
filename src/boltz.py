# Boltz processing script
import re
import subprocess
from pathlib import Path
from typing import List, Dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from .config import PipelineConfig
import yaml


class FastaToYamlProcessor:
    """Handles FASTA to YAML conversion logic."""
    def __init__(self, config: PipelineConfig):
        self.config = config

    def parse_chains_from_record(self, record: SeqRecord) -> List[Dict]:
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

    def extract_score_from_record(self, record: SeqRecord) -> float:
        """Extract score from FASTA record description."""
        match = re.search(r"score=([0-9.]+)", record.description)
        return float(match.group(1)) if match else 0.0

    def process_fasta_to_yaml(self, fasta_file: Path) -> List[Path]:
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
        if len(homomers) > self.config.k:
            sorted_indices = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
            top_k_indices = sorted_indices[: self.config.k]
            homomers = [homomers[i] for i in top_k_indices]

        # Create YAML files
        yaml_files = []
        for i, homomer in enumerate(homomers):
            yaml_data = {"sequences": homomer}
            yaml_file_path = self.config.yaml_output_dir / f"{fasta_file.stem}_seq_{i}.yaml"
            with open(yaml_file_path, "w") as yaml_file:
                yaml.dump(yaml_data, yaml_file, default_flow_style=False)
            yaml_files.append(yaml_file_path)

        return yaml_files

    def convert_all_sequences_to_yaml(self) -> List[Path]:
        """Convert all ProteinMPNN sequences to Boltz YAML format."""
        fasta_files = [f for f in self.config.mpnn_seqs_dir.iterdir() if f.suffix in [".fasta", ".fa"]]

        all_yaml_files = []
        for fasta_file in fasta_files:
            yaml_files = self.process_fasta_to_yaml(fasta_file)
            all_yaml_files.extend(yaml_files)

        return all_yaml_files


def run_boltz_processes(yaml_folder: str, max_parallel_samples: int = 2, out_dir: str = None):
    """Run the Boltz folding processes."""

    process = subprocess.Popen(
            [
                "boltz",
                "predict",
                str(yaml_folder),
                "--max_parallel_samples",
                str(max_parallel_samples),
                "--out_dir",
                str(out_dir),
            ],
        )
    ret = process.wait()
    if ret != 0:
        raise RuntimeError(f"Boltz folding failed with return code {ret}")

