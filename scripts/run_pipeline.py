# Script for running ProteinMPNN and Boltz folding model

import argparse
import sys

from src.pmpnn import run_pmpnn_processes
from src.boltz import run_boltz_processes

from src.evaluator import Evaluator
from src.config import PipelineConfig
from src.file_handling import FileHandler
from src.boltz import FastaToYamlProcessor

from src.config import create_config


class ValidationPipeline:
    def __init__(self, config: PipelineConfig):
        self.config = config

        self.evaluator = Evaluator(config)
        self.file_handler = FileHandler(config)
        self.fasta_processor = FastaToYamlProcessor(config)
    
    def run_protein_mpnn(self):
        """Run ProteinMPNN to generate sequences."""

        self.file_handler.move_processed_backbones(df=self.evaluator.df)

        print("Running ProteinMPNN...")
        run_pmpnn_processes(
            input_path=self.config.backbone_pdbs_dir,
            output_dir=self.config.mpnn_output_dir,
            symmetry=True,
            seqs=self.config.num_seq_per_target,
            sampling_temp=self.config.sampling_temp,
            use_soluble_model=self.config.use_soluble_model,
        )

        self.evaluator.evaluate_mpnn_outputs()
        self.file_handler.move_all_pdbs_to_processed_dir()

    def convert_sequences_to_yaml(self):
        """Convert ProteinMPNN sequences to YAML format for Boltz."""
        self.fasta_processor.convert_all_sequences_to_yaml()

        self.file_handler.move_all_seqs_to_processed_dir()


    def run_boltz_folding(self):
        """Run Boltz folding model."""
        self.file_handler.move_processed_yamls()

        print("Running Boltz folding model...")
        run_boltz_processes(
            yaml_folder=str(self.config.yaml_output_dir),
            max_parallel_samples=self.config.max_parallel_samples,
            out_dir=str(self.config.boltz_output_dir),
        )

        self.evaluator.evaluate_boltz_outputs()
        self.file_handler.move_all_yamls_to_processed_dir()


    def run_pipeline(self):
        """Main pipeline function."""
        # Step 1: Run ProteinMPNN
        self.run_protein_mpnn()

        # Step 2: Convert FASTA to YAML format for Boltz --> eventually build in Foldseek here at some point
        self.convert_sequences_to_yaml() 

        # Step 3: Run Boltz folding
        self.run_boltz_folding()

        print("Validation pipeline completed successfully.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ProteinMPNN and Boltz folding model")
    parser.add_argument("--folder_path", type=str, required=True, help="Path to input PDB files")

    args = parser.parse_args()

    config = create_config(args.folder_path)

    pipeline = ValidationPipeline(config)
    pipeline.run_pipeline()

