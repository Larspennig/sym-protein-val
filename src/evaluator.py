import pandas as pd
from pathlib import Path
import shutil
from Bio import SeqIO
from .config import PipelineConfig


class Evaluator:
    """Evaluator class to manage and process MPNN and Boltz outputs."""
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.folder_path = config.base_folder  # Keep for backward compatibility

        self.load_or_initialize_dataframe()

    def load_or_initialize_dataframe(self):
        """Load existing dataframe or create new one."""
        if self.config.dataframe_path.exists():
            self.df = pd.read_parquet(self.config.dataframe_path)
        else:
            self.df = pd.DataFrame(
                columns=[
                    "seq",
                    "seq_name",
                    "seq_file",
                    "backbone_file",
                    "boltz_output_file",
                    "confidence_score",
                    "complex_plddt",
                    "complex_iplddt",
                    "complex_ipde",
                    "pmpnn_description",
                ]
            )
            # save empty DataFrame to parquet
            self.df.to_parquet(self.config.dataframe_path, index=False)

    def evaluate_mpnn_outputs(self):
        """Evaluate MPNN outputs and save to dataframe."""
        out_seqs = [f for f in self.config.mpnn_seqs_dir.iterdir() if f.suffix in [".fasta", ".fa"]]

        backbones = []
        new_rows = []

        for seq_file in out_seqs:
            seqs = list(SeqIO.parse(str(seq_file), "fasta"))
            for i, seq in enumerate(seqs[1:]):  # Skip the first sequence as it only contains the dummy sequence
                seq_str = str(seq.seq)
                seq_name = seq_file.stem + f"_seq_{i}"
                # Check if sequence already exists in dataframe
                if seq_str in self.df['seq'].values:
                    print(f"Sequence {seq_name} already exists in dataframe, skipping.")
                    continue

                # Create a new row in the dataframe
                new_row = {
                    "seq": seq_str,
                    "seq_file": seq_file.name,
                    "seq_name": seq_name,
                    "backbone_file": seq_file.stem + ".pdb",
                    "pmpnn_description": seq.description, # @TODO parse the description
                }
                backbones.append(seq_file.stem + ".pdb")
                new_rows.append(new_row)
        
        # Add new rows to the dataframe
        if new_rows:
            self.df = pd.concat([self.df, pd.DataFrame(new_rows)], ignore_index=True)

        # Add fasta file to the dataframe
        self.df.to_parquet(self.config.dataframe_path, index=False)

    def evaluate_boltz_outputs(self):
        """Evaluate Boltz outputs and update dataframe."""
        sample_paths = list(self.config.boltz_predictions_dir.glob(self.config.boltz_confidence_pattern))
        
        for sample in sample_paths:
            seq_name = sample.stem.split("confidence_")[1].split("_model")[0]

            # Check if sequence exists in dataframe
            mask = self.df["seq_name"] == seq_name

            if not mask.any():
                # Create new row, load JSON and store all values as lists
                print(f"Boltz sample {seq_name} not found in dataframe, creating new row.")
                try:
                    sample_df = pd.read_json(sample)
                except ValueError as e:
                    print(f"Error reading {seq_name}: {e}")
                    continue

                sample_row = {
                    "seq_name": seq_name,
                    "boltz_output_file": str(sample),
                }
                # Store all values as lists
                for key in sample_df.keys():
                    sample_row[key] = sample_df[key].tolist()

                self.df = pd.concat([self.df, pd.DataFrame([sample_row])], ignore_index=True)

            else:
                # Update existing row
                row_index = mask.idxmax()
                has_boltz_results = pd.notna(self.df.loc[row_index, "boltz_output_file"])

                if has_boltz_results:
                    print(f"Sample {seq_name} already has Boltz results, skipping.")
                    continue
                else:
                    print(f"Sample {seq_name} exists in dataframe, adding boltz values.")
                    try:
                        sample_df = pd.read_json(sample)
                    except ValueError as e:
                        print(f"Error reading {seq_name}: {e}")
                        continue

                    # Store all values as lists
                    for key in sample_df.keys():
                        if key not in self.df.columns:
                            self.df[key] = None
                        self.df.at[row_index, key] = sample_df[key].tolist()

                    self.df.loc[row_index, "boltz_output_file"] = str(sample)

        # update df parquet
        self.df.to_parquet(self.config.dataframe_path, index=False)
