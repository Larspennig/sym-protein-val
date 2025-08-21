import pandas as pd
from pathlib import Path
import os
import shutil
from Bio import SeqIO


class Evaluator:
    """Evaluator class to manage and process MPNN and Boltz outputs."""
    def __init__(self, folder_path):
        self.folder_path = folder_path

        self.load_or_initialize_dataframe()

    def load_or_initialize_dataframe(self):
        df_path = Path(self.folder_path) / "dataframe.parquet"
        if df_path.exists():
            self.df = pd.read_parquet(df_path)
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
                ]
            )
            # save empty DataFrame to parquet
            self.df.to_parquet(df_path, index=False)

    def get_backbones_needing_mpnn(self, backbone_pdb_path=None):
        """Check which backbones need MPNN processing."""

        # create processed folder if it doesn't exist
        new_folder = Path(backbone_pdb_path) / "pdbs"
        proc_folder = Path(backbone_pdb_path) / "processed"
        proc_folder.mkdir(parents=True, exist_ok=True)

        if not any(new_folder.iterdir()):
            print("No new backbones provided.")
            return

        for file in new_folder.glob("*.pdb"):
            # Check if the file is already processed
            if file.name in self.df['backbone_file'].values:
                print(f"File {file.name} already processed, skipping.")
                if file in proc_folder.iterdir():
                    file.unlink()  # Remove file since it's already processed
                    print(f"Removed {file.name} from new folder as it was already processed.")
                else:
                    shutil.move(file, proc_folder / file.name)
                    print(f"Moved {file.name} to processed folder.")

        # # Move file out of the folder if it was already processed
        # for fasta_file in Path(mpnn_fasta_path).glob("*.fasta"):
        #     # Move processed fasta files to the processed folder
        #     destination = proc_folder / fasta_file.name
        #     if fasta_file in self.df['seq_file'].values:
        #         if destination.exists():
        #             print(f"File {fasta_file.name} already exists in processed folder, skipping move.")
        #             fasta_file.unlink()  # Remove original since it's already processed
        #         else:
        #             shutil.move(fasta_file, destination)
        #             print(f"Moved {fasta_file.name} to processed folder.")


    def get_yamls_needing_boltz(self, yaml_folder):
        # Check which yamls have been processed

        pred_folder = Path(self.folder_path) / "boltz_files" / "boltz_out" / "boltz_results_boltz_yaml" / "predictions"
        pred_folder.mkdir(parents=True, exist_ok=True)

        # list all subfolder names
        subfolders = [f.name for f in pred_folder.iterdir() if f.is_dir()]
        processed_yaml_files = [f for f in Path(yaml_folder).iterdir() if f.stem in subfolders]

        # Move processed yaml files to the processed folder
        for yaml_file in processed_yaml_files:
            destination = Path(self.folder_path) / "boltz_files" / "processed" / yaml_file.name
            if destination.exists():
                print(f"File {yaml_file.name} already exists in processed folder, skipping move.")
                yaml_file.unlink()
            else:
                if yaml_file.exists():
                    shutil.move(yaml_file, destination)
                    print(f"Moved {yaml_file.name} to processed folder.")


    def move_seqs_to_processed_folder(self, fasta_files, proc_folder):
        for file in fasta_files:
            if file.name in self.df['seq_file'].values:
                print(f"File {file.name} already processed.")
                if proc_folder / file.name in proc_folder.iterdir():
                    print(f"File {file.name} already exists in processed folder, skipping move.")
                    file.unlink()
                else:
                    shutil.move(file, proc_folder / file.name)
                    print(f"Moved {file.name} to processed folder.")


    def evaluate_mpnn_outputs(self, backbone_path, mpnn_path):
        """Evaluate MPNN outputs and save to dataframe."""
        backbone_path = Path(backbone_path)
        
        seqs_path = Path(mpnn_path) / "seqs"
        # proc_folder = Path(mpnn_path) / "processed"
        # if not seqs_path.exists():
        #     print(f"No sequences found in {seqs_path}.")
        #     return
        # if not proc_folder.exists():
        #     proc_folder.mkdir(parents=True, exist_ok=True)

        out_seqs = [f for f in seqs_path.iterdir() if f.suffix in [".fasta", ".fa"]]

        backbones = []

        for seq_file in out_seqs:
            seqs = list(SeqIO.parse(str(seq_file), "fasta"))
            for i, seq in enumerate(seqs):
                seq_str = str(seq.seq)
                # Check if sequence already exists in dataframe
                if seq_str in self.df['seq'].values:
                    print(f"Sequence {seq_str} already exists in dataframe, skipping.")
                    continue

                # Create a new row in the dataframe
                new_row = {
                    "seq": seq_str,
                    "seq_file": seq_file.name,
                    "seq_name": seq_file.stem + f"_seq_{i}",
                    "backbone_file": seq_file.stem + ".pdb",
                    "pmpnn_description": seq.description, # @TODO parse the description
                }
                backbones.append(seq_file.stem + ".pdb")
                self.df = pd.concat([self.df, pd.DataFrame([new_row])], ignore_index=True)

        # Move processed pdbs to the processed folder
        for pdb_file in backbones:
            pdb_path = backbone_path / "pdbs" / pdb_file
            destination = backbone_path / "processed" / pdb_file
            if destination.exists() and pdb_path.exists():
                pdb_path.unlink()  # Remove file since it's already processed
                print(f"Removed {pdb_file} from new folder as it was already processed.")
            else:
                if pdb_path.exists():
                    shutil.move(backbone_path / "pdbs" / pdb_file, backbone_path / "processed" / pdb_file)

        # Add fasta file to the dataframe
        self.df.to_parquet(os.path.join(self.folder_path, "dataframe.parquet"))

    def evaluate_boltz_outputs(self, boltz_folder):
        sample_paths = list(Path(self.folder_path,"boltz_files","boltz_out","boltz_results_boltz_yaml","predictions" ).glob("**/confidence*.json"))
        
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
        self.df.to_parquet(os.path.join(self.folder_path, "dataframe.parquet"))

        print("Sample evaluation saved. Folding done.")

