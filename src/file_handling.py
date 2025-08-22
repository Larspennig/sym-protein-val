import shutil
from .config import PipelineConfig
from pathlib import Path

class FileHandler:
    def __init__(self, config: PipelineConfig):
        self.config = config

        # Extract commonly used paths at initialization
        self.backbone_pdbs_dir = config.backbone_pdbs_dir
        self.processed_pdbs_dir = config.processed_pdbs_dir
        self.yaml_output_dir = config.yaml_output_dir
        self.yaml_processed_dir = config.yaml_processed_dir
        self.boltz_predictions_dir = config.boltz_predictions_dir
        self.mpnn_seqs_dir = config.mpnn_seqs_dir
        self.mpnn_processed_dir = config.mpnn_processed_dir

    def _move_or_remove_file(self, source_path: Path, dest_path: Path) -> None:
        """Move file to destination or remove if destination exists."""
        if dest_path.exists():
            source_path.unlink()
            print(f"Removed {source_path.name} (already exists in destination).")
        else:
            shutil.move(source_path, dest_path)
            print(f"Moved {source_path.name} to {dest_path.parent.name}/ folder.")

    def move_all_yamls_to_processed_dir(self):
        """Moves all YAML files to the processed folder."""
        yaml_files = list(self.yaml_output_dir.glob("*.yaml"))
        if not yaml_files:
            print("No new YAML files provided.")
            return
        for yaml_file in yaml_files:
            destination = self.yaml_processed_dir / yaml_file.name
            self._move_or_remove_file(yaml_file, destination)

    def move_all_pdbs_to_processed_dir(self):
        """Moves all pdb files to the processed folder."""

        backbone_files = list(self.backbone_pdbs_dir.glob("*.pdb"))

        # Move processed pdbs to the processed folder
        for pdb_file in backbone_files:
            destination = self.processed_pdbs_dir / pdb_file.name
            self._move_or_remove_file(pdb_file, destination)

    def move_all_seqs_to_processed_dir(self):
        """Move all sequences to the processed folder."""
        fasta_files = list(self.mpnn_seqs_dir.glob("*.fa"))
        proc_folder = self.mpnn_processed_dir

        for fasta in fasta_files:
            destination = proc_folder / fasta.name
            self._move_or_remove_file(fasta, destination)

    def move_processed_backbones(self, df):
        """Check which backbones need MPNN processing."""

        backbone_files = list(self.backbone_pdbs_dir.glob("*.pdb"))

        if not backbone_files:
            print("No new backbones provided.")
            return

        for pdb in backbone_files:
            # Check if the file is already processed
            if pdb.name in df["backbone_file"].values:
                print(f"File {pdb.name} already processed, skipping.")
                destination = self.processed_pdbs_dir / pdb.name
                self._move_or_remove_file(pdb, destination)

    def move_processed_yamls(self):
        # Check which yamls have been processed
        yaml_files = list(self.yaml_output_dir.glob("*.yaml"))

        if not yaml_files:
            print("No new YAML files provided.")
            return
        
        # list all subfolder names and check if they match the yaml files
        subfolders = [f.name for f in self.boltz_predictions_dir.iterdir() if f.is_dir()]
        processed_yaml_files = [f for f in yaml_files if f.stem in subfolders]

        # Move processed yaml files to the processed folder
        for yaml_file in processed_yaml_files:
            destination = self.yaml_processed_dir / yaml_file.name
            self._move_or_remove_file(yaml_file, destination)
