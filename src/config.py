"""Configuration management for sym-protein-val pipeline."""

from pathlib import Path
from dataclasses import dataclass


@dataclass
class PipelineConfig:
    """Centralized configuration for all pipeline paths and settings."""
    
    # Base paths
    base_folder: str

    #MPNN settings
    num_seq_per_target: int = 8
    sampling_temp: float = 0.2
    use_soluble_model: bool = True
    k: int = 4

    # Boltz options
    max_parallel_samples: int = 2

    def __post_init__(self):
        """Initialize computed paths after creation."""
        self.base_folder = Path(self.base_folder)

    
    # Input paths
    @property
    def backbones_dir(self) -> Path:
        """Directory containing input PDB files."""
        return self.base_folder / "backbone_pdbs"

    @property
    def backbone_pdbs_dir(self) -> Path:
        """Directory containing input PDB files."""
        return self.backbones_dir / "pdbs"
    
    @property
    def processed_pdbs_dir(self) -> Path:
        """Directory for processed PDB files."""
        return self.backbones_dir / "processed"

    # MPNN paths
    @property
    def mpnn_output_dir(self) -> Path:
        """ProteinMPNN output directory."""
        return self.base_folder / "pmpnn"
    
    @property
    def mpnn_seqs_dir(self) -> Path:
        """ProteinMPNN sequence outputs."""
        return self.mpnn_output_dir / "seqs"

    @property
    def mpnn_processed_dir(self) -> Path:
        """Processed MPNN sequences."""
        return self.mpnn_output_dir / "processed"

    @property
    def parsed_chains(self) -> Path:
        """Parsed chain JSONL file."""
        return self.mpnn_output_dir / "parsed_pdbs.jsonl"

    @property
    def tied_positions(self) -> Path:
        """Tied positions JSONL file."""
        return self.mpnn_output_dir / "tied_pdbs.jsonl"

    # Boltz paths
    @property
    def boltz_base_dir(self) -> Path:
        """Base Boltz output directory."""
        return self.base_folder / "boltz_files"
    
    @property
    def yaml_output_dir(self) -> Path:
        """YAML files for Boltz input."""
        return self.boltz_base_dir / "boltz_yaml"

    @property
    def yaml_processed_dir(self) -> Path:
        """Directory containing YAML files."""
        return self.boltz_base_dir / "processed"
    
    @property
    def boltz_output_dir(self) -> Path:
        """Boltz output directory."""
        return self.boltz_base_dir / "boltz_out"

    @property
    def boltz_results_dir(self) -> Path:
        """Boltz results directory."""
        return self.boltz_output_dir / "boltz_results_boltz_yaml"

    @property
    def boltz_predictions_dir(self) -> Path:
        """Boltz predictions directory."""
        return self.boltz_results_dir / "predictions"

    @property
    def boltz_confidence_pattern(self) -> str:
        """Pattern for confidence JSON files."""
        return "**/confidence*.json"
    
    # Data management paths
    @property
    def dataframe_path(self) -> Path:
        """Main dataframe parquet file."""
        return self.base_folder / "dataframe.parquet"
    

    def create_directories(self) -> None:
        """Create all necessary directories."""
        directories = [
            self.backbone_pdbs_dir,
            self.processed_pdbs_dir,
            self.mpnn_output_dir,
            self.mpnn_processed_dir,
            self.mpnn_seqs_dir,
            self.yaml_output_dir,
            self.yaml_processed_dir,
            self.boltz_base_dir,
            self.boltz_output_dir,
            self.boltz_results_dir,
            self.boltz_predictions_dir
        ]

        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)


def create_config(base_folder: str, **kwargs) -> PipelineConfig:
    """Factory function to create configuration with validation."""
    config = PipelineConfig(
        base_folder=base_folder,
        **kwargs
    )
    config.create_directories()
    return config
