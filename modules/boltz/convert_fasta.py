import argparse
import yaml
import pathlib
from pathlib import Path
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Convert FASTA files to a YAML format as required by boltz.")
    parser.add_argument("--input_fasta_folder", type=str, help="Path to the input FASTA folder.")
    parser.add_argument("--output_yaml_folder", type=str, help="Path to the output YAML folder.")

    args = parser.parse_args()

    print(f"Input FASTA folder: {args.input_fasta_folder}")
    print(f"Output YAML folder: {args.output_yaml_folder}")

    if not Path(args.input_fasta_folder).exists():
        raise FileNotFoundError(f"The input folder {args.input_fasta_folder} does not exist.")

    if not Path(args.output_yaml_folder).exists():
        Path(args.output_yaml_folder).mkdir(parents=True, exist_ok=True)
        print(f"Created output folder: {args.output_yaml_folder}")

    fasta_files = [f for f in Path(args.input_fasta_folder).iterdir() if f.suffix in ['.fasta', '.fa']]

    if not fasta_files:
        raise FileNotFoundError("No FASTA files found in the input folder.")

    for fasta_file in fasta_files:
        # load the FASTA file
        records = list(SeqIO.parse(fasta_file, "fasta"))

        # Assume that we only use the monomer (the first record)
        if not records:
            print(f"No sequences found in {fasta_file}. Skipping.")
            continue
        # create yaml
        sequences = []
        entity_type = "protein"
        for record in records:
            entry = {
                entity_type:
                    {'id': record.id,
                    'sequence': str(record.seq),
                    'msa': "empty",
                    'cyclic': False, # This assumes the protein is cyclic on monomer level
                }
            }
            sequences.append(entry)

        # Save the YAML file
        yaml_data = {"sequences": sequences}

        with open(Path(args.output_yaml_folder) / f"{fasta_file.stem}.yaml", 'w') as yaml_file:
            yaml.dump(yaml_data, yaml_file, default_flow_style=False)



if __name__ == "__main__":
    main()