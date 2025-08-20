import argparse
import os
from pdb import pm
import subprocess


def run_pmpnn_processes(input_path, output_dir, symmetry=True, seqs=10, sampling_temp=0.2, use_soluble_model=False):
    """
    Run subprocesses to process PDB files and run a prediction model.

    Args:
    input_path (str): Path to the input directory containing PDB files.
    output_dir (str): Path where the outputs will be saved.
    """
    # Ensure the output directory exists
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Set up the path for the output JSONL file
    path_for_parsed_chains = os.path.join(output_dir, "parsed_pdbs.jsonl")
    path_for_tied_positions = os.path.join(output_dir, "tied_pdbs.jsonl")

    # Run the parsing script
    subprocess.run(
        [
            "python",
            "./modules/ProteinMPNN/helper_scripts/parse_multiple_chains.py",
            "--input_path",
            input_path,
            "--output_path",
            path_for_parsed_chains,
        ]
    )

    # Set parameters for proteinmpnn
    if symmetry:
        # Set up tied positions if symmetry is True
        subprocess.run(
            [
                "python",
                "./modules/ProteinMPNN/helper_scripts/make_tied_positions_dict.py",
                "--input_path",
                path_for_parsed_chains,
                "--output_path",
                path_for_tied_positions,
                "--homooligomer",
                "1",
            ]
        )

        # Set pmpnn args with symmetry
        pmpnn_args = [
            "python",
            "./modules/ProteinMPNN/protein_mpnn_run.py",
            "--jsonl_path",
            path_for_parsed_chains,
            "--tied_positions_jsonl",
            path_for_tied_positions,
            "--out_folder",
            output_dir,
            "--num_seq_per_target",
            f"{seqs}",
            "--sampling_temp",
            f"{sampling_temp}",
            "--seed",
            "37",
            "--batch_size",
            f"{min(seqs, 10)}",
        ]
    else:
        # Set pmpnn args
        pmpnn_args = [
            "python",
            "./modules/ProteinMPNN/protein_mpnn_run.py",
            "--jsonl_path",
            path_for_parsed_chains,
            "--out_folder",
            output_dir,
            "--num_seq_per_target",
            f"{seqs}",
            "--sampling_temp",
            "0.2",
            "--seed",
            "37",
            "--batch_size",
            f"{min(seqs,10)}",
        ]
    if use_soluble_model:
        pmpnn_args.append("--use_soluble_model")

    gpu_id = os.environ.get("CUDA_VISIBLE_DEVICES", "0")
    if gpu_id:
        pmpnn_args.extend(["--device", str(gpu_id)])

    # Run the ProteinMPNN prediction script
    subprocess.run(pmpnn_args)


def main(input_path, output_path, symmetry=True, seqs=10, sampling_temp=0.2, use_soluble_model=False):
    """
    Main function to process PDB files and predict using a machine learning model.

    Args:
    input_path (str): Directory containing the input PDB files.
    output_path (str): Directory where the outputs will be saved.
    """
    run_pmpnn_processes(
        input_path,
        output_path,
        symmetry=symmetry,
        seqs=seqs,
        sampling_temp=sampling_temp,
        use_soluble_model=use_soluble_model,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB files and predict outcomes.")
    parser.add_argument("input_path", type=str, help="Directory containing the input PDB files.")
    parser.add_argument("output_path", type=str, help="Directory where the outputs will be saved.")
    parser.add_argument("symmetry", type=bool, help="Whether pdb is symmetric protein or not", default=True)
    parser.add_argument("seqs", type=int, help="Number of generated seqs per pdb.", default=10)
    parser.add_argument("sampling_temp", type=float, help="Sampling temperature for the prediction model.", default=0.2)
    parser.add_argument(
        "--use_soluble_model",
        action="store_true",
        default=False,
        help="Flag to load ProteinMPNN weights trained on soluble proteins.",
    )

    args = parser.parse_args()

    # Call main function using arguments from the command line
    main(
        args.input_path,
        args.output_path,
        symmetry=args.symmetry,
        seqs=args.seqs,
        sampling_temp=args.sampling_temp,
        use_soluble_model=args.use_soluble_model,
    )
