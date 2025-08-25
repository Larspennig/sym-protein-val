#!/bin/bash

#SBATCH -N 1 
#SBATCH --job-name=SymVal
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu_a100
#SBATCH -t 2-00:00 
#SBATCH -o "./slurm_outputs/slurm-%j.out"
#SBATCH --mem=32000 
#SBATCH --qos=default
#SBATCH --cpus-per-task=4
 
cd ${SLURM_SUBMIT_DIR}
echo Starting job ${SLURM_JOBID}w
echo SLURM assigned me these nodes:
squeue -j ${SLURM_JOBID} -O nodelist | tail -n +2

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate sym_val
export PYTHONPATH=${SLURM_SUBMIT_DIR}:$PYTHONPATH

python scripts/run_pipeline.py --folder_path ff_output/C_test/