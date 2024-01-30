#!/bin/bash

#SBATCH --partition=p.psycl 
#SBATCH --cpus-per-task=20 
#SBATCH --mem=200G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

source ~/.bashrc

conda activate sc-norm-sct

# Run Python script
# -u option required to save print statements in out file
python -u 02_sct_norm.py

conda deactivate
