#!/bin/bash

#SBATCH --partition=p.psycl 
#SBATCH --cpus-per-task=20 
#SBATCH --mem=200G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

source ~/.bashrc

conda activate sc-downsampling

# Run Rscript from command line
Rscript 00_downsampling_perCell.R

conda deactivate
