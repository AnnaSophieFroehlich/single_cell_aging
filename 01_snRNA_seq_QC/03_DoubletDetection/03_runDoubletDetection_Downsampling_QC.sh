#!/bin/bash

#SBATCH --partition=p.psycl
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --time=11-00:00:00

source ~/.bashrc

conda activate sc-mar2021

# Run Jupyter notebook from command line
jupyter nbconvert --to notebook --execute 02_DoubletDetection_Downsampling_QC.ipynb --inplace --ExecutePreprocessor.timeout=-1

conda deactivate
