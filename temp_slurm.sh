#!/bin/bash
#SBATCH --job-name=my_r_job
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=01:00:00
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

module load R   # or your cluster’s way of loading R

Rscript my_script.R