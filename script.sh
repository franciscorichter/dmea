#!/bin/bash
#SBATCH --job-name=R_job
#SBATCH --time=00:11:00
#SBATCH --ntasks=1
#SBATCH --mem=1000
#SBATCH --partition=short
pwd
module load R/3.3.1-foss-2016a
module list
Rscript test1.R
