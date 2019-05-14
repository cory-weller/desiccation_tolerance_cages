#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 8G
#SBATCH --time 0-05:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab



module load gcc
module load R/3.5.1

Rscript ./getFeatureCounts.Rscript
