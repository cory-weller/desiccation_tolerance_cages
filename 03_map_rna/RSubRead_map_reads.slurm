#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 5
#SBATCH --mem 40G
#SBATCH --time 0-12:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv

cd ${projectDirectory}

Rscript - <<EOF


EOF
