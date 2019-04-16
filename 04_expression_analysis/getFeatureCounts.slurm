source ../03_map_rna/slurm.config

sbatch << EOF
#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 8G
#SBATCH --time 0-05:00:00
#SBATCH --partition standard
#SBATCH --account ${allocation}
#SBATCH --mail-type ${mail_notifications}
#SBATCH --mail-user ${userEmail}



module load gcc
module load R/3.5.1

Rscript ./getFeatureCounts.Rscript
