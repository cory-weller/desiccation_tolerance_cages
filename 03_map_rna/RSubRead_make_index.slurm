#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 40G
#SBATCH --time 0-2:00:00
#SBATCH --partition largemem
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv

module load R/3.5.1

Rscript - <<EOF
#!/usr/bin/env Rscript

library(Rsubread)

buildindex(basename="dm3_index", reference="./SNPs/dm3.fa")

EOF
