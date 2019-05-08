#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 12G
#SBATCH --time 0-01:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab


module load gcc
module load R/3.5.1

Rscript - <<EOF
#!/usr/bin/env Rscript

library(GWASTools)
library(data.table)
library(SNPRelate)
library(gdsfmt)

vcf_filename <- "filtered_estimates/filtered.all.vcf"
gds_filename <- "filtered_estimates/filtered.all.vcf.gds"

snpgdsVCF2GDS(vcf_filename, gds_filename, method=c("biallelic.only"), snpfirstdim = FALSE)

EOF
