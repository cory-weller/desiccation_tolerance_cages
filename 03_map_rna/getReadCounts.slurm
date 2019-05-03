#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 8G
#SBATCH --time 0-01:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv

module load samtools
module load gatk

mapper=${1}

filestem=$(sed -n ${SLURM_ARRAY_TASK_ID}p filestems.txt)

echo "Tabulating allele-specific reads with GATK ASEReadCounter"
gatk ASEReadCounter \
--input mapped_reads/${mapper}/${filestem}.rg.sort.bam \
--output ../04_expression_analysis/ASE/${mapper}/${filestem}.readcounts \
--variant SNPs/dm3.variants \
--reference SNPs/dm3.fa && \
gzip ../04_expression_analysis/ASE/${mapper}/${filestem}.readcounts

echo "done"
