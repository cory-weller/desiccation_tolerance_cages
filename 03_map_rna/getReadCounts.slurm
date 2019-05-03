#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 5
#SBATCH --mem 40G
#SBATCH --time 0-01:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv

#module load samtools
module load gatk

mapper=${1}

function getReadCounts {
  mapper=${1}
  filename=${2}
  echo "Tabulating allele-specific reads with GATK ASEReadCounter"
  gatk ASEReadCounter \
  --input mapped_reads/${mapper}/${filename} \
  --output ../04_expression_analysis/ASE/${mapper}/${filename%.rg.sort.bam}.readcounts \
  --variant SNPs/dm3.variants \
  --reference SNPs/dm3.fa && \
  gzip ../04_expression_analysis/ASE/${mapper}/${filename%.rg.sort.bam}.readcounts
}

export -f getReadCounts

parallel -j 5 getReadCounts ::: ${mapper} ::: $(ls ./mapped_reads/${mapper}/ | grep -v "bai")

echo "done"
