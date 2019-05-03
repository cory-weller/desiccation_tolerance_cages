#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 8G
#SBATCH --time 0-24:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab

#add_read_groups.slurm

# type should be iMapSplice or RSubRead
type=${1}

module load samtools/1.9
module load gatk

cd mapped_reads/${type}

while read filename; do
  filestem=${filename%.[sb]am}
  # If .sam (instead of .bam), convert to .bam
  if [ ${filename: -4} == ".sam" ]; then
    samtools view -b ${filename} > ${filename%.sam}.bam && rm ${filename}
  fi

  samtools sort ${filestem}.bam > ${filestem}.sort.bam
  samtools index ${filestem}.sort.bam

  gatk AddOrReplaceReadGroups \
  --I=${filestem}.sort.bam \
  --O=${filestem}.rg.bam \
  --RGID=${filestem} \
  --RGLB=indoorCages \
  --RGPL=illumina \
  --RGPU=indoorCages \
  --RGSM=${filestem}

  samtools sort ${filestem}.rg.bam > ${filestem}.rg.sort.bam && \
  samtools index ${filestem}.rg.sort.bam && \
  rm ${filestem}.bam ${filestem}.rg.bam ${filestem}.sort.bam ${filestem}.sort.bam.bai

done < <(ls *.[sb]am)
