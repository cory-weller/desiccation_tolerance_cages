#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 24G
#SBATCH --time 0-4:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv

module load singularity

# extract this run's filestem as nth line from RNAreadFileStemList
# where n is the nth run within the SLURM job array

filestem=$(sed -n ${SLURM_ARRAY_TASK_ID}p filestems.txt)
tmpdir=/dev/shm/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}
working_directory=$(pwd)
mkdir -p ${tmpdir}

echo "working on sample ${filestem}"

gunzip -c ./rna_reads/${filestem}_R1_001.fastq.gz > ${tmpdir}/${filestem}_1.fastq || \
{ echo "could not open ./rna_reads/${filestem}_R1_001.fastq.gz"; exit 1; }

gunzip -c ./rna_reads/${filestem}_R2_001.fastq.gz > ${tmpdir}/${filestem}_2.fastq || \
{ echo "could not open ./rna_reads/${filestem}_R2_001.fastq.gz"; exit 1; }


singularity exec iMapSplice.simg /iMapSplice/code/bin/iMapSplice-unphased \
-P ./SNPs/dgrp2.snps \
-Q ./SNPmerIndex/ \
-G ./globalIndex/ \
-L ./localIndex/ \
-1 ${tmpdir}/${filestem}_1.fastq \
-2 ${tmpdir}/${filestem}_2.fastq \
-T 4 \
-O ${tmpdir} || { echo "iMapSplice mapping failed"; rm -rf  ${tmpdir}; exit 1; }

echo "adding read groups and sorting"


function addReadGroups {
  filename=${1}
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
}

export -f addReadGroups

module load samtools
module load gatk

cd $tmpdir && rm -rf *.fastq phase1_output phase2_output

mv output.sam ${filestem}.sam

addReadGroups ${filestem}.sam

mv ${filestem}.rg.sort.bam* ${working_directory}/mapped_reads/iMapSplice/ && \
cd ${working_directory} && rm -rf ${tmpdir} && \
echo "done"
