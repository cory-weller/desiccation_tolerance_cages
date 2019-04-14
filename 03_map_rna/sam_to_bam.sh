#!/usr/bin/env bash

# To convert unconverted .sam files to .bam
module load samtools/1.9

cd mapped_reads/

while read filename; do
  samtools view -b ${filename} | samtools sort - > ${filename%.sam}.bam && samtools index ${filename%.sam}.bam
done < <(ls *.sam)
