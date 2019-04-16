# Retrieving data

run `init.sh` to retrieve RNA sequencing reads and other required supplemental files, including:
  * singularity image for pre-compiled iMapSplice
  * a table of SNPs within the DGRP
  * reference FASTA files for dm3 genome, chromosomes 2L, 2R, 3L, 3R and X
  * a gene annotation file for dm3
  * a table of file stems identifying `.fastq` files

## Edit the `slurm.config` file
This file contains variables that will be used when generating required index files and mapping the sequencing reads. Edit any variables that are specific to your intended tasks.

## Build index files
While in the main project directory, run `bash ./iMapSplice_make_index.slurm` to extract variables from the `slurm.config` file that you edited and submit the job to index files. It should take a couple hours. If you add your email and keep `mail_notifications` to `ALL` then you should receive an email when the job completes.

## Build SNPmers
Run `bash ./iMapSplice_make_SNPmers.slurm` to generate the second requirement prior to mapping reads.

## Map sequencing reads
Map reads with `bash iMapSplice_map_reads.slurm`

## Convert from `.sam` to `.bam` format
Run `bash sam_to_bam.slurm` to convert to `.bam` format.

## Get read counts
Run `bash getReadCounts.slurm` to add read groups, sort and index `.bam` files, and run the `ASEReadCounter` tool from `GATK` on each sample. Output is placed in `../04_expression_analysis/ASE`.
