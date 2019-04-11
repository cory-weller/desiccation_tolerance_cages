# Retrieving data

run `init.sh` to retrieve RNA sequencing reads and other required supplemental files, including:
  * singularity image for pre-compiled iMapSplice
  * a table of SNPs within the DGRP
  * reference FASTA files for dm3 genome, chromosomes 2L, 2R, 3L, 3R and X
  * a gene annotation file for dm3
  * a table of file stems identifying `.fastq` files

## Edit the `slurm.config` file
This file contains variables that will be used when generating required index files and mapping the sequencing reads. Edit any variables that are specific to your intended tasks.

## Build the required index files
While in the main project directory, run `bash ./iMapSplice_map_reads.slurm` to extract variables from the `slurm.config` file that you edited and submit the job to index files. It should take a couple hours, minimum. If you add your email and keep `mail_notifications` to `ALL` then you should receive an email when the job completes.

## Map sequencing reads

Test if everything is configured appropriately by running `iMapSplice_map_reads.slurm` as-is, which submits a job array of size 1. It will only try to map the first sample's reads, i.e. the first line within your file of file stems.

```bash
bash ./iMapSplice_map_reads.slurm
```

If it works appropriately, you can edit the `sbatch` command inside the `iMapSplice_map_reads.slurm`, so that the array spans the number of samples to map. E.g., with 250 samples, you would edit the option to be `--array=1-250%10` to run a total of 250 mapping jobs, with at most 10 running simultaneously.
