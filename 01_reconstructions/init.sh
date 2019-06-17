#!/usr/bin/env bash

# Retrieve and unzip harp, htslib, GATK, etc
wget -L https://virginia.box.com/shared/static/ntqgqmkl5xtw10a3awvyoazptkbhyfnx.gz -O etc.tar.gz
tar -zxvf etc.tar.gz && rm etc.tar.gz

# Retrieve and untar low-coverage DNA reads
wget -L https://virginia.box.com/shared/static/hne9mcpr3nd3yur11aeg9zwfc7m5cytz.tar -O low_coverage_dna.tar && \
tar -xf low_coverage_dna.tar -d ./dna_reads/ && \
rm low_coverage_dna.tar

# Retrieve reference genome from FlyBase
wget ftp://ftp.flybase.net/releases/FB2014_03/dmel_r5.57/fasta/dmel-all-chromosome-r5.57.fasta.gz -O dm3.fasta.gz && \
gunzip dm3.fasta.gz

# Retrieve dgrp2.vcf

# polarize VCF

Rscript polarizeVCF.Rscript dgrp2.vcf haplotypes.polarized && gzip haplotypes.polarized
