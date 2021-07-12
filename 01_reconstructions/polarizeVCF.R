#!/usr/bin/env Rscript

# argument 1 = infilename (as .vcf.gz)
# argument 2 = outfilename (will be output as uncompressed .vcf)

# VCF file should already be curated to only include bi-allelic SNPs (no indels, duplications, CNV)
# and individuals should have homozygous genotypes (fully inbred)
# any heterozygous sites will be set to NA

library(data.table)

importVCF <- function(filename) {
	DT <- fread(cmd="zcat ", filename, sep=" "), na.strings=c(".","./.", "./1", "./0", "0/.", "1/.", "0/1", "1/0"), showProgress=FALSE, header=TRUE, skip="#CHROM")
    polarizeVCF(DT)
    return(DT)
}

polarizeVCF <- function(DT) {
    # Polarizes a VCF data.table to -1 for alt and +1 for ref allele
    # save column names of lineIDs for subset with .SD
    cols <- colnames(DT)[10:ncol(DT)]

	    # convert any NA values to 0
    for (j in cols) {
        set(DT, which(DT[[j]]=="0/0"), j, "0")
    }
    for (j in cols) {
        set(DT, which(DT[[j]]=="1/1"), j, "1")
    }

    # convert anything else to 0
    for (j in cols) {
        set(DT, which(! DT[[j]] %in% c("0","1")), j, ".")
    }

    invisible(DT[])
}


# Parse arguments
	args <- commandArgs(trailingOnly=TRUE)
	infilename <- as.character(args[[1]])
	outfilename <- as.character(args[[2]])

# Load VCF file containing founder Haplotypes, and polarize VCF such that reference = 0; alt = -1; missing data = 0
	vcf.DT <- importVCF(infilename)


# Write table
	write.table(vcf.DT, outfilename, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
