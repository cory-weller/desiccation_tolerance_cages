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

library(data.table)
library(foreach)

#some functions

reconstruction_dir = "../estimate_genomes/"

filestems <- readLines("filestems.txt")
chromosomes <- c("2L", "2R", "3L", "3R", "X")


#imports haplotype path files and pulls individual id out of file name
import_vcf <- function(reconstruction_dir, filestem, chromosome){
    command <- paste(reconstruction_dir, "/", filestem, ".", chromosome, ".estimate.vcf", sep="")
    fread(command, header=FALSE, col.names=c("CHROM","POS", "GENO"))
}


import_path <- function(reconstruction_dir, filestem_extended, chromosome){
    command <- paste("unzip -p ", reconstruction_dir, "/", filestem_extended, ".zip ", filestem_extended, ".", chromosome, ".estimate.haps", sep="")
    fread(cmd=command, header=TRUE)
}

import_filter <- function(reconstruction_dir, filestem){
    command <- paste(reconstruction_dir, "/", filestem, ".paths.filter", sep="")
    fread(command, header=TRUE)
}

filterInner <- function(vcf, vcf_filter) {
    setkey(vcf_filter, chromosome, start, stop)

    setnames(vcf, "CHROM", "chromosome")
    vcf[, "stop" := POS]
    setnames(vcf, "POS", "start")

    vcf[, id := 1:.N]

    ids_within_vcf_filter_ranges <- unique(foverlaps(vcf, vcf_filter[par=="unknown"], type="within")[!is.na(start), id])

    vcf <- vcf[! id %in% ids_within_vcf_filter_ranges]

    setnames(vcf, "start", "POS")
    setnames(vcf, "chromosome", "CHROM")

    vcf[, c("stop","id") := NULL]
    vcf[, GENO := paste(A1, A2, sep="/")]
    vcf[, c("A1","A2") := NULL]
    return(vcf)
}



filterOuter <- function(filestem, reconstruction_dir, chromosomes) {
    vcf_filter <- import_filter(reconstruction_dir, filestem)
    vcf <- foreach(chromosome.i=chromosomes, .combine="rbind", .errorhandling="remove") %do% {
        tmp <- import_vcf(reconstruction_dir, filestem, chromosome.i)
        tmp[, c("A1","A2") := tstrsplit(GENO, split="/")]
        return(tmp[,c("CHROM","POS","A1","A2")])
    }
    filterInner(vcf, vcf_filter)
}

for(filestem in filestems) {
    vcf.filtered <- filterOuter(filestem, reconstruction_dir, chromosomes)
    write.table(vcf.filtered, file=paste(filestem, ".filtered.vcf", sep=""), quote=F, col.names=T, row.names=F, sep="\t")
}

EOF
