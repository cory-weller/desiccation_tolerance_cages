#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 12G
#SBATCH --time 0-05:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab


module load gcc
module load R/3.5.1

Rscript - <<EOF
#!/usr/bin/env Rscript

library(GWASTools)
library(data.table)
library(SNPRelate)
library(gdsfmt)

gds_filename <- "filtered_estimates/filtered.all.vcf.gds"


print("opening GDS!")
geno <- snpgdsOpen(gds_filename, allow.fork=TRUE, readonly = TRUE)

#read genotypes and filters
print("filtering snps")

genotyped_samples <- read.gdsn(index.gdsn(geno, "sample.id"))


snpset <- snpgdsLDpruning(geno,
                          ld.threshold=0.2,
                          slide.max.bp = 5000,
                          autosome.only=FALSE,
                          missing.rate=.15,
                          maf=.05)


# Make chr2 matrix
chr.to.use <- c("chr2L", "chr2R")
print(paste0("using chromosomes: ", paste0(chr.to.use, collapse=", ")))
snpset.use <- unlist(snpset[chr.to.use])
print(paste0("nSNPs: ", length(snpset.use)))
grm_raw <- snpgdsGRM(geno, method="Eigenstrat", snp.id=snpset.use, sample.id=genotyped_samples, autosome.only=FALSE, num.thread=1)
grm <- grm_raw\$grm
rownames(grm) <- grm_raw\$sample.id
colnames(grm) <- grm_raw\$sample.id

saveRDS(grm, file="chr2_only.GRM.rds")
remove(grm)
remove(grm_raw)

# Make chr2-LOCO matrix
chr.to.use <- c("chr3L", "chr3R", "chrX")
print(paste0("using chromosomes: ", paste0(chr.to.use, collapse=", ")))
snpset.use <- unlist(snpset[chr.to.use])
print(paste0("nSNPs: ", length(snpset.use)))
grm_raw <- snpgdsGRM(geno, method="Eigenstrat", snp.id=snpset.use, sample.id=genotyped_samples, autosome.only=FALSE, num.thread=1)
grm <- grm_raw\$grm
rownames(grm) <- grm_raw\$sample.id
colnames(grm) <- grm_raw\$sample.id

saveRDS(grm, file="chr2_excluded.GRM.rds")
remove(grm)
remove(grm_raw)

# Make chr3 matrix
chr.to.use <- c("chr3L", "chr3R")
print(paste0("using chromosomes: ", paste0(chr.to.use, collapse=", ")))
snpset.use <- unlist(snpset[chr.to.use])
print(paste0("nSNPs: ", length(snpset.use)))
grm_raw <- snpgdsGRM(geno, method="Eigenstrat", snp.id=snpset.use, sample.id=genotyped_samples, autosome.only=FALSE, num.thread=1)
grm <- grm_raw\$grm
rownames(grm) <- grm_raw\$sample.id
colnames(grm) <- grm_raw\$sample.id

saveRDS(grm, file="chr3_only.GRM.rds")
remove(grm)
remove(grm_raw)

# Make chr3-LOCO matrix
chr.to.use <- c("chr2L", "chr2R", "chrX")
print(paste0("using chromosomes: ", paste0(chr.to.use, collapse=", ")))
snpset.use <- unlist(snpset[chr.to.use])
print(paste0("nSNPs: ", length(snpset.use)))
grm_raw <- snpgdsGRM(geno, method="Eigenstrat", snp.id=snpset.use, sample.id=genotyped_samples, autosome.only=FALSE, num.thread=1)
grm <- grm_raw\$grm
rownames(grm) <- grm_raw\$sample.id
colnames(grm) <- grm_raw\$sample.id

saveRDS(grm, file="chr3_excluded.GRM.rds")
remove(grm)
remove(grm_raw)


# Make chrX matrix
chr.to.use <- c("chrX")
print(paste0("using chromosomes: ", paste0(chr.to.use, collapse=", ")))
snpset.use <- unlist(snpset[chr.to.use])
print(paste0("nSNPs: ", length(snpset.use)))
grm_raw <- snpgdsGRM(geno, method="Eigenstrat", snp.id=snpset.use, sample.id=genotyped_samples, autosome.only=FALSE, num.thread=1)
grm <- grm_raw\$grm
rownames(grm) <- grm_raw\$sample.id
colnames(grm) <- grm_raw\$sample.id

saveRDS(grm, file="chrx_only.GRM.rds")
remove(grm)
remove(grm_raw)

# Make chrX-LOCO matrix
chr.to.use <- c("chr2L", "chr2R", "chr3L", "chr3R")
print(paste0("using chromosomes: ", paste0(chr.to.use, collapse=", ")))
snpset.use <- unlist(snpset[chr.to.use])
print(paste0("nSNPs: ", length(snpset.use)))
grm_raw <- snpgdsGRM(geno, method="Eigenstrat", snp.id=snpset.use, sample.id=genotyped_samples, autosome.only=FALSE, num.thread=1)
grm <- grm_raw\$grm
rownames(grm) <- grm_raw\$sample.id
colnames(grm) <- grm_raw\$sample.id

saveRDS(grm, file="chrX_excluded.GRM.rds")
remove(grm)
remove(grm_raw)


# Make allChr matrix

chr.to.use <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")
print(paste0("using chromosomes: ", paste0(chr.to.use, collapse=", ")))
snpset.use <- unlist(snpset[chr.to.use])
print(paste0("nSNPs: ", length(snpset.use)))
grm_raw <- snpgdsGRM(geno, method="Eigenstrat", snp.id=snpset.use, sample.id=genotyped_samples, autosome.only=FALSE, num.thread=1)
grm <- grm_raw\$grm
rownames(grm) <- grm_raw\$sample.id
colnames(grm) <- grm_raw\$sample.id

saveRDS(grm, file="chr_all.GRM.rds")
remove(grm)
remove(grm_raw)


snpgdsClose(geno)
