#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 24G
#SBATCH --time 0-1:30:00
#SBATCH --partition standard
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv@virginia.edu

module load gcc
module load R/3.5.1

iteration=$SLURM_ARRAY_TASK_ID

Rscript - $iteration <<EOF
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
iteration <- args[1]
print(iteration)

library(GWASTools, quietly=TRUE)
library(GENESIS, quietly=TRUE)
library(data.table, quietly=TRUE)
library(SNPRelate, quietly=TRUE)
library(gdsfmt, quietly=TRUE)
library(foreach, quietly=TRUE)
library(doMC, quietly=TRUE)
registerDoMC(core=4)

sessionInfo()

# close gds files if they are still open
showfile.gds(closeall=TRUE)


gds_filename <- "../01_reconstructions/filtered_estimates/filtered.all.vcf.gds"

geno<-snpgdsOpen(gds_filename, allow.fork=TRUE, readonly = TRUE)
genotyped_samples<- read.gdsn(index.gdsn(geno, "sample.id"))
snpgdsClose(geno)

# Organize phenotype annotated data.frame
phenotypes <- data.table(scanID=genotyped_samples)
phenotypes[, sampleTime := tstrsplit(scanID, "-")[1]]
phenotypes[, cage := tstrsplit(scanID, "-")[2]]


# permute sample time within cage
phenotypes[, sampleTimePermuted := sample(sampleTime, size=.N, replace=F), by=cage]
phenotypes[, sampleTime := sampleTimePermuted]
phenotypes[, sampleTimePermuted := NULL]

phenotypes[, pheno := ifelse(sampleTime=="23HR", 1, 0)]

phenotypes[, sampleTime := as.factor(sampleTime)]
phenotypes[, cage := as.factor(cage)]
scanAnnot <- ScanAnnotationDataFrame(phenotypes)


geno <- GdsGenotypeReader(filename = gds_filename)
genoData <- GenotypeData(geno)
chrom <- getChromosome(geno)
chrom_variant_key <- data.table(variant.id=getSnpID(geno), chr=getChromosome(geno))
setkey(chrom_variant_key, variant.id)


for(chromosome in c("chr2", "chr3", "chrX")) {

    print(chromosome)

    # Exclude problematic individual "4CB-3-426" and cage 4 (no paired data)
    ids.to.use <- phenotypes[cage != 4 & scanID != "4CB-3-426",scanID]

    # load LOCO GRM and SNP ids for this chromosome
    if(grepl("2", chromosome)) {
        loco.grm <- readRDS("../01_reconstructions/GRM/chr2_excluded.GRM.rds")
        snpIDs <- getSnpID(geno, index=(chrom %in% c("2L","2R")))
    } else if(grepl("3", chromosome)) {
        loco.grm <- readRDS("../01_reconstructions/GRM/chr3_excluded.GRM.rds")
        snpIDs <- getSnpID(geno, index=(chrom %in% c("3L","3R")))
    } else if(grepl("X", chromosome)) {
        loco.grm <- readRDS("../01_reconstructions/GRM/chrX_excluded.GRM.rds")
        snpIDs <- getSnpID(geno, index=(chrom %in% c("X")))
    }


    nullmod <- fitNullModel(scanAnnot,
                       outcome = "pheno",
                       cov.mat = loco.grm,
                       covars = "cage",
                       family = binomial,
                       sample.id = ids.to.use,
                   drop.zeros=F)

    iterator <- GenotypeBlockIterator(genoData, snpInclude=snpIDs)

    model <- assocTestSingle(gdsobj = iterator,
                         null.model = nullmod,
                         test = "Score")

     dat <- as.data.table(model)[!is.na(Score)]
     dat[, chr := NULL]
     setkey(dat, variant.id)
     dat <- merge(dat, chrom_variant_key)
    gz1 <- gzfile(paste("permutations/", iteration, ".", chromosome, ".gwas.gz", sep=""), "w")
     write.table(dat[n.obs >= 400][, c("chr","pos","n.obs","Score.pval")], gz1, row.names=F, col.names=F, sep="\t", quote=F, append=T)
    close(gz1)
    remove(dat)
    remove(loco.grm)
    remove(snpIDs)
    gc()
}

close(geno)
EOF
