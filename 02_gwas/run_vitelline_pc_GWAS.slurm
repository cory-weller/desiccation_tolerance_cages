#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 5
#SBATCH --mem 40G
#SBATCH --time 0-12:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv@virginia.edu

module load gcc
module load R/3.5.1

Rscript - <<EOF
#!/usr/bin/env Rscript

library(GWASTools, quietly=TRUE)
library(GENESIS, quietly=TRUE)
library(data.table, quietly=TRUE)
library(SNPRelate, quietly=TRUE)
library(gdsfmt, quietly=TRUE)
library(doMC, quietly=TRUE)
registerDoMC(5)
library(foreach, quietly=TRUE)

sessionInfo()

# close gds files if they are still open
showfile.gds(closeall=TRUE)


gds_filename <- "../01_reconstructions/filtered_estimates/filtered.all.vcf.gds"

geno<-snpgdsOpen(gds_filename, allow.fork=TRUE, readonly = TRUE)
genotyped_samples<- read.gdsn(index.gdsn(geno, "sample.id"))
snpgdsClose(geno)

# Organize phenotype annotated data.frame
phenotypes <- data.table(scanID=genotyped_samples)
vitelline_pc <- fread('../04_expression_analysis/DE/vitelline_pca.tab')
vitelline_pc[, scanID := gsub("\\\.", "-", rn)]
vitelline_pc <- vitelline_pc[,c("scanID","PC1")]
phenotypes <- merge(phenotypes, vitelline_pc, by.x="scanID", by.y="scanID")
phenotypes[, sampleTime := tstrsplit(scanID, "-")[1]]
phenotypes[, cage := tstrsplit(scanID, "-")[2]]
phenotypes[, pheno := ifelse(sampleTime=="23HR", 1, 0)]
phenotypes[, pheno := PC1]
phenotypes[, PC1 := NULL]

phenotypes[, sampleTime := as.factor(sampleTime)]
phenotypes[, cage := as.factor(cage)]
phenotypes <- phenotypes[sampleTime=="4CB"]
phenotypes[,sampleTime := NULL]
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
                       cov.mat = loco.grm[phenotypes\$scanID,phenotypes\$scanID],
                       covars = "cage",
                       family = "gaussian",
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
     write.table(dat[, c("chr","pos","n.obs","Score","Score.SE","Score.Stat","Score.pval")],
        file=paste(chromosome, ".vitelline.gwas", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
    remove(dat)
    remove(loco.grm)
    remove(snpIDs)
    gc()

}

close(geno)
