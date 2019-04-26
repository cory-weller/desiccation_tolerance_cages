#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 50G
#SBATCH --time 0-12:00:00
#SBATCH --partition largemem
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


gds_filename <- "../01_reconstruction/filtered_estimates/filtered.all.vcf.gds"

geno<-snpgdsOpen(gds_filename, allow.fork=TRUE, readonly = TRUE)
genotyped_samples<- read.gdsn(index.gdsn(geno, "sample.id"))
snpgdsClose(geno)

# Organize phenotype annotated data.frame
phenotypes <- data.table(scanID=genotyped_samples)
phenotypes[, sampleTime := tstrsplit(scanID, "-")[1]]
phenotypes[, cage := tstrsplit(scanID, "-")[2]]
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
        loco.grm <- readRDS("chr2_excluded.grm.rds")
        snpIDs <- getSnpID(geno, index=(chrom %in% c("2L","2R")))
    } else if(grepl("3", chromosome)) {
        loco.grm <- readRDS("chr3_excluded.grm.rds")
        snpIDs <- getSnpID(geno, index=(chrom %in% c("3L","3R")))
    } else if(grepl("X", chromosome)) {
        loco.grm <- readRDS("chrX_excluded.grm.rds")
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
     write.table(dat[, c("chr","pos","n.obs","Score","Score.SE","Score.Stat","Score.pval")],
        file=paste(chromosome, ".gwas", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
    remove(dat)
    remove(loco.grm)
    remove(snpIDs)
    gc()

}

}


  #if a permutation, randomly scramble sample ids within a swarm. if perm is 0, ids will stay the same




  # pass=filters[adaptive_perm_filter=="PASS", snp.id]

















  assoc.list<-foreach(i=c("2", "3", "X"))%do%{
    if(i=="2"){
      chr.i="chr2L" }
    if(i=="3"){
      chr.i="chr3L"}
    if(i=="X"){
      chr.i="chrX"}
    print(chr.i)
    if(pop!="both"){

    } else {
      nullmod <- fitNullMM(scanData = scanAnnot,
                           outcome = phenotype,
                           covMatList = loco.list[[chr.i]],
                           covars = c("temp.rack.cal", "generation", "photoperiod", "swarm"),
                           family = binomial,
                           dropZeros=F)

    }

  ### do association w/ parallelization
    chrom <- getChromosome(geno)

    if(i=="2"){
      snpIDs <- getSnpID(geno, index=(chrom=="2L"|chrom=="2R"))}
      chromosome <- "chr2"
    if(i=="3"){
      snpIDs <- getSnpID(geno, index=(chrom=="3L"|chrom=="3R"))}
    if(i=="X"){
      snpIDs <- getSnpID(geno, index=(chrom=="X"))}


    snpIDs <- intersect(snpIDs, pass)

getSnpID(geno, index=(chrom=="2L"|chrom=="2R"))


      model <- assocTestMM(genoData = genoData,
                         nullMMobj = nullmod,
                         test = "Score",
                     )
      return(as.data.table(model))
    }
    return(rbindlist(assoc))
  }








  close(geno)
