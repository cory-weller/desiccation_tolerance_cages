#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 20G
#SBATCH --time 0-01:00:00
#SBATCH --partition largemem
#SBATCH --account berglandlab

minimum_BaseMean=1
minimum_l2fc=0

minimum_MAF="0.1"
minimum_fraction_genotyped="0.8"

# Print header to single file for each chromosome
grep -m 1 -e "^#CHROM" ../../01_reconstructions/filtered_estimates/filtered.all.vcf | sed 's/#//g' | tee "2L.vcf" "2R.vcf" "3L.vcf" "3R.vcf" "X.vcf"

# Split VCF into separate chromosomes
awk -F "\t" '{if($1!~/^#/) {s=$1".vcf";  print >> s}}' ../../01_reconstructions/filtered_estimates/filtered.all.vcf

# Get filtered gene IDs
Rscript - ${minimum_BaseMean} ${minimum_l2fc} <<EOF
#!/usr/bin/env Rscript
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
minimum_BaseMean <- as.numeric(args[1])
minimum_l2fc <- as.numeric(args[2])

DEseq_output <- fread('../DE/dds.shrink.dat')
pass_filter <- DEseq_output[baseMean >= minimum_BaseMean & abs(log2FoldChange) > minimum_l2fc]

writeLines(pass_filter[,FlyBaseID], con="filtered_gene_ids.txt")
EOF

# Format gene expression file
Rscript - <<EOF
#!/usr/bin/env Rscript

chromosomes <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(foreach)

dat <- fread('../DE/dm3.featureCountMatrix')
dat <- dat[GeneID != ""]

# convert "." in column names to "-"
setnames(dat, gsub("[.]", "-", colnames(dat)))

lengths <- fread('../DE/dm3.genelength')
gc_to_fb <- fread('../dm3.GC_to_FB.txt', header=TRUE, col.names=c("GC_Annotation_ID", "current_id", "FB_id", "gene_symbol"), key=c("GC_Annotation_ID"))

# add FB gene IDs
lengths <- merge(lengths, gc_to_fb, by.x="AnnotationID", by.y="GC_Annotation_ID")
lengths <- lengths[, c("FB_id","gene_symbol","Length")]

dat <- merge(dat, lengths, by.x="GeneID", by.y="FB_id")

# get sample ID names and write to file
sampleNames <- c(colnames(dat)[grepl("4CB", colnames(dat))], colnames(dat)[grepl("23HR", colnames(dat))])
sampleNames <- data.table(id=sampleNames)
sampleNames[, n := tstrsplit(id, split="-")[3]]
sampleNames[, n := as.numeric(n)]
setkey(sampleNames, n)
writeLines(sampleNames[,id], con="samples.txt")

geneLengths <- dat[,Length]

TPM_matrix <- foreach(sample=sampleNames[,id], .combine="cbind") %do% {
    counts <- dat[[sample]]
    rates <- counts / geneLengths
    TPM <- data.table((rates / sum(rates)) * 1e6)
    setnames(TPM, sample)
    return(TPM)
}

TPM_matrix <- cbind("id" = dat[,GeneID], TPM_matrix

# load filtered genes (by baseMean expression)
pass_filter <- readLines('filtered_gene_ids.txt')
TPM_matrix <- TPM_matrix[GeneID %in% pass_filter]

fwrite(TPM_matrix, file="TPM_expression.txt", quote=F, na="NA", col.names=T, row.names=F, sep="\t")
EOF


# Generate Genotypes file
Rscript - ${minimum_MAF} ${minimum_fraction_genotyped} <<EOF
#!/usr/bin/env Rscript
library(data.table)

chromosomes <- c("2L","2R","3L","3R","X")

args <- commandArgs(trailingOnly=TRUE)

minimum_MAF <- as.numeric(args[1])
minimum_fraction_genotyped <- as.numeric(args[2])

convert_vcf_to_ref_dosage <- function(DT, sampleNames) {
  cols <- colnames(DT)[3:ncol(DT)]
  for(col in cols) set(DT, i=which(DT[[col]]=="0/0"), j=col, value="2")
  for(col in cols) set(DT, i=which(DT[[col]]=="1/0"), j=col, value="1")
  for(col in cols) set(DT, i=which(DT[[col]]=="0/1"), j=col, value="1")
  for(col in cols) set(DT, i=which(DT[[col]]=="1/1"), j=col, value="0")

  # filter by missing data
  DT[, N_missing_genotype := apply(.SD, 1, function(x) sum(is.na(x))), .SD=sampleNames]

  # filter by MAF

}

sampleNames <- readLines('samples.txt')

for(chromosome in chromosomes) {
  inFile <- paste(chromosome, ".vcf", sep="")
  outFile <- paste(chromosome, ".genotypes", sep="")
  nCol <- as.numeric(system(paste("head -n 1 ", inFile, " | wc -w", sep=""), intern=TRUE))

  # Read in VCF file for given chromosome, converting any genotype with "." into NA
  dat <- fread(inFile, na.strings=c("./.", "./0", "./1", "0/.", "1/."), select=c("CHROM", "POS", sampleNames), colClasses=rep(c("character"), nCol))

  # Convert VCF genotype code into reference allele dosage
  convert_vcf_to_ref_dosage(dat, sampleNames)

  # Convert sample IDs
  ids <- paste(dat[,CHROM], dat[,POS], sep="_")
  dat[, c("POS","CHROM") := NULL]

  fwrite(cbind("id"=ids, dat), file=outFile, quote=F, sep="\t", col.names=T, row.names=F, na="NA")
  remove(dat)
  gc()
}
EOF


# Get filtered SNPs by minimum MAF
Rscript - ${minimum_MAF} ${chromosome} <<EOF
#!/usr/bin/env Rscript
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
minimum_MAF <- as.numeric(args[1])
chromosome <- args[1]

genotypes <- fread(paste(chromosome, ".


DEseq_output <- fread('../DE/dds.shrink.dat')
pass_filter <- DEseq_output[baseMean >= minimum_BaseMean & abs(log2FoldChange) > minimum_l2fc]

writeLines(pass_filter[,FlyBaseID], con="filtered_gene_ids.txt")
EOF


# Generate gene location (annotation) file
Rscript - <<EOF
#!/usr/bin/env Rscript
library(data.table)

gtf <- fread('../dm3.gtf')
gtf <- gtf[V3=="transcript"]
gtf[, AnnotationID := tstrsplit(V9, '"')[2]]
gtf <- gtf[, list("s1"=min(V4), "s2"=max(V5)), by=list(V1, AnnotationID)]
gtf[, chr := tstrsplit(V1, "hr")[2]]
gtf <- gtf[chr %in% c("2L","2R","3L","3R","X")]

gc_to_fb <- fread('../dm3.GC_to_FB.txt', header=TRUE, col.names=c("GC_Annotation_ID", "current_id", "FB_id", "gene_symbol"), key=c("GC_Annotation_ID"))

dat <- merge(gtf, gc_to_fb, by.x="AnnotationID", by.y="GC_Annotation_ID")
setnames(dat, "FB_id", "geneid")
dat.out <- dat[, c("geneid", "chr", "s1", "s2")]

pass_filter <- readLines('filtered_gene_ids.txt')
dat.out <- dat.out[geneid %in% pass_filter]

for(chromosome in c("2L","2R","3L","3R","x")) {
  outFile <- paste(chromosome, ".geneloc.txt", sep="")
  fwrite(dat.out[chr==chromosome], outFile, quote=F, row.names=F, col.names=T, sep="\t")
}
EOF


# Generate SNP file
Rscript - 2L 2R 3L 3R X <<EOF
#!/usr/bin/env Rscript
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

for(chromosome in args) {
  inFile <- paste(chromosome, ".vcf", sep="")
  outFile <- paste(chromosome, ".snps", sep="")
  dat <- fread(cmd=paste("cut -f1,2 ", inFile, sep=""))
  dat[, snp := paste(CHROM, POS, sep="_")]
  setnames(dat, "CHROM", "chr")
  setnames(dat, "POS", "pos")
  setcolorder(dat, c("snp","chr","pos"))
  fwrite(dat, outFile, quote=F, row.names=F, col.names=T, sep="\t")
}
EOF


# Generate covariates file
Rscript - <<EOF
#!/usr/bin/env Rscript

dt <- data.table(sample=readLines('samples.txt'))
dt[, treatment := ifelse(grepl("4CB", sample), "0", "1")]
dt[, cage := tstrsplit(sample, "-")[2]]

dt.out <- cbind(var0=c("id","treatment","cage"), t(dt))

fwrite(dt.out, file="covariates.txt", row.names=FALSE, col.names=FALSE, quote=F, na="NA", sep="\t")
EOF

# Save order of sample IDs
