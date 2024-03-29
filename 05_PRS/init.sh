#!/usr/bin/env bash

# Retrieve PRSice
wget https://github.com/choishingwan/PRSice/releases/download/2.1.11/PRSice_linux.zip
unzip PRSice_linux.zip PRSice.R
unzip PRSice_linux.zip PRSice_linux


# Retrieve PLINK
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip
unzip plink_linux_x86_64_20190304.zip plink


# Generate PCs as covariates
Rscript - <<EOF
library(data.table)
dat <- readRDS('chrX.LOCO.grm.rds')
vec <- eigen(as.matrix(dat))
PCS <- data.table(vec$vectors[,1:20])
setnames(PCS, paste("PC", 1:20, sep=""))

PCS[, IID := rownames(dat)]
PCS[, FID := rownames(dat)]
setcolorder(PCS, c("FID", "IID", paste("PC",1:20, sep="")))
fwrite(PCS, file="PCS.cov", quote=F, sep="\t", row.names=F, col.names=T)
EOF

# Create nice chromosome codes
sed 's/2L/1/g' ../01_reconstructions/post_filter/filtered.all.vcf | sed 's/2R/2/g' | sed 's/3L/3/g' | sed 's/3R/4/g' | sed 's/X/5/g' > convertedChrs.vcf

# Fill  ID column with chrom_bp (instead of NA)
sed -E -i "s/(^[0-9])\t([0-9]+)\tNA/\1\t\2\t\1_\2/g" convertedChrs.vcf


# Convert vcf to plink binary, excluding X chromosome (not supported by PRSice)
./plink --vcf convertedChrs.vcf --make-bed --allow-extra-chr --vcf-half-call missing --out plinkFmt

# Edit FAM file
Rscript - <<EOF
library(data.table)
fam <- fread('plinkFmt.fam')
setnames(fam, c("FID", "IID", "IIDF", "IIDM", "SEX", "PHENO"))
fam[FID %like% "23HR", PHENO := 2]
fam[FID %like% "4CB", PHENO := 1]
fam[, SEX := 2]
fwrite(fam, file="plinkFmt.fam", quote=F, sep="\t", row.names=F, col.names=F)
EOF


# Create Sites table fo radding A1 and A2 info later
cut -f 1,2,4,5 convertedChrs.vcf > sites.vcf

# Genesis reports log odds ratio; create odds ratio by 10^(Score)
Rscript - <<EOF
library(data.table)
gwas <- fread('../02_gwas/allChr.gwas')

setnames(gwas, "Score", "LogOR")
setnames(gwas, "Score.SE", "OR.se")
setnames(gwas, "Score.Stat", "Z")
setnames(gwas, "Score.pval", "P")
setnames(gwas, "chr", "CHR")
setnames(gwas, "pos", "POS")

gwas[CHR == "2L", CHR := "1"]
gwas[CHR == "2R", CHR := "2"]
gwas[CHR == "3L", CHR := "3"]
gwas[CHR == "3R", CHR := "4"]
gwas[CHR == "X", CHR := "5"]

# Convert Genesis log-odds ratio to
gwas[, OR := exp(LogOR)]
gwas[, SNP_ID := paste(CHR, POS, sep="_")]


# Merge in ref/alt alleles
sites <- fread('sites.vcf', skip="#CHROM")
setnames(sites, "#CHROM", "CHROM")

sites[, SNP_ID := paste(CHROM, POS, sep="_")]
sites[, c("CHROM", "POS") := NULL]
setkey(sites, SNP_ID)
setkey(gwas, SNP_ID)
dat.all <- merge(sites, gwas)

# Reorder by chromosome then position
setkey(dat.all, CHR, POS)

# Filter low observations
dat.all <- dat.all[n.obs >= 300]
setnames(dat.all, "REF", "A1")
setnames(dat.all, "ALT", "A2")
setnames(dat.all, "POS", "BP")
setnames(dat.all, "SNP_ID", "SNP")

# Write output, excluding X chromosome
fwrite(dat.all[, c("SNP","CHR","BP","A1","A2","P","OR","Z")][CHR != "5"], file="BASE.assoc", sep=" ", quote=F, col.names=T, row.names=F)
EOF

#cat <(egrep '^4CB' plinkFmt.fam | cut -f 1 | awk '{OFS="\t";}{print $1,"1"}') <(egrep '^23HR' plinkFmt.fam | cut -f 1 | awk '{OFS="\t";}{print $1,"2"}') | sort | awk '{print $1,$1,$2}' > PRSice.phenos

Rscript ./PRSice.R --dir . \
    --prsice ./PRSice_linux \
    --base BASE.assoc \
    --target plinkFmt \
    --thread 4 \
    --stat OR \
    --binary-target T \
    --clump-kb 100 \
    --cov-file PCS.cov \
    --cov-col @PC[1-20] \
    --print-snp

# Subset only RNA reads individuals

egrep -v "NA$" plinkFmt.fam |cut -d " " -f 1 > samplefile

# Add the following to the header:
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
#
sed -i 's/NA/./g' convertedChrs.vcf

module load bcftools
module load samtools

bgzip convertedChrs.vcf
tabix -p vcf convertedChrs.vcf.gz




bcftools view -S samplefile convertedChrs.vcf.gz > subset.vcf

./plink --vcf subset.vcf --make-bed --allow-extra-chr --vcf-half-call missing --out subset


Rscript ./PRSice.R --dir . \
    --prsice ./PRSice_linux \
    --base BASE.assoc \
    --target subset \
    --stat Z \
    --binary-target F \
    --beta \
    --thread 4 \
    --clump-kb 100 \
    --print-snp

Rscript - <<EOF
#!/usr/bin/env Rscript
library(ggplot2)
library(ggthemes)

# Load in output of PRSice
prs <- fread('PRSice.best')
prs[, group := tstrsplit(FID, split="-")[1]]
prs[, cage := tstrsplit(FID, split="-")[2]]
prs[, treatment := ifelse(group=="4CB", "Pre-starvation", "Post-starvation")]
prs[, outcome := ifelse(treatment=="Pre-starvation", 0, 1)]

# Load in expression phenotypes
pheno <- fread('subset.fam')

merged <- merge(pheno, prs, by.x="V1", by.y="FID")




logitModel <- glm(outcome ~ PRS, data=merged, family=binomial(link="logit"))
predicted <- plogis(predict(logitModel, merged))
merged[, predicted := predicted]
PRSs <- merged$PRS

merged[, percentile := ecdf(PRSs)(PRS), by=V1]
ggplot(merged, aes(x=factor(outcome), y=percentile)) + geom_boxplot() + geom_jitter(alpha=0.3)

g1 <- ggplot(merged, aes(x=PRS, y=V6, color=treatment, shape=cage)) + geom_point() + theme_few(16) + labs(x="Polygenic Risk Score", y="PC1 of Expression Data")
g2 <- modelPrediction <- ggplot(merged[,c("PRS","predicted")], aes(x=PRS, y=predicted)) + geom_line() + theme_few(16) + labs(x="Polygenic Risk Score", y="Predicted Treatment")
ggsave(g1, file="PRS.svg")
ggsave(g2, file="model.svg")

EOF
