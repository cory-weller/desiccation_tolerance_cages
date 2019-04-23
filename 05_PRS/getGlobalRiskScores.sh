#!/usr/bin/env bash

"/mnt/icy_1/inbredLines/metaData/popInfo.unix.csv"

vcf="inbred.ug.filter.vsqr.recode.noREP.noINDEL.vcf"

# convert to nice (numeric) chromosome codes
sed 's/^2L/1/g' $vcf | sed 's/^2R/2/g' | sed 's/^3L/3/g' | sed 's/^3R/4/g' | sed 's/^X/5/g' > globalLines.vcf

# Fill  ID column with chrom_bp (instead of .)
sed -E -i "s/(^[0-9])\t([0-9]+)\t./\1\t\2\t\1_\2/g" globalLines.vcf

grep -v 1_10503132 globalLines.vcf > globalLines.fixed.vcf
egrep -v "[ACTG],[ACTG]" globalLines.fixed.vcf > globalLines.biallelic.vcf

../plink --vcf globalLines.biallelic.vcf --make-bed --allow-extra-chr --vcf-half-call missing --out globalLines --const-fid

# One line was missing data... exclude it
# Error: Line 393406 of .vcf file has fewer tokens than expected.
# In fread("globalLines.vcf") :
#  Discarded single-line footer: <<1     10503132        1_10503132      A       T       752.88  PASS    AC=3;AF=4.348e-03;AN=690;BaseQRankSum=-2.086;DP=5024;Dels=0.00;FS=0.000;HaplotypeScore=6.4370;InbreedingCoeff=0.1198;MLEAC=2;MLEAF=2.899e-03;MQ=59.98;MQ0=1;MQRankSum=-2.056;QD=27.88;ReadPosRankSum=0.070;SOR=0.871;VQSLOD=4.34;culprit=MQ     GT:AD:DP:GQ:PL  ./.:.:.:.:.     0/0:4,0:4:12:0,12,157   0/0:9,0:9:27:0,27,351   0/0:17,0:17:51:0,51,657 0/0:4,0:4:12:0,12,151   0/0:8,0:8:24>>
# grep -n 1_10503132 = line 393406

../plink --vcf globalLines.fixed.vcf --make-bed --allow-extra-chr --vcf-half-call missing --out globalLines --const-fid


../plink --bfile plinkFmt --bmerge globalLines --make-bed --out merged


# replace -9 with pre-starvation "2" phenotype in .fam file
Rscript - <<EOF
library(data.table)

dat <- fread('globalLines.fam')
dat[, V6 := "2"]
dat <- dat[,c("V2", "V2", "V3", "V4", "V5", "V6")]
fwrite(dat, file="globalLines.fam", quote=F, sep=" ", row.names=F, col.names=F)
EOF

Rscript ../PRSice.R --dir . \
    --prsice ../PRSice_linux \
    --base ../BASE.assoc \
    --target globalLines \
    --thread 4 \
    --stat OR \
    --binary-target T \
    --clump-kb 100 \
    --print-snp
