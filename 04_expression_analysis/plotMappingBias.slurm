#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 8G
#SBATCH --time 0-1:00:00
#SBATCH --partition largemem
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv@virginia.edu

#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(doMC)
registerDoMC(cores=5)


# # create separate filtered chromosome files to be decent size to load into R
# grep -e "^#CHROM" -e "^2L" filtered.all.vcf > filtered.2L.vcf
# grep -e "^#CHROM" -e "^2R" filtered.all.vcf > filtered.2R.vcf
# grep -e "^#CHROM" -e "^3L" filtered.all.vcf > filtered.3L.vcf
# grep -e "^#CHROM" -e "^3R" filtered.all.vcf > filtered.3R.vcf
# grep -e "^#CHROM" -e "^X" filtered.all.vcf > filtered.X.vcf


for i in $(seq 10 680); do
for i in $(seq 10 11); do

done

function getHet {
  header=($(grep -m 1 "^#CHROM" filtered.all.vcf))
  i=${1}
  sampID=${header[$i]}
  cut -f 1,2,$i filtered.all.vcf | grep -e "0/1" -e "1/0" > $sampID.het
}

export -f getHet

parallel -j 1 getHet ::: 13

iMapSplice_files <- list.files(path="./ASE/iMapSplice/", pattern="*.readcounts.gz", full.names=T)

RSubRead_files <- list.files(path="./ASE/RSubRead/", pattern="*.readcounts.gz")


importReadCounts <- function(filename) {
  rc <- fread(cmd=paste('zcat ', filename, sep=""), header=T, select=c(1,2,6,7))
}

vcf <- fread('../01_reconstructions/filtered_estimates/filtered.all.vcf')
sampleID <- "23HR-1-483"
command <- paste("bcftools query -s ", sampleID, " -f '%CHROM\t%POS[\t%GT]\n' ../01_reconstructions/filtered_estimates/filtered.all.vcf.gz")
bcftools query -s 23HR-1-483 -f '%CHROM\t%POS[\t%GT]\n' filtered.all.vcf.gz

o <- fread(cmd=command, header=FALSE, col.names=c("CHROM","POS","GT"))

dat <- fread(cmd="cut -f 1,2,20 filtered.all.vcf | grep -e '^#CHROM' -e '0/1' -e '1/0' ", header=F, col.names=c("CHROM","POS","GT"), skip="#CHROM")

get_heterozygous_sites <- function(ID) {


}

dat <- importReadCounts(iMapSplice_files[1])
dat[,total := refCount + altCount]

dat[, percent_ref := refCount/total]


tabix filtered.all.vcf.gz 2L > filtered.all.2L.vcf
tabix filtered.all.vcf.gz 2R > filtered.all.2R.vcf
tabix filtered.all.vcf.gz 2L > filtered.all.2L.vcf
tabix filtered.all.vcf.gz 2L > filtered.all.2L.vcf


dat[,c("V2","V587")][V587 %like% "0/1" | V587 %like% "1/0"]

bcftools query -s 23HR-1-483 filtered.all.vcf.gz
