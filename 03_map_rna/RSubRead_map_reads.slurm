#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 5
#SBATCH --mem 40G
#SBATCH --time 0-12:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv

module load R/3.5.1

Rscript - <<EOF
#!/usr/bin/env Rscript

library(Rsubread)

# Index pre-built during init.sh
# buildindex(basename="dm3_index", reference="./SNPs/dm3.fa")

args <- commandArgs(trailingOnly=TRUE)
file_id <- args[1]

readDir <- "./rna_reads/"
indexLoc <- "./SNPs/dm3_index"

filestems <- readLines('filestems.txt')
print(filestems)

for(filestem in filestems) {
  read1 <- paste("./rna_reads/", filestem, "_R1_001.fastq.gz", sep="")
  read2 <- paste("./rna_reads/", filestem, "_R2_001.fastq.gz", sep="")
  align(index=indexLoc, readfile1=read1, readfile2=read2, input_format="gzFASTQ", output_format="BAM", output_file=paste("./mapped_reads/RSubRead/", filestem, ".bam", sep=""), nthreads=5)
  command <- paste("rm ./mapped_reads/RSubRead/", filestem, ".bam.indel.vcf", , sep="")
  system(command)
}

EOF
