#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 50G
#SBATCH -t 0-12:00:00
#SBATCH -p largemem
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv@virginia.edu

module load gcc
module load R/3.5.1

Rscript - <<EOF

library(data.table)
library(Rsubread)

expression_matrix_filename <- "dm3.counts.mat"

DESeq_filename <- "DESeq.RDS"

("Structuring expression count matrix...")
counts <- as.matrix(read.csv(expression_matrix_filename, sep="\t", row.names="GeneID", check.names=F))


# exclude bad samples (columns) with 0 feature counts for all genes

bad_samples <- melt(as.data.table(counts), measure.vars=colnames(counts))[, sum(value), by=variable][V1==0][,variable]

if(length(bad_samples) > 0) {
    print("The following sample has 0 feature counts and will be excluded: ")
    cat(paste0(bad_samples, sep="\n"))

    counts <- counts[, !(colnames(counts) %in% bad_samples)]
}

# Structuring column data table
coldata <- data.table(id=colnames(counts))
coldata[, treatment := tstrsplit(id, split="[.]")[1]]
coldata[, treatment := ifelse(treatment=="23HR", "post_starvation", "pre_starvation")]
coldata[, cage := tstrsplit(id, split="[.]")[2]]
coldata[, type := "paired-end"]
coldata[, id := NULL]
coldata <- as.data.frame(coldata)
rownames(coldata) <- colnames(counts)


library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~cage+treatment)

# Filter out rows (genes) with counts fewer than 10 across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# DESeq with "poscounts" in place of ratio for dispersion estimates,
# as we have zero-count values across all rows (genes)
#dds <- DESeq(dds, sfType="poscounts")

print("Running DESeq")

dds <- DESeq(dds)

print(paste("Saving DESeq object to", DESeq_filename))
saveRDS(dds, file=DESeq_filename)

EOF
