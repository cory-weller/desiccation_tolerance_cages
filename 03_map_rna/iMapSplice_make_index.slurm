#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 50G
#SBATCH --time 0-12:00:00
#SBATCH --partition largemem
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv

module load singularity

# Part 1: generate reference index
singularity exec iMapSplice.simg /iMapSplice/code/bin/buildWholeGenome ./reference_genome/ ./globalIndex/ || { echo "buildWholeGenome failed"; exit 1; }
singularity exec iMapSplice.simg /iMapSplice/code/bin/build2ndLevelIndex ./globalIndex/ ./localIndex/ || { echo "build2ndLevelIndex failed"; exit 1; }

# Part 2: generate SNP-mers
singularity exec iMapSplice.simg /iMapSplice/code/bin/getSNPmer-unphased \
"./globalIndex/ \
"./geneAnnotationFile/dm3.imapsplice.pruned.gaf \
"./SNPs/dgrp2.snps \
"./SNPmer/ \
201 \
5 || { echo "getSNPmer-unphased failed"; exit 1; }

# Part 3: build SNP-mer Index
singularity exec iMapSplice.simg /iMapSplice/code/bin/buildSNPmerIndex \
./SNPmer/SNPinAnn.fa \
./SNPmerIndex || { echo "buildSNPmerIndex failed"; exit 1; }

echo "done!"
# mapping reads can now be done in a separate job for each sequenced individual
