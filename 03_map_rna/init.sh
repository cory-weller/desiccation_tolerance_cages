#!/usr/bin/env bash

# Retrieve 15 separate chunks of tar archive (5GB file limit on Box)
wget -L https://virginia.box.com/shared/static/gziz8l406cdm4rv3x5lzwbchm59qgt0w.00 -O cDNA.tar.00
wget -L https://virginia.box.com/shared/static/vntszvl4lcgz5zk7mmakpog8quv89xix.01 -O cDNA.tar.01
wget -L https://virginia.box.com/shared/static/3g3rofgdrr840llnd8jxk2yd0aeac3bc.02 -O cDNA.tar.02
wget -L https://virginia.box.com/shared/static/5pwqe0xw3va3jyw1k29chz64b4cm2v6t.03 -O cDNA.tar.03
wget -L https://virginia.box.com/shared/static/dijv93z0gaocsmmb0ikv59ecgljz1zk8.04 -O cDNA.tar.04
wget -L https://virginia.box.com/shared/static/o9ab6qh7lizudu22qtg5b77pt2ihk1iv.05 -O cDNA.tar.05
wget -L https://virginia.box.com/shared/static/mhji937bqpmuvlycur999ie6a6y7mj9p.06 -O cDNA.tar.06
wget -L https://virginia.box.com/shared/static/232dh6eim3obejw6mvlnsqg2vais3a35.07 -O cDNA.tar.07
wget -L https://virginia.box.com/shared/static/jnfr5dam6lnkpia4v529ud1tj14ghss5.08 -O cDNA.tar.08
wget -L https://virginia.box.com/shared/static/8oyp3rxrli2d622yqbn93eqp7v4511lj.09 -O cDNA.tar.09
wget -L https://virginia.box.com/shared/static/jl4ij6f1rdz95cn7dt42yj9kv9d46pi1.10 -O cDNA.tar.10
wget -L https://virginia.box.com/shared/static/df8h8hxz5johit5nf7k8ucwpnvyrrwb1.11 -O cDNA.tar.11
wget -L https://virginia.box.com/shared/static/j5m97bw8pz5ndx8ps4k1uoacilhm3uat.12 -O cDNA.tar.12
wget -L https://virginia.box.com/shared/static/n8ip3ypco6wn13e8bieugf6hwzxgz8ce.13 -O cDNA.tar.13
wget -L https://virginia.box.com/shared/static/103b6l1djw5n9xui8enrvgv4litdicl9.14 -O cDNA.tar.14

# Combine separate parts & untar
cat cDNA.tar.* > cDNA.tar && rm cDNA.tar.* && tar -xf cDNA.tar -d ./_rna/ && rm cDNA.tar

# Generate list of filestems
cd ./rna_reads/ && ls *.fastq.gz | cut -d "_" -f 1 | sort -u > ../filestems.txt && cd ..

# Move reference genome to proper folder, reformat, and split into separate files for eahc chromosome
cp ../01_reconstructions/dm3.fasta ./reference_genome && \
cd ./reference_genome && \
sed -i 's/^>/>chr/g' dm3.fasta && \
bash ../splitFastaFiles.sh dm3.fasta && \
rm chr2LHet.fa chr2RHet.fa chr3LHet.fa chr3RHet.fa chrXHet.fa chrU.fa chrUextra.fa chrYHet.fa chr4.fa chrdmel_mitochondrion_genome.fa dm3.fasta && \
cd ..

# Retrieve singularity image file
module load singularity
singularity pull -n iMapSplice.simg shub://cory-weller/iMapSplice.simg

# Retrieve reference genome files
# wget -O refgenome.zip -L https://virginia.box.com/shared/static/4pwlmpjjzzhihm8gb4h60p5cus590x7a.zip && \
# unzip refgenome.zip -d ./reference_genome/ && rm refgenome.zip && cat ./reference_genome/*.fa > ./SNPs/dm3.fa
# Above code retrieves the pre-made version form UVA box; below version recreates it from the files already extant within directory
cat ./reference_genome/chr2L.fa ./reference_genome/chr2R.fa ./reference_genome/chr3L.fa ./reference_genome/chr3R.fa ./reference_genome/chrX.fa > ./SNPs/dm3.fa
module load samtools
module load gatk
if [ ! -f ./SNPs/dm3.fa.fai ]; then
  gatk CreateSequenceDictionary --REFERENCE ./SNPs/dm3.fa
  samtools faidx ./SNPs/dm3.fa
fi

# Retrieve dgrp SNP table (for iMapSplice mapping)
wget -O dgrp2.snps.gz -L https://virginia.box.com/shared/static/5ia84k3fc531e5f0rzxmxrs8ayc5zwre.gz && \
gunzip -c dgrp2.snps.gz > ./SNPs/dgrp2.snps && rm dgrp2.snps.gz

# Retrieve dgrp variant table (for ASEReadCounter)
wget -O dm3.variants.gz -L https://virginia.box.com/shared/static/tbqyp2j8e5nds14s78hhrcg7bdpauffc.gz && \
gunzip -c dm3.variants.gz | sed 's/FORMAT$/FORMAT\thet/g'  | sed 's/GT$/GT\t0\/1/g' > ./SNPs/dm3.variants && rm dm3.variants.gz
if [ ! -f ./SNPs/dm3.variants.idx ]; then
  gatk IndexFeatureFile --feature-file .//SNPs/dm3.variants
fi

# Generate iMapSplice format gene annotation file
zcat ../dm3.genepred.gz | tail -n +2 | cut -f 2- > dm3.genepred
sed 's/$/\tINFO1\tINFO2/g' dm3.genepred > dm3.imapsplice.gaf
bash ./pruneGAF.sh dm3.imapsplice.gaf && mv dm3.imapsplice.pruned.gaf ./geneAnnotationFile/

# Build Rsubread index
sbatch  <<EOF
#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 20G
#SBATCH --time 0-00:10:00
#SBATCH --partition largemem
#SBATCH --account berglandlab

cd ./SNPs/
module load R/3.5.1
Rscript - <<EOF
#!/usr/bin/env Rscript
library(Rsubread)
buildindex(basename="dm3_index", reference="dm3.fa")
EOF
