#!/usr/bin/env bash

# Remove bin and header row
zcat ../dm3.genepred.gz | cut -f 2- | tail -n +2 > dm3.genepred

# Generate gtf file from dm3.genepred:
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
chmod +x ./genePredToGtf
./genePredToGtf file dm3.genepred dm3.gtf

# Fix gene IDs to be identical across isoforms
sed -i 's/\(gene_id "CG[^-]*\)-R./\1/g' dm3.gtf
