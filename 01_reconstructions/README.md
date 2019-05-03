# Reconstructions of whole genomes from low coverage DNA sequencing

Download reads and retrieve supplemental executables/libraries with `init.sh`

Map reads using `mapReads.slurm`

Run reconstructions with `reconstruct.slurm`

Calculate which regions to be filtered / excluded from downstream analysis with `cleanup.slurm`

Conduct the aforementioned filtering with `filter.slurm` to output filtered `.vcf` file
