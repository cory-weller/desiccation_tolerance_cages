#!/usr/bin/env bash

# argument should be iMapSplice or RSubRead
type=${1}

source ../03_map_rna/slurm.config

nSamples=$(wc -l ../03_map_rna/filestems.txt | cut -d " " -f 1)

sbatch --array=1-${nSamples}%5 <<EOF
#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 8G
#SBATCH --time 0-01:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab
#SBATCH --mail-type ALL
#SBATCH --mail-user caw5cv

module load samtools
module load gatk

filestem=\$(sed -n \${SLURM_ARRAY_TASK_ID}p ../03_map_rna/filestems.txt)

echo "Tabulating allele-specific reads with GATK ASEReadCounter"
gatk ASEReadCounter \
--input "./ASE/\${type}/\${filestem}.rg.sort.bam" \
--output "./ASE/\${type}/\${filestem}.readcounts" \
--variant ../03_map_rna/SNPs/dm3.variants \
--reference ../03_map_rna/SNPs/dm3.fa && \

echo "done"

EOF
