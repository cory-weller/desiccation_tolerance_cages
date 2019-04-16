#!/usr/bin/env bash

source slurm.config


sbatch <<EOF

#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 8G
#SBATCH --time 0-24:00:00
#SBATCH --partition standard
#SBATCH --account ${allocation}

module load samtools/1.9

cd mapped_reads/

while read filename; do
  samtools view -b \${filename} | samtools sort - > \${filename%.sam}.bam && samtools index \${filename%.sam}.bam
done < <(ls *.sam)

EOF
