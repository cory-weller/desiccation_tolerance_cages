#!/usr/bin/env bash

source slurm.config

nSamples=$(wc -l filestems.txt | cut -d " " -f 1)

sbatch --array=1-${nSamples}%5 <<EOF
#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 8G
#SBATCH --time 0-01:00:00
#SBATCH --partition ${mappingQueueRequest}
#SBATCH --account ${allocation}
#SBATCH --mail-type ${mail_notifications}
#SBATCH --mail-user ${userEmail}

module load samtools
module load gatk

filestem=\$(sed -n \${SLURM_ARRAY_TASK_ID}p ${RNAreadFileStemList})

gatk AddOrReplaceReadGroups \
--I=mapped_reads/\${filestem}.bam \
--O=mapped_reads/\${filestem}.rg.bam \
--RGID=\${filestem} \
--RGLB=indoorCages \
--RGPL=illumina \
--RGPU=indoorCages \
--RGSM=\${filestem}

echo "Sorting and indexing .bam"
samtools sort mapped_reads/\${filestem}.rg.bam > mapped_reads/\${filestem}.rg.sort.bam && \
samtools index mapped_reads/\${filestem}.rg.sort.bam


if [ ! -f ../SNPs/dm3.variants.idx ]; then
  gatk IndexFeatureFile --feature-file SNPs/dm3.variants
fi

if [ ! -f ../SNPs/dm3.fa.fai ]; then
  gatk CreateSequenceDictionary --REFERENCE SNPs/dm3.fa
  samtools faidx SNPs/dm3.fa
fi

echo "Tabulating allele-specific reads with GATK ASEReadCounter"
gatk ASEReadCounter \
--input mapped_reads/\${filestem}.rg.sort.bam \
--output ../04_expression_analysis/ASE/\${filestem}.readcounts \
--variant SNPs/dm3.variants \
--reference SNPs/dm3.fa && \

if [ -f mapped_reads/\${filestem}.rg.sort.bam ]; then
  rm mapped_reads/\${filestem}.rg.bam mapped_reads/\${filestem}.bam mapped_reads/\${filestem}.bam.bai
fi

echo "done"

EOF
