#!/usr/bin/env bash

source ./slurm.config

nSamples=$(wc -l filestems.txt | cut -d " " -f 1)

sbatch --array=5-${nSamples}%5 <<EOF
#!/usr/bin/env bash
#SBATCH --nodes ${mappingNodesRequest}
#SBATCH --ntasks-per-node ${mappingCoresRequest}
#SBATCH --mem ${mappingMemRequest}
#SBATCH --time ${mappingTimeRequest}
#SBATCH --partition ${mappingQueueRequest}
#SBATCH --account ${allocation}
#SBATCH --mail-type ${mail_notifications}
#SBATCH --mail-user ${userEmail}

cd ${projectDirectory}

module load singularity

# extract this run's filestem as nth line from RNAreadFileStemList
# where n is the nth run within the SLURM job array

filestem=\$(sed -n \${SLURM_ARRAY_TASK_ID}p ${RNAreadFileStemList})

echo "working on sample \${filestem}"

gunzip -c ${RNAreadsDir}/\${filestem}_R1_001.fastq.gz > ${projectDirectory}/.tmp/\${filestem}_1.fastq || \
{ echo "could not open ${RNAreadsDir}/\${filestem}_R1_001.fastq.gz"; exit 1; }

gunzip -c ${RNAreadsDir}/\${filestem}_R2_001.fastq.gz > ${projectDirectory}/.tmp/\${filestem}_2.fastq || \
{ echo "could not open ${RNAreadsDir}/\${filestem}_R1_001.fastq.gz"; exit 1; }


singularity exec ${projectDirectory}/iMapSplice.simg /iMapSplice/code/bin/iMapSplice-unphased \
-P ${projectDirectory}/SNPs/dgrp2.snps \
-Q ${projectDirectory}/SNPmerIndex/ \
-G ${projectDirectory}/globalIndex/ \
-L ${projectDirectory}/localIndex/ \
-1 ${projectDirectory}/.tmp/\${filestem}_1.fastq \
-2 ${projectDirectory}/.tmp/\${filestem}_2.fastq \
-T ${mappingCoresRequest} \
-O ${projectDirectory}/.tmp/\${filestem}/ || { echo "iMapSplice mapping failed"; rm -rf  ${projectDirectory}/.tmp/\${filestem}/; rm ${projectDirectory}/.tmp/\${filestem}*.fastq; exit 1; }

echo "cleaning up"

# Move output to separate folder
mv ${projectDirectory}/.tmp/\${filestem}/output.sam ${projectDirectory}/mapped_reads/\${filestem}.sam

# Remove temporary output folder
rm -rf ${projectDirectory}/.tmp/\${filestem}/

# Remove temporary,unzipped reads
rm ${projectDirectory}/.tmp/\${filestem}_1.fastq
rm ${projectDirectory}/.tmp/\${filestem}_2.fastq

echo "done"

EOF
