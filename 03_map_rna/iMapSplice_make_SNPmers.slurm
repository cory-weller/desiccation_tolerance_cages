#!/usr/bin/env bash

source ./slurm.config

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --nodes ${indexingNodesRequest}
#SBATCH --ntasks-per-node ${indexingCoresRequest}
#SBATCH --mem ${indexingMemRequest}
#SBATCH --time ${indexingTimeRequest}
#SBATCH --partition ${indexingQueueRequest}
#SBATCH --account ${allocation}
#SBATCH --mail-type ${mail_notifications}
#SBATCH --mail-user ${userEmail}

cd ${projectDirectory}

module load singularity

singularity exec iMapSplice.simg /iMapSplice/code/bin/getSNPmer-unphased \
"${projectDirectory}"/globalIndex/ \
"${projectDirectory}"/geneAnnotationFile/${gene_annotation_file_name} \
"${projectDirectory}"/SNPs/${snp_file_name} \
"${projectDirectory}"/SNPmer/ \
${snpmer_seq_length} \
${max_snpmer_samples} || { echo "getSNPmer-unphased failed"; exit 1; }

singularity exec iMapSplice.simg /iMapSplice/code/bin/buildSNPmerIndex \
"${projectDirectory}"/SNPmer/SNPinAnn.fa \
"${projectDirectory}"/SNPmerIndex || { echo "buildSNPmerIndex failed"; exit 1; }

# mapping reads can now be done in a separate job for each sequenced individual

EOF
