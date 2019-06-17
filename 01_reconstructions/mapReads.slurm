#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-1:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

# window and step size for HARP
window_step=100000
window_width=100000

topDirectory=$(pwd)

readDirectory="./dna_reads"

referenceGenome=${topDirectory}/data_external/all_dmel.fasta
gatk=${topDirectory}/etc/GenomeAnalysisTK.jar
picard=${topDirectory}/etc/picard.jar
harp=${topDirectory}/etc/harp

# If no file of unique identifying IDs, create it
if [[ ! -f filestems.txt ]]; then
    ls ${readDirectory}/*.fastq.gz | xargs -n1 basename | cut -d "_" -f 1 | sort -u > filestems.txt
fi


# fileID="23HR-1-481_S289"
fileID=$(sed -n ${SLURM_ARRAY_TASK_ID}p filestems.txt)

chromosomes=( 2L 2R 3L 3R X )



module load gcc
module load R/3.5.1

# map reads with Rsubread
Rscript - <<EOF
#!/usr/bin/env Rscript

library(Rsubread)

# index already built
# buildindex(basename="dmel", reference="../data_external/all_dmel.fasta")

file_id <- "${fileID}"

readDir <- "${readDirectory}"
indexLoc <- "${topDirectory}/Proj_Map_Reads/dmel"

read1 <- paste(readDir, file_id, "_R1_001.fastq.gz", sep="")
read2 <- paste(readDir, file_id, "_R2_001.fastq.gz", sep="")

align(index=indexLoc, readfile1=read1, readfile2=read2, input_format="gzFASTQ", output_format="BAM", output_file=paste(file_id, ".Rsubread.bam", sep=""), nthreads=1)
EOF

echo "read finished aligning"

# Add read groups
module load samtools/1.9
# module load htslib/1.9


# if SNP VCFs not already created
# for chromosome in ${chromosomes[@]}; do
#     head -n 100 /scratch/caw5cv/IndoorData/data_raw/dgrp2.vcf | grep "^#" > /scratch/caw5cv/IndoorData/data_raw/dgrp2.${chromosome}.SNPs.vcf
#     grep "^${chromosome}" /scratch/caw5cv/IndoorData/data_raw/dgrp2.vcf | grep "SNP" >> /scratch/caw5cv/IndoorData/data_raw/dgrp2.${chromosome}.SNPs.vcf
# done

if [[ ${seqType}=="DNA" ]]; then
    java -jar ${picard} AddOrReplaceReadGroups \
           I=${fileID}.Rsubread.bam \
           O=/dev/stdout \
           RGID=${fileID} \
           RGLB="indoorCageData" \
           RGPU="indoorCageData" \
           RGPL=illumina \
           RGSM=${fileID} \
           QUIET=true | samtools sort - | \
           samtools rmdup - - > ${fileID}.bam && \
           rm ${fileID}.Rsubread.bam

elif [[ ${seqType}=="RNA" ]]; then
    outSuffix="bam"
    outFile=${fileID}.${outSuffix}
    java -jar ${picard} AddOrReplaceReadGroups \
           I=${fileID}.Rsubread.bam \
           O=/dev/stdout \
           RGID=${fileID} \
           RGLB="indoorCageData" \
           RGPU="indoorCageData" \
           RGPL=illumina \
           RGSM=${fileID} \
           QUIET=true | samtools sort - > ${fileID}.bam && \
           rm ${fileID}.Rsubread.bam
fi

samtools index  ${fileID}.bam
if [[ ! -f ${fileID}.bam.bai ]]; then
    exit 1
fi


for chromosome in ${chromosomes[@]}; do
java \
-Xmx6G \
-jar $gatk \
-R ${referenceGenome} \
-T ASEReadCounter \
-I ${fileID}.bam \
-o ${fileID}.${chromosome}.readcounts \
-sites /scratch/caw5cv/IndoorData/data_raw/dgrp2.${chromosome}.SNPs.vcf
done


# calculate some bam file statistics
samtools idxstats ${fileID}.bam > ${fileID}.stats

for chromosome in ${chromosomes[@]}; do

    echo working on chr $chromosome for individual ${fileID}
    length=$(awk -v chr=$chromosome '$1==chr {print $2}' ${fileID}.stats)

    # Run harp like
    echo running harp like
    $harp like \
    --bam $fileID.bam \
    --region $chromosome:1-$length \
    --refseq $referenceGenome \
    --snps $topDirectory/data_raw/$chromosome.priors.csv \
    --stem $fileID.$chromosome

    # Run harp freq
    echo running harp freq
    $harp freq \
    --bam $fileID.bam \
    --region $chromosome:1-$length \
    --refseq $referenceGenome \
    --snps $topDirectory/data_raw/$chromosome.priors.csv \
    --stem $fileID.$chromosome \
    --window_step $window_step \
    --window_width $window_width \
    --em_min_freq_cutoff 0.0001

    # clean up extra files and directories
    rm -rf ${fileID}.${chromosome}.output
    rm ${fileID}.${chromosome}.hlk
done

rm ${fileID}.stats

# Calculate most likely parents


# BEGIN calculate most likely parents

module load R/3.3.0

Rscript - <<EOF > $ind_id.mlp
#!/usr/bin/env Rscript

library(data.table)

harp.freqs <- fread("${fileID}.${chromosome}.freqs", header=FALSE, showProgress=FALSE, na.strings="-nan")
founders.list <- fread('founders.txt', header=FALSE)   # This file may need to be created
setnames(harp.freqs, c("chromosome","start","stop", founders.list[,V1]))

harp.freqs.long <- melt(harp.freqs, measure.vars = colnames(harp.freqs)[4:length(colnames(harp.freqs))], variable="lineID", value="freq")
harp.freqs.long[, q99 := quantile(freq, 0.99, na.rm=TRUE), by=chromosome]

mlp <- harp.freqs.long[, .N, by=list(chromosome, lineID, freq>=q99)][freq==TRUE][order(chromosome,-N)][,c("chromosome","lineID","N")]
write.table(mlp, file="", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# output: $ind_id.mlp
EOF

# archive readcounts after HARP
zip -Tm ${fileID}.zip ${fileID}.*.readcounts ${fileID}.*.freqs ${fileID}.Rsubread.bam.indel.vcf
