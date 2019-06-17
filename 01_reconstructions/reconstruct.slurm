#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 12G
#SBATCH -t 0-2:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

# The parameter file is a tab-delimited file containing a list of individuals to reconstruct within a SLURM job array
topDirectory=$(pwd)
nGenerations=5
chromosomes=( "2L" "2R" "3L" "3R" "X" )
readsDir="${topDirectory}/dna_reads/"

# If file  with each sample ID doesn't exist, create it
if [ ! -f filestems.txt ]; then
  ls $readsDir | cut -d "_" -f 1 | sort -u > filestems.txt
fi

seqIDfile="${readsDir}/filestems.txt"
ind_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${seqIDfile})

workDir="${topDirectory}/tmp_reconstruct/${ind_id}/"
bedFile="${topDirectory}/etc/recombination.bed"

mkdir -p $workDir && cd $workDir


#RABBIT Parameters
eps="0.005"
epsF="0.005"
RABBITmodel="jointModel"
RABBITestfun="origViterbiDecoding"
RABBITpackageLocation="${topDirectory}/etc/RABBIT/RABBIT_Packages/"
topNfounders=16     #topNfounders indicates that X most likely founders are chosen for performing higher-resolution imputation
maxNSites=5000      # Maximum number of SNPs per chromosome to use for imputation


# BEGIN preparation of RABBIT input file

Rscript - <<EOF

#!/usr/bin/env Rscript

library(data.table)
library(foreach)


# in R, combine information from three sources:
# 1. Readcounts, to find out which SNPs is the whole set we're considering
# 2. population VCF file, for all variants
# 3. alleleFreqs, for information content
# 4. recombination .bed file

# 3. Of remaining SNPs, merge in read depth and allele frequency
# 4. Subset by top N, according to criteria deemed important
# 5. Re-organize by chromosome + POS
# 6. Iterate over columns and print to file (effectively transposing without holding in memory)


# 1. read in population VCF, readcounts, and freqs for this individual

ind_id <- "${ind_id}"
topNfounders <- ${topNfounders}
maxNSites <- ${maxNSites}


founders.list <- fread('founders.txt', header=FALSE, showProgress=FALSE)

vcf <- fread("zcat haplotypes.vcf.gz", skip="#CHROM", header=TRUE, na.strings="./.", showProgress=FALSE)
vcf <- vcf[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", founders.list[,V1]), with=FALSE]
setnames(vcf, "#CHROM", "CHROM")
setkey(vcf, "CHROM", "POS", "ID")

readcounts <- fread(paste("zcat ", ind_id, ".readcounts.gz", sep=""), header=TRUE, showProgress=FALSE)
setkey(readcounts, "contig", "position", "variantID")

mlp <- fread(paste(ind_id, ".mlp", sep=""), header=TRUE, showProgress=FALSE)
mlp[,rank := frank(-N, ties.method="random"), by=chromosome]
mlp.chosen <- mlp[rank <= topNfounders]


# Read in .bed file containing recombination rates
# Rquired column headers are chr (chromosome); start; stop; c (recombination rate, in units of cM/Mb)
bed <- fread("recombination.bed", header=TRUE, showProgress=FALSE)
# Key .bed file by chromosome and window start
setkey(bed, chr, start)

# Generate column for cM (centiMorgans)
bed[, cM := c*((stop-start)/1e6)]

# Calculate cumulative cM at the END of each window
bed[, cumulative_cM := cumsum(cM), by=chr]

# Generate functions (cM ~ BP) to translate base pair position (in VCF) to units of cM
recombination_function <- new.env()

# Create recombination function based on .bed file
for(chromosome in unique(bed[,chr])) {
    recombination_function[[as.character(chromosome)]] <- approxfun(c(0, bed[chr==chromosome][,stop]), c(0,bed[chr==chromosome][,cumulative_cM]))
}

#options(scipen=999)

for(chr.i in unique(bed[,chr])) {
    cat(chr.i)
    cat("\n")
    chosen.founders <- sort(mlp.chosen[chromosome==chr.i][,lineID])

    # If no most likely founders for this chromosome, then skip. This should not happen.
    if(length(chosen.founders)==0) {
        cat("No most likely founders for chromosome ")
        cat(chr.i)
        cat("\n")
        next
    }

    # if only one founder was chosen for the entire chromosome (e.g., homozygous for entire chromosome for a single founder)
    if(length(chosen.founders)==1) {
        writeLines(paste("Nonrecombinant Homozygous ", chosen.founders, sep=""), con=paste(ind_id, ".", chr.i, ".RABBIT.out.csv", sep=""))
        next
    }

    vcf.sub <- copy(vcf[CHROM==chr.i][,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", chosen.founders), with=FALSE])
    # if in 0/0 & 1/1 format instead of 0 or 1:
    if(any(vcf.sub=="0/0") | any(vcf.sub=="1/1")) {
        # Convert to factor with level 1 = "0/0", level 2 = "1/1"
        vcf.sub[, (chosen.founders) := lapply(.SD, factor, levels=c("0/0","1/1")), .SDcols=chosen.founders]
        # Convert to numeric
        vcf.sub[, (chosen.founders) := lapply(.SD, as.numeric), .SDcols=chosen.founders]
        # Subtract 1, such that "0/0" is now 0, "1/1" is now 1
        vcf.sub[, (chosen.founders) := lapply(.SD, "-", 1), .SDcols=chosen.founders]
    }

    vcf.sub[, nRef := apply(.SD, 1, function(x) sum(x == 0, na.rm=TRUE)), .SDcols=chosen.founders]
    vcf.sub[, nAlt := apply(.SD, 1, function(x) sum(x == 1, na.rm=TRUE)), .SDcols=chosen.founders]
    vcf.sub[, refFreq := nRef/(nRef+nAlt)]
    vcf.sub[, altFreq := nAlt/(nRef+nAlt)]
    setkey(vcf.sub, ID)

    # Merge in sequenced allele freq
    vcf.sub.merge <- merge(vcf.sub, readcounts, by.x="ID", by.y="variantID")[refFreq != 1 & altFreq != 1 & otherBases == 0 & improperPairs == 0]
    vcf.sub.merge[, c("totalCount","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs") := NULL]
    vcf.sub.merge[, "imputed_ind" := ifelse(refCount > 0 & altCount > 0, "12",
                               ifelse(altCount==0, "1N",
                               ifelse(refCount==0, "2N", "NN")))]

    # Select up to maximum number of SNPs
    vcf.sub.merge[,indx := 1:.N]

    # Retain ALL heterozygous sites
    vcf.sub.merge[refCount>0 & altCount>0, marker := TRUE]
    vcf.sub.merge[, freq := min(refFreq, altFreq), by=indx]
    nSites <- dim(vcf.sub.merge[is.na(marker)])[1]
    samplePercent <- 0.1
    sampleNSites <- trunc(nSites * samplePercent)
    retainmentPercent <- 0.002
    retainNSites <- trunc(nSites * retainmentPercent)



    if(dim(vcf.sub.merge)[1] <= maxNSites) {
        vcf.sub.merge[, marker := TRUE]
    } else {
        while(dim(vcf.sub.merge[marker == TRUE])[1] < maxNSites ) {
            indicesToMark <- vcf.sub.merge[is.na(marker)][sample(.N, size=sampleNSites, replace=TRUE)][order(freq)][1:retainNSites][,indx]
            vcf.sub.merge[indicesToMark, marker := TRUE]
        }
    }

    markers <- vcf.sub.merge[marker==TRUE][sample(.N, size=min(maxNSites, .N), replace=FALSE)]
    setkey(markers, "CHROM", "POS")
    markers[, c("contig", "position", "refAllele", "altAllele", "refCount", "altCount", "freq", "indx", "marker", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "refFreq", "altFreq", "nRef", "nAlt") := NULL]
    setnames(markers, "CHROM", "Chromosome")
    setnames(markers, "ID", "SNP")
    markers[, cM := recombination_function[[chr.i]](POS)][]
    markers[, POS := NULL]

    markers[, (chosen.founders) := lapply(.SD, "+", 1), .SDcols=chosen.founders]
    line.names=colnames(markers)[3:(length(colnames(markers))-1)]
    nFoundersUsed <- length(colnames(markers))-4

    setcolorder(markers, c("SNP","Chromosome","cM", line.names))
    writeLines(paste("#founders,",nFoundersUsed, sep=""), con=paste(ind_id, ".", chr.i, ".RABBIT.in", sep=""))
    write.table(t(as.matrix(markers)), file=paste(ind_id, ".", chr.i, ".RABBIT.in", sep=""), quote=FALSE, row.names=TRUE, col.names=FALSE, sep=",", na="N", append=TRUE)
}
EOF




for chromosome in ${chromosomes[@]}; do
{

python - <<EOF > $ind_id.$chromosome.RABBIT.m

print """SetDirectory["%s"]""" % "${RABBITpackageLocation}"
print """Needs["MagicReconstruct\`"]"""
print """SetDirectory["%s"]""" % "$(pwd)"
print """popScheme = Table["RM1-E", {%s}]""" % "${nGenerations}"
print 'epsF = %s' % "${epsF}"
print 'eps = %s' % "${eps}"
print 'model = "%s"' % "${RABBITmodel}"
print 'estfun = "%s"' % "${RABBITestfun}"
print 'inputfile = "%s"' % "$(pwd)/${ind_id}.${chromosome}.RABBIT.in"
print 'resultFile = "%s.txt"' % "$(pwd)/${ind_id}.${chromosome}.RABBIT.out"
print """magicReconstruct[inputfile, model, epsF, eps, popScheme, resultFile, HMMMethod -> estfun, PrintTimeElapsed -> True]"""
print 'summaryFile = StringDrop[resultFile, -4] <> ".csv"'
print 'saveAsSummaryMR[resultFile, summaryFile]'
print 'Exit'
EOF



module load mathematica

math -noprompt -script $ind_id.$chromosome.RABBIT.m

module remove mathematica

# END RABBIT
#
#
#
#
#
# BEGIN conversion of RABBIT output to haplotype map

echo -e "chromosome\tstart\tstop\tpar1\tpar2" > "${ind_id}.${chromosome}.estimate.haps"

filename=${ind_id}.${chromosome}.RABBIT.out.csv

nLines=$(wc -l $filename | awk '{print $1}')
if [[ ${nLines} -eq 1 ]]; then
    par=$(awk '{print $3}' $filename)
    chromosome=$(echo $filename | cut -d "." -f 2)
    echo -e "${chromosome}\t1\t25000000\t${par}\t${par}" >> "${ind_id}.${chromosome}.estimate.haps"
else
python - ${filename} <<EOF >> "${ind_id}.${chromosome}.estimate.haps"
#!/usr/bin/env python

import re
import sys

filename = sys.argv[1]
ind_id, chromosome = filename.split(".")[0:2]

def getLines(filename):
    mylines = []
    with open(filename, 'r') as infile:
        for line in infile:
            if not line.startswith(('Chromosome', 'cM')):
                mylines.append(line.rstrip())
                if 'ViterbiPath' in line:
                    viterbiLine = len(mylines)
        return [mylines, viterbiLine]

def getDiplotypes(mylines):
    diplotypes = {}
    pattern = re.compile('^diplotype[0-9]+')
    for line in mylines:
        if re.match(pattern, line):
            diplotype, parentCodes, founders = line.split(', ')
            diplotype = diplotype.split('diplotype')[-1]
            founders = founders.split('---')
            diplotypes[int(diplotype)] = founders
    return diplotypes

def getNucleotidePositions(mylines):
    positions = []
    for line in mylines:
        if line.startswith('SNP'):
            #return [x.split('_')[0:2] for x in line.rstrip().split(', ')] if returning Chromosome and position
            return [int(x.split('_')[1]) for x in line.rstrip().split(', ')[1:]] #exclude 1st item in line because it's simply 'SNP'

def getSampleDiplotypes(mylines, viterbiLine):
    paths = {}
    for i in range(viterbiLine, len(mylines)):
        sample, path = mylines[i].rstrip().split(', ')
        path = path.split('-')
        paths[sample] = path
    return paths

def phasePaths(paths):
    phasedPaths = {}
    for path in paths:
        entry = paths[path]
        pairs = []
        for i in range(0,len(entry),2):
            pairs.append(entry[i:i+2])
        for pair in pairs[:-1]:
            pair[1] = [x for x in diplotypes[int(pair[1])]]
        phasedPaths[path] = pairs
    return phasedPaths

def convertToPositions(phasedPaths):
    convertedPaths = {}
    for path in phasedPaths:
        segments = []
        for segment in phasedPaths[path][:-1]:
            segments.append([positions[int(segment[0]) - 1], segment[1]])
        for i in range(0,len(segments)-1):
            segments[i][0] = [int(segments[i][0]), int(segments[i+1][0])-1]
        segments[-1][0] = [int(segments[-1][0]), chromosomeLength]
        convertedPaths[path] = segments
    return convertedPaths

mylines, viterbiLine = getLines(filename)

#viterbiPath = [int(x) for x in mylines[viterbiLine].split()[1].split("-")]

# Get dictionary of all possible diplotype pairs and their numeric code
diplotypes = getDiplotypes(mylines)

# Get nucleotide positions from SNP ID
positions = getNucleotidePositions(mylines)

# Get the viterbi path from RABBIT output
paths = getSampleDiplotypes(mylines, viterbiLine)

phasedPaths = phasePaths(paths)

chromosomeLength = int(positions[-1])

convertedPaths = convertToPositions(phasedPaths)

for i in convertedPaths:
    for j in convertedPaths[i]:
        print '\t'.join([chromosome] + [str(x) for x in j[0]]) + '\t' + '\t'.join([str(x) for x in j[1]])

EOF
fi

# END conversion of RABBIT output to haplotype map
#
#
#
#
#
# BEGIN conversion of haplotype map to VCF file



# Set up associative array of lineID : column, to be used by TABIX
vals=($(zgrep -m 1 "^#CHROM" "${topDirectory}/data_raw/dgrp2.${chromosome}.SNPs.polarized.vcf.gz"))
declare -A founderIndices

for i in $(seq 0 "${#vals[@]}"); do
    let j=$i+1
    founderIndices[[${vals[$i]}]]=$j
done

# Read through diploid paths, extracting genotypes with TABIX, appending to estimate.vcf file
# Here is where you would change the bgzipped vcf filename to whatever one has all the sites you're trying to pull out based on the path
while read chr start stop par1 par2; do
    col1=${founderIndices[[$par1]]}
    col2=${founderIndices[[$par2]]}
    "${topDirectory}"/etc/htslib/bin/tabix "${topDirectory}/data_raw/dgrp2.${chr}.SNPs.polarized.vcf.gz" ${chr}:${start}-${stop} | awk -v col1=$col1 -v col2=$col2 '{print $1,$2,$col1,$col2}' >> ${ind_id}.${chromosome}.tmp.vcf
done < <(awk 'NR > 1 {print}' ${ind_id}.${chromosome}.estimate.haps)


awk '{print $1,$2,$3"/"$4}' ${ind_id}.${chromosome}.tmp.vcf > ${ind_id}.${chromosome}.estimate.vcf && \
rm ${ind_id}.${chromosome}.tmp.vcf
} || continue

done

zip -Tm "${ind_id}".zip "${ind_id}".*
