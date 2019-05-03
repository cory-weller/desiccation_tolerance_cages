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
workDir="${topDirectory}/Proj_Reconstruction/${reconstructionGroup}/"
bedFile="${topDirectory}/data_external/recombination.bed"

ind_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${parameterFile})

mkdir -p $workDir && cd $workDir

#RABBIT Parameters
eps="0.005"
epsF="0.005"
RABBITmodel="jointModel"
RABBITestfun="origViterbiDecoding"
RABBITpackageLocation="${topDirectory}/etc/RABBIT/RABBIT_Packages/"
topNfounders=16     #topNfounders indicates that X most likely founders are chosen for performing higher-resolution imputation
maxNSites=5000      # Maximum number of SNPs per chromosome to use for imputation


ind_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${seqIDfile})


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
