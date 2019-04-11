#!/usr/bin/env bash

# DESCRIPTION: Creates a separate .fa file for each contig/chromosome within the input
# EXAMPLE USAGE: bash splitFasta.sh reference.fa

inFasta=${1}

filestem=${inFasta%.*}

nMatches=$(grep -E "^>" ${inFasta} | awk -F "(^>|\t| )" '{print $2}' | grep ${filestem} | wc -l)
if [[ ${nMatches} > 0 ]]; then
    echo "warning: ${inFasta} contains a header with the first word ${filestem}."
    echo "warning: exiting to prevent over-writing ${inFasta}."
    echo "warning: rename reference to not be any of the first words in the header:"
    grep -E "^>" ${inFasta} | awk -F "(^>|\t| )" '{print $2}'
    echo ""
    exit 1
fi



# DETAILS: Functions by setting the field separator to a regular expression capturing ">", tab, and space.
# For each row in the file, if it is a header (begins with ">"), update the value of variable s
# to be equal to the first word in the header, and write a header line only containing that first word.
# Every line is written to a file <first_word_in_header>.fa, splitting the file into separate contigs.

awk -F "(^>|\t| )" '{if($0 ~ /^>/) {s=$2".fa";  print ">"$2 > s} else print > s}' ${inFasta}

