#!/usr/bin/env bash
echo "missing iMapSplice mapped reads"
grep -n -v -f <(ls ./mapped_reads/iMapSplice/*.bam | cut -d "/" -f 4 | cut -d "." -f 1 | sed 's/$/$/g') filestems.txt


echo "missing RSubRead mapped reads"
grep -n -v -f <(ls ./mapped_reads/RSubRead/*.bam | cut -d "/" -f 4 | cut -d "." -f 1 | sed 's/$/$/g') filestems.txt
