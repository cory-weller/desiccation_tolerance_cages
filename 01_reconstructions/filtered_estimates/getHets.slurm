#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 10
#SBATCH --mem 50G
#SBATCH --time 0-2:00:00
#SBATCH --partition standard
#SBATCH --account berglandlab


function getHet {
  header=($(grep -m 1 "^#CHROM" filtered.all.vcf))
  i=${1}
  sampID=${header[$i]}
  cut -f 1,2,$i filtered.all.vcf | grep -e "0/1" -e "1/0" > $sampID.het
}

export -f getHet

parallel -j 10 getHet ::: $(seq 11 680)
