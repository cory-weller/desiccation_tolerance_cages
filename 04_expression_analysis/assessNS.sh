# Download snpEff using wget

# assess nonsynonymous - synonyous
wget -L http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip && unzip snpEff_latest_core.zip

cd snpEff/
mkdir dm3
cp ../dm3.gtf .


java -jar snpEff.jar build -gtf22 -v dm3


java -jar snpEff.jar build -gtf22 -v dm3

# annotate vcf file
java -Xmx4g -jar snpEff.jar -v dm3 ../../03_map_rna/SNPs/dm3.variants > dm3.snpeff.vcf

# remove superfluous info from annotated vcf
cut -d "|" -f 1-4 dm3.snpeff.vcf | less | tr "|" "\t" | sed 's/-R.$//g' > dm3.snpeff.edit.vcf

# in R:
