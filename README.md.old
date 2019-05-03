# Starvation Tolerance

## Mapping RNA sequence data

I mapped RNA sequencing reads using [*iMapSplice*](https://github.com/LiuBioinfo/iMapSplice) against the dm3 reference genome.
The gene annotation file `dm3.genepred` contains gene prediction table of *D. melanogaster* FlyBase genes from UCSC table browser.

  * clade: Insect
  * genome: D.melanogaster
  * assembly: Apr. 2006 (BDGP R5/dm3)
  * group: Genes and Gene predictions
  * track: FlyBase Genes
  * table: FlyBaseGene
  * region: genome

Download `dm3.genepred` from pre-filled Table Browser link [here](https://genome.ucsc.edu/cgi-bin/hgTables?&clade=insect&org=D.+melanoqgaster&db=dm3&hgta_group=genes&hgta_track=flyBaseGene&hgta_table=flyBaseGene&hgta_regionType=genome&hgta_outputType=primaryTable&hgta_outFileName=dm3.genepred).

Perform a couple of fixes and  excluded the first column (bin) and header row:

```
# Remove bin and header row
zcat dm3.genepred.gz | cut -f 2- | tail -n +2 > dm3.genepred

# Generate gtf file from dm3.genepred:
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
chmod +x ./genePredToGtf
./genePredToGtf file dm3.genepred dm3.gtf

# Fix gene IDs to be identical across isoforms
sed -i 's/\(gene_id "CG[^-]*\)-R./\1/g' dm3.gtf
```

iMapSplice gene annotation file format requirements, from Xinan Liu of U Kentucky

| column | value |
---------|--------
| 1 | id |
| 2 | chromosome |
| 3 | strand |
| 4 | transcription start pos |
| 5 | transcription end pos |
| 6 | coding region start pos |
| 7 | coding region end pos |
| 8 | exon count |
| 9 | exon start positions |
| 10 | exon end positions |
| 11| additional information 1 |
| 12 | additional information 2 |

To create the iMapSplice required format, I padded the required additional information columns

```
sed 's/$/\tINFO1\tINFO2/g' dm3.genepred > dm3.imapsplice.gaf
```

And then removed any rows that didn't contain chromosomes I was mapping to

```
bash pruneGAF.sh dm3.imapsplice.gaf
```
Which generates `dm3.imapsplice.pruned.gaf`


To convert CG annotation IDs to Flybase IDs, I pulled the annotation IDs (removing isoform suffixes) from the first column of `dm3.genepred`:

```
<dm3.genepred cut -f 1 | cut -d "-" -f 1 | sort -u > dm3.genepred.AnnoIDs
```

`dm3.genepred` contents (21243 rows):
```
CG1674-RB       chr4    +       251355  266500  252579  266389  11      251355,252560,252904,254890,255489,257020,257894,260939,263891,264259,265805,   251521,252603,253474,254971,255570,257101,258185,261024,264211,264374,266500,
CG1674-RC       chr4    +       252055  266500  252579  266389  9       252055,252560,254890,257020,257894,260939,263891,264259,265805, 252253,252603,254971,257101,258185,261024,264211,264374,266500,
CG1674-RA       chr4    +       252108  265205  252579  264217  12      252108,252560,252904,254890,257020,257894,259747,260939,261321,263086,263891,265043,    252230,252603,253474,254971,257101,258185,259864,261024,261459,263305,264374,265205,
CG1710-RC       chr4    +       380577  395610  381316  393018  13      380577,381212,382929,383296,383571,384119,384493,386299,390132,390435,391716,392438,392875,     380601,381632,383239,383508,383760,384265,385567,387388,390370,390590,391983,392802,395610,
CG1710-RB       chr4    +       380587  395610  381316  393018  13      380587,381212,382929,383296,383571,384119,384493,386299,390132,390435,391716,392438,392875,     380606,381632,383239,383508,383760,384265,385567,387388,390370,390590,391983,392802,395610,
...
```

the resulting `dm3.genepred.AnnoIDs` contents (14058 rows):
```
CG00000
CG10000
CG10001
CG10002
CG10005
...
```

and used the [FlyBase batch converter](https://flybase.org/convert/id) to generate a conversion table `FlyBase_IDs.txt` for CG annotation IDs and corresponding FlyBase IDs. I renamed the conversion table `dm3.CG_to_FB.txt`

`dm3.CG_to_FB.txt` contents (14552 rows, which is less than `dm3.genepred.AnnoIDs` because some current IDs are mapped to by multiple submitted IDs):
```
# submitted_id  current_id      converted_id    current_symbol
CG4190  FBgn0001229     FBgn0001229     Hsp67Bc
CG2986  FBgn0015521     FBgn0015521     RpS21
CG18628 FBgn0036091     FBgn0036091     CG18628
CG11880 FBgn0039637     FBgn0039637     Ctl2
...
```


## Generating Feature Counts
 I generated feature counts using the subread aligner, implemented in the R package RSubRead. The subread aligner requires a gtf format file, which I made from the gene prediction table from UCSC via the `genePredToGft` command-line tool, available [here](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf).

```
# Retrieve tool if not already in directory
if [[ ! -f genePredToGtf ]]; then
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
  chmod 111 genePredToGtf
fi

# exclude first row and column and perform gtf conversion
tail -n +2 dm3.genepred | cut -f 2- > dm3 &&
./genePredToGtf file dm3 dm3.gtf
rm dm3
```

I then used `dm3.gtf` to generate feature counts (see `getFeatureCounts.Rscript`), which produces
`gene_lengths.txt`, a file containing Annotation ID and gene length
`d
