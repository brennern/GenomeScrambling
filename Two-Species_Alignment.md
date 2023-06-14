# Two-Species Alignment
## Data Mining
Choose two species to gather genome data from using the NCBI Genome database (https://www.ncbi.nlm.nih.gov/datasets/genome/). 
Then, for each species, view the 'legacy Assembly page', click the 'FTP directory for RefSeq/GenBank assembly' link, and copy the ...genomic.fna.gz file link.

## Producing the Alignment
Load both the bioinfo-ugrp-modules and Nextflow2 modules.

```shell
ml bioinfo-ugrp-modules
ml Nextflow2
```

Create the following command using the ...genomic.fna.gz files from both species you want to compare.

```shell
nextflow run oist/plessy_pairwiseGenomeComparison -r main \
  --target https://ftp.ncbi.nlm.nih.gov/SPECIES1_genomic.fna.gz 
  --query https://ftp.ncbi.nlm.nih.gov/SPECIES2_genomic.fna.gz 
  -profile oist 
  -w /path/to/your/workspace
```

The results will be stored in `results/last`. To visualize plots, use XQuartz and the command eog. The ...postmasked.maf.gz file may be of interest for your futher analysis.
```shell
cd results/last
eog query.08.plot.png
```


