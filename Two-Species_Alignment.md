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

## Analysis in R using `GenomicBreaks`
Description

Load in these packages:
```r
library('BiocManager')
library('plyranges')
library('remotes')
remotes::install_github("oist/GenomicBreaks", repos=BiocManager::repositories())
library('GenomicBreaks') |> suppressPackageStartupMessages()
```

Load the ...maf.gz file:
```r
data <- load_genomic_breaks("~/localfolder/yourfile.maf.gz")
```

Synteny
```r
synteny_index(data)
synteny_index(swap(data))
```

Correlation
```r
correlation_index(data)
correlation_index(swap(data))
```

Gene order conservation
```r
GOC(data)
GOC(swap(data))
```

Strand randomisation index
```r
strand_randomisation_index(data)
```

Coalescing alignments
```r
coa <- coalesce_contigs(data)
length(data)
length(coa)
```

Genome Plots
```r
plotApairOfChrs(data, main = "Lokiarchaetoa / Odinarchaeota")

data |> forceSeqLengths() |> reverse(query = TRUE) |>
  plotApairOfChrs(main = "Lokiarchaetoa / Odinarchaeota (rev-complemented)")
  
plotApairOfChrs(coa, main = "Lokiarchaetoa / Odinarchaeota")

makeOxfordPlots(data, col = "strand") +
  scale_x_continuous() + scale_y_continuous() +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle("Lokiarchaetoa / Odinarchaeota `Oxford` plot")
```

Calculate Percent Identity:
```r
x <- data
width(x)
width(x$query)
x$matches / width(x)
sum(x$matches) / sum(width(x))
```



