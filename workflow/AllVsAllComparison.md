# All Vs. All Genome Alignments
This workflow will allow you to collect all available genomes in a specific taxon, conduct an all vs. all comparison by aligning all genomes against each other, and finally analyze the data in R.

Most scripts we will use require the following modules:
```shell
ml bioinfo-ugrp-modules
ml Nextflow2
```

## Data Mining
To collect all available RefSeq genomes from a particular taxon, first enter the taxon name on the ___ website. If the taxon name is available for search, utilize the `makeNcbiTable.sh` script in the `scripts` folder and change each TAXON_NAME to your desired taxon name. For example, if the taxon name on the NCBI Genome Database is `Thermococci`, replace each TAXON_NAME in the script with `Thermococci`.

## All Vs. All Alignments
Once you have your `.tsv` file containing your genome data, run the `AllVsAllComparison.sh` script. Again, change the TAXON_NAME.tsv to match your `.tsv` file name. Also, make sure to change the output directory in the script.

## `GenomicBreaks` Input File
After running the `AllVsAllComparison.sh` script, you should have a multitude of `results_...` files in your results directory. The next script, named `makeGBreaksInputFile.sh` is going to format the data for the analysis in R. This will produce an input file in the `.tsv` format, upon which the script will prompt you to run a `nextflow` command containing that file. Here is the command:

```shell
nextflow run oist/plessy_nf_GenomicBreaks -profile oist -r main --input input.tsv --skel 'https://raw.githubusercontent.com/oist/GenomicBreaks/95cad1b661ff756f22e7e2794b79f0d4b48dc3fc/inst/rmarkdown/templates/countFeatures/skeleton/skeleton.Rmd' -w /flash/LuscombeU/yourname/thefolderyoulike\n"
```

## R Analysis
Now, you should have a results directory called `results/genomicbreaks` within your working directory. The `genomicbreaks` directory will contain all of the results, which include `.html` and `.yaml` files. 

First, load the following packages in R:
```r
library('GenomicBreaks') |> suppressPackageStartupMessages()
library("ggplot2")
library("pheatmap")
```
For the analysis, we will only want to include the `.yaml` files. Load the results and filter for `.yaml` files.
```r
params <- list()
params$resultsDir <- '/your/working/directory/results/genomicbreaks'
yamlFiles <- list.files(params$resultsDir, pattern = "*.yaml", full.names = TRUE)
names(yamlFiles) <- yamlFiles |> basename() |> sub(pat = ".yaml", rep="")
```



