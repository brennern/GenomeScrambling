# GenomeScrambling
Repository for workflows and data utilized in the Luscombe Unit's 'Scrambling in the Tree of Life' project.

## Data Mining

Species genome data was collected via the usage of the makeNcbiTable.sh, AllVsAllComparison.sh (Nextflow pairwise alignment pipeline), and makeGBreaksInputFile.sh (GenomicBreaks pipeline) scripts.

``` shell
ml bioinfo-ugrp-modules
ml Nextflow2

./makeNcbiTable.sh
./AllVsAllComparison.sh
./makeGBreaksInputFile.sh

nextflow run oist/plessy_nf_GenomicBreaks -profile oist -r main --input input.tsv --skel 'https://raw.githubusercontent.com/oist/GenomicBreaks/95cad1b661ff756f22e7e2794b79f0d4b48dc3fc/inst/rmarkdown/templates/countFeatures/skeleton/skeleton.Rmd' -w /flash/LuscombeU/yourname/thefolderyoulike\n"
```

## R Analysis

Upon completion of the All Vs. All genome comparisons, the resulting .yaml files were analyzed in R with the ScrambledTreeBuilder Package. 

``` r
library(ScrambledTreeBuilder)

resultsDir <- system.file("extdata/PairwiseComparisons", package = "ScrambledTreeBuilder")
yamlFileData <- list.files(resultsDir, pattern = "*.yaml.bz2", full.names = TRUE)
names(yamlFileData) <- yamlFileData |> basename() |> sub(pat = ".yaml.bz2", rep="")

exDataFrame <- formatStats(yamlFileData)

valuesToBuildTheTree <- "percent_identity_global"
treeMatrix <- makeMatrix(exDataFrame, valuesToBuildTheTree, 100, 50)
valuesToPlaceOnLabels <- "index_avg_strandRand"
valueMatrix <- makeMatrix(exDataFrame, valuesToPlaceOnLabels, 1, 0.5)

HClust <- hclust(dist(treeMatrix), method = "complete")
Tibble <- tidytree::as_tibble(tidytree::as.phylo(HClust))
tibbleWithValue <- makeValueTibble(Tibble, valueMatrix)

Tree <- visualizeTree(tibbleWithValue, tibbleWithValue$value)

Tree + 
  ggplot2::ggtitle(paste("Tree built with", valuesToBuildTheTree, "and labelled with", valuesToPlaceOnLabels)) + 
  viridis::scale_color_viridis(name = valuesToPlaceOnLabels)
```
