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
library('ggplot2')
library('pheatmap')
library('plotly')
```

For the analysis, we will only want to include the `.yaml` files. Load the results and filter for `.yaml` files.
```r
params <- list()
params$resultsDir <- '/your/working/directory/results/genomicbreaks'
yamlFiles <- list.files(params$resultsDir, pattern = "*.yaml", full.names = TRUE)
names(yamlFiles) <- yamlFiles |> basename() |> sub(pat = ".yaml", rep="")
```

Create this getStats function which will be used in the next step.
```r
getStats <- function(file) {
  y <- yaml::read_yaml(file) |> yaml::yaml.load()
  unlist(y)
}
```

Additional formatting:
```r
df <- do.call(rbind, lapply(yamlFiles, getStats)) |> as.data.frame()
df <- df[,colSums(df, na.rm = TRUE) !=0]
df$species1 <- strsplit(rownames(df), "___") |> lapply(\(.) .[1]) |> unlist()
df$species2 <- strsplit(rownames(df), "___") |> lapply(\(.) .[2]) |> unlist()
df <- df[df$species1 != df$species2,]
```

Create this makeMatrix function.
```r
makeMatrix <- function(df, column, defaultDiagonal = 100, defaultValue = NA) {
  species <- unique(df$species2)
  m <- matrix(defaultValue, nrow=length(species), ncol=length(species))
  colnames(m) <- rownames(m) <- species
  for (i in 1:length(species)) {
    m[i,i] <- defaultDiagonal
  }
  for (i in 1:nrow(df)) {
    s1 <- df[i, "species1"]
    s2 <- df[i, "species2"]
    if(s1 %in% species)
      m[s1, s2] <- df[i, column]
  }
  m
}
```

Analyze percent identities and mismatches:
```r
df$percent_identity_global <- df$matches_number_Total / df$aligned_length_Total * 100
df$percent_mismatches_global <- df$mismatches_number_Total / df$aligned_length_Total * 100
ggplot(df) + geom_point() + aes(percent_identity_global, percent_mismatches_global)
```

Percent identity global heatmaps:
```r
m <- makeMatrix(df, "percent_identity_global")
pheatmap::pheatmap(as.matrix(cluster::daisy(m)))
m <- makeMatrix(df, "percent_identity_global", 100, 50)
pheatmap::pheatmap(as.matrix(m), sym=T)
```

Fraction of genome aligned:
```r
df$fraction_genome_aligned_target <- df$aligned_target_Total / df$guessed_target_length * 100
df$fraction_genome_aligned_query  <- df$aligned_query_Total  / df$guessed_query_length  * 100
df$fraction_genome_aligned_avg    <- (df$fraction_genome_aligned_target + df$fraction_genome_aligned_query) / 2

ggplot(df) + geom_point() + aes(percent_mismatches_global, fraction_genome_aligned_avg, col = percent_identity_global)
ggplot(df) + geom_point() + aes(percent_identity_global,   fraction_genome_aligned_avg, col = percent_mismatches_global)
```

Fraction of genome chained:
```r
df$fraction_genome_chained_target <- df$chain_target_Total / df$guessed_target_length * 100
df$fraction_genome_chained_query  <- df$chain_query_Total  / df$guessed_query_length  * 100
df$fraction_genome_chained_avg    <- (df$fraction_genome_chained_target + df$fraction_genome_chained_query) / 2

ggplot(df) + geom_point() + aes(percent_mismatches_global, fraction_genome_chained_avg, col = percent_identity_global)
```

Width of aligned regions:
```r
##Aligned Width 1
ggplot(df) + geom_point() + aes(percent_identity_global, aligned_length_Mean, col = percent_mismatches_global)
ggplot(df) + geom_point() + aes(percent_mismatches_global, aligned_length_Mean, col = percent_identity_global)

##Aligned Width 2
ggplot(df) + geom_point() + aes(percent_identity_global, aligned_target_Mean / guessed_target_length, col = percent_mismatches_global) + scale_y_log10()
ggplot(df) + geom_point() + aes(percent_identity_global, aligned_query_Mean  / guessed_query_length, col = percent_mismatches_global) + scale_y_log10()
ggplot(df) + geom_point() + aes(percent_identity_global, aligned_length_Mean  / (guessed_target_length + guessed_query_length) / 2, col = percent_mismatches_global)  + scale_y_log10()
```

Length of the chains:
```r
ggplot(df) + geom_point() + aes(percent_identity_global, chain_target_Mean, col = percent_mismatches_global) + scale_y_log10()
ggplot(df) + geom_point() + aes(percent_mismatches_global, chain_target_Mean, col = percent_identity_global) + scale_y_log10()

ggplot(df) + geom_point() + aes(percent_identity_global, chain_target_Mean / guessed_target_length, col = percent_mismatches_global) + scale_y_log10()
ggplot(df) + geom_point() + aes(percent_identity_global, chain_query_Mean  / guessed_query_length, col = percent_mismatches_global) + scale_y_log10()
```

Calculate averages for synteny, correlations, strand randomisation indices, etc.
```r
df$index_avg_synteny      <- ( df$index_synteny_target + df$index_synteny_query ) / 2
df$index_avg_correlation  <- ( df$index_correlation_target + df$index_correlation_query ) / 2
df$index_avg_GOCvicinity4 <- ( df$index_GOCvicinity4_target + df$index_GOCvicinity4_query ) / 2
df$index_avg_strandRand   <- ( df$index_strandRand_target + df$index_strandRand_query ) / 2

df[,grepl("index_avg", colnames(df))] |> pairs()
```

Percent Identity Vs. Strand Randomisation Index Graph:
```r
df$lab <- ifelse(df$species1 > df$species2,
       paste(df$species1, df$species2, sep = "\n"),
       paste(df$species2, df$species1, sep = "\n"))

df |>
  #dplyr::filter(index_avg_synteny > 0.25) |>
  dplyr::group_by(lab) |>
  dplyr::summarise(percent_identity_global = mean(percent_identity_global),
                   index_avg_strandRand = mean(index_avg_strandRand),
                   lab = lab) |>
  ggplot() +
    geom_point() +
    aes(percent_identity_global, index_avg_strandRand, label=lab) +
    scale_x_continuous("Synteny index") +
    scale_y_continuous("Strand randomisation index") +
    theme_bw() +
    ggtitle("fix me later") +
    geom_text() -> gg
    
plotly::ggplotly(gg)
```

Similarity Vs. Synteny Index Graph:
```r
ggplot(df) + theme_bw() +
  aes(percent_identity_global, index_avg_synteny, col = percent_mismatches_global) +
  xlab("Identity between aligned regions (%)") +
  ylab("Synteny index") +
  geom_point()
```

Similarity Vs. Correlation Index Graph:
```r
ggplot(df) + theme_bw() +
  aes(percent_identity_global, index_avg_strandRand, col = percent_mismatches_global) +
  xlab("Identity between aligned regions (%)") +
  ylab("Strand randomisation index") +
  geom_point() + geom_smooth()
```

Gene Order Correlation Graph:
```r
ggplot(df) +
  aes(percent_identity_global, index_avg_GOCvicinity4) +
  aes(color = species1) +
  xlab("Similarity between aligned regions") +
  ylab("abs(correlation index)") +
  geom_point()
```


