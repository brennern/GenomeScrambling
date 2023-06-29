# Extracting Percent IDs from Dendrogram Nodes
## Creating the Dendrogram
First, load the following libraries.
```r
library("GenomicBreaks") |> suppressPackageStartupMessages()
library("plotly")
library("dplyr")
library("dendextend")
library("BiocManager")
library("ggtree")
library("ape")
library("phytools")
library("tidytree")
library("ggplot2")
library("TreeTools")
```

Next, load in the data just like in the AllVsAllComparison.md workflow.
```r
params <- list()
params$resultsDir <- '/flash/LuscombeU/noa/halobacteria/results/genomicbreaks'
yamlFiles <- list.files(params$resultsDir, pattern = "*.yaml", full.names = TRUE)
names(yamlFiles) <- yamlFiles |> basename() |> sub(pat = ".yaml", rep="")

getStats <- function(file) {
  y <- yaml::read_yaml(file) |> yaml::yaml.load()
  unlist(y)
}

df <- do.call(rbind, lapply(yamlFiles, getStats)) |> as.data.frame()
df <- df[,colSums(df, na.rm = TRUE) !=0]
df$species1 <- strsplit(rownames(df), "___") |> lapply(\(.) .[1]) |> unlist()
df$species2 <- strsplit(rownames(df), "___") |> lapply(\(.) .[2]) |> unlist()
df <- df[df$species1 != df$species2,]
df$percent_identity_global   <- df$matches_number_Total    / df$aligned_length_Total * 100

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

Make the dendrogram.
```r
m <- makeMatrix(df, "percent_identity_global", 100, 50)
hclust <- hclust(dist(m), method = "complete")
dend <- as.dendrogram(hclust)
```

## 'extractPercentID' Function
We will create a function that iterates through the dendrogram, calculating the percent identity between every internal node
with the consideration of their "children".
```r
extractPercentID <- function(node, tibble, matrix, fun = mean) {
  children <- child(tibble, node)
  stopifnot(nrow(children) == 2)
  left_side_node  <- children[1, "node", drop = TRUE]
  right_side_node <- children[2, "node", drop = TRUE]
  tipLabels <- function(node, tibble) {
    tipLabel <- tibble[node,"label", drop = TRUE]
    if(!is.na(tipLabel)) return(tipLabel)
    offspring(tibble, node)$label |> Filter(complete.cases, x=_)
  }
  left_side_species  <- tipLabels(left_side_node,  tibble)
  right_side_species <- tipLabels(right_side_node, tibble)
  comparison <- matrix[left_side_species, right_side_species, drop=F]
  percentID <- fun(comparison)
  percentID
}
```

Then, we will apply this function to the dendrogram.
```r
as_tibble(as.phylo(dend)) -> td
pIDs <- unique(td$parent) |> sort() |> purrr::set_names() |> sapply(extractPercentID, td, m)
td$percentID <- NA
td[names(pIDs), "percentID"] <- unname(pIDs)
```

Finally, utilize 'ggtree' to plot the percent identities on the dendrogram.
```r
ggtree(as.treedata(td), branch.length='none') + 
  geom_tiplab(as_ylab=TRUE, color="purple", size = 7, ) + 
  geom_label(aes(label=round(percentID, digit = 2), color = percentID), label.size = 0.25, size = 2.5, na.rm = TRUE) +
  scale_color_viridis_c()
```
