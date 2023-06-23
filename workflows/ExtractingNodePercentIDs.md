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

## 'extractPercentID' and 'addLabels' Functions
We will create a function that iterates through the dendrogram, calculating the percent identity between every internal node
with the consideration of their "children".
```r
extractPercentID <- function(tree, matrix, ...) {
  left_side <- tree[[1]]
  right_side <- tree[[2]]
  left_side_species <- labels(left_side)
  right_side_species <- labels(right_side)
  comparison <- matrix[left_side_species, right_side_species, drop=F]
  percent_ID <- mean(comparison)
  percent_ID
}
```

Then, we will create a function which creates a percent identity attribute following each calculation, providing a label to each node.
```r
addLabels <- function(tree, fun, matrix) {
  if(is.null(attr(tree, "leaf"))) {
    attr(tree, "percentID") <- fun(tree, matrix)
  }
  tree
}
```

Finally, utilize 'dendrapply' to iterate both functions to the dendrogram object.
```r
dendrapply(dend, addLabels, extractPercentID, m)
```




