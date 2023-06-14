# All Vs. All Genome Alignments
This workflow will allow you to collect all available genomes in a specific taxon, conduct an all vs. all comparison by aligning all genomes against each other, and finally analyze the data in R.

## Data Mining
To collect all available RefSeq genomes from a particular taxon, first enter the taxon name on the ___ website. If the taxon name is available for search, utilize the `makeNcbiTable.sh` script in the `scripts` folder and change each TAXON_NAME to your desired taxon name. For example, if the taxon name on the NCBI Genome Database is `Thermococci`, replace each TAXON_NAME in the script with `Thermococci`.

## All Vs. All Alignments
Once you have your `.tsv` file containing your genome data, run the `AllVsAllComparison.sh` script. Again, change the TAXON_NAME.tsv to match your `.tsv` file name. Also, make sure to change the output directory in the script.

## `GenomicBreaks` Input File
After running the `AllVsAllComparison.sh` script, you should have a multitude of `results_...` files in your results directory.
