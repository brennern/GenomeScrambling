#!/bin/sh

# We use the -resume option to avoid re-downloading all the files
# But we purge the work directory from all the non-staging files
# to save some space.

cut -f 1,2 $1.tsv  |
    sed 1d |
    while read ID FILE
    do
      printf "$ID\n"
      nextflow run oist/plessy_pairwiseGenomeComparison \
	  -r lastal--split \
          --input $1.tsv \
          -profile oist \
          -w /flash/LuscombeU/$(whoami)/nf_work_$1 \
          --target $FILE \
	  --targetName $ID \
          -resume \
          --skip_dotplot_1 \
          --skip_dotplot_2 \
          --skip_dotplot_3 \
          --skip_m2m
      rm -rf /flash/LuscombeU/$(whoami)/nf_work_$1/??
      mv results results_$ID
  done
