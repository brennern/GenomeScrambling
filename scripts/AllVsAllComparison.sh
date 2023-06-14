#By Dr. Michael Mansfield and Dr. Charles Plessy at OIST

cut -f 1,2 ../TAXON_NAME.tsv |
  sed 1d |
  while read ID FILE
  do
    nextflow run -r v5.2.0 oist/plessy_pairwiseGenomeComparison \
      --input ../TAXON_NAME.tsv \
      -profile oist \
      -w /your/directory/to/store/results \
      --target $FILE \
      -resume
    mv results results_$ID
  done
