#By Dr. Michael Mansfield and Dr. Charles Plessy at OIST

cut -f 1,2 ../Enterobacteriaceae.tsv |
  sed 1d |
  while read ID FILE
  do
    nextflow run -r v5.2.0 oist/plessy_pairwiseGenomeComparison \
      --input ../Enterobacteriaceae.tsv \
      -profile oist \
      -w /flash/LuscombeU/deletemeCharlesPlessy/nf_work_Enterobacteriaceae \
      --target $FILE \
      -resume
    mv results results_$ID
  done
