#!/bin/sh

#By Dr. Charles Plessy at OIST

if [ $# -eq 0 ]
then
  printf "Usage: $0 path/to/dir/containing/results_dirs > input.tsv\n"
  exit 1
fi
printf 'id\tfile\n'
for d in $(ls -d $(realpath "$1")/results_*)
do
  DIR=$(basename $d)
  SPEC1=${DIR##results_}
  for FILE in $(ls $d/last/*07.postmasked.maf.gz)
  do
    SPEC2=$(basename $FILE .07.postmasked.maf.gz)
    printf "${SPEC1}___${SPEC2}\t$FILE\n"
  done
done
printf "Use this input file with the nf-GenomicBreaks pipeline\n" > /dev/stderr
printf "Also do not forget to adjust the -w option to your case\n" > /dev/stderr
printf "nextflow run oist/plessy_nf_GenomicBreaks -profile oist -r main --input input.tsv --skel 'https://raw.githubusercontent.com/oist/GenomicBreaks/95cad1b661ff756f22e7e2794b79f0d4b48dc3fc/inst/rmarkdown/templates/countFeatures/skeleton/skeleton.Rmd' -w /flash/LuscombeU/yourname/thefolderyoulike\n" > /dev/stderr
