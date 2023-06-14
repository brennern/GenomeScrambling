#!/bin/bash

# This snippet will not work unless NCBI Datasets and the NCBI EUtils are installed.

if ! command -v datasets >/dev/null 2>&1; then
  printf "datasets command not found\n"
  printf "See https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/\n"
  printf 'Or `ml bioinfo-ugrp-modules Other/ncbi-datasets-cli`\n'
  exit 1
fi

if ! command -v dataformat >/dev/null 2>&1; then
  printf "dataformat command not found\n"
  printf "See https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/\n"
  printf 'Or `ml bioinfo-ugrp-modules Other/ncbi-datasets-cli`\n'
  exit 1
fi

if ! command -v esearch >/dev/null 2>&1; then
  printf "esearch command not found\n"
  printf "See https://www.ncbi.nlm.nih.gov/books/NBK179288/\n"
  printf 'Or `ml bioinfo-ugrp-modules DebianMed ncbi-entrez-direct`\n'
  exit 1
fi

if ! command -v efetch >/dev/null 2>&1; then
  printf "efetch command not found\n"
  printf "See https://www.ncbi.nlm.nih.gov/books/NBK179288/\n"
  printf 'Or `ml bioinfo-ugrp-modules DebianMed ncbi-entrez-direct`\n'
  exit 1
fi

if ! command -v esummary >/dev/null 2>&1; then
  printf "esummary command not found\n"
  printf "See https://www.ncbi.nlm.nih.gov/books/NBK179288/\n"
  printf 'Or `ml bioinfo-ugrp-modules DebianMed ncbi-entrez-direct`\n'
  exit 1
fi

if ! command -v xtract >/dev/null 2>&1; then
  printf "xtract command not found\n"
  printf "See https://www.ncbi.nlm.nih.gov/books/NBK179288/\n"
  printf 'Or `ml bioinfo-ugrp-modules DebianMed ncbi-entrez-direct`\n'
  exit 1
fi

datasets summary genome taxon 'TAXON_NAME' \
	--assembly-source refseq \
	--as-json-lines \
	--assembly-level complete \
	--reference |
		dataformat tsv genome \
			--fields accession,assminfo-name,annotinfo-name,annotinfo-release-date,organism-name,organism-tax-id > ncbi_datasets.tsv

printf "id\tfile\tGenomeID\tBinomial\tAccNum\tTaxID\tPMID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tSubrank1\tSubrank2\n" > TAXON_NAME.tsv

cat ncbi_datasets.tsv | grep -v '^Assembly' | while IFS=$'\t' read ASSEMBLYACC ASSEMBLYNAME ANNOTATIONNAME ANNOTATIONRELEASE ORGNAME TAXID;
do
	# Split up the assembly name into chunks to stitch together the link to the NCBI FTP server programmatically.
	# So far, this works for the ~200 genomes here. It may not work for others.
	# Alternatively you can get this from an esearch against the NCBI Assembly database into an esummary.
	ACC1=$(echo "${ASSEMBLYACC}" | cut -d "_" -f 1) ;
	ACC2=$(echo "${ASSEMBLYACC}" | cut -d "_" -f 2 | cut -c 1,2,3) ;
	ACC3=$(echo "${ASSEMBLYACC}" | cut -d "_" -f 2 | cut -c 4,5,6) ;
	ACC4=$(echo "${ASSEMBLYACC}" | cut -d "_" -f 2 | cut -c 7,8,9) ;
	LINK=$(echo "https://ftp.ncbi.nlm.nih.gov/genomes/all/"${ACC1}"/"${ACC2}"/"${ACC3}"/"${ACC4}"/"${ASSEMBLYACC}"_"${ASSEMBLYNAME}"/"${ASSEMBLYACC}"_"${ASSEMBLYNAME}"_genomic.fna.gz");

	# Make a UNIX-safe name out of the NCBI Datasets "Organism name" column (otherwise "binomials" can contain e.g. ":", "/", "(", etc.)
	SAFEID=$(echo "${ORGNAME}" | sed 's/ /_/g'  | sed 's/[^[:alnum:]+._-]//g')

	# Retrieve PMIDs for assembly.
	# Note that this may be prohibitively large for some accessions. Depends on how NCBI has curated that record. I haven't tried too many of them.
	# It collapses all of the IDs and separates them by comma (the sed expression at the end).
	PMIDS=$(esearch -db assembly -query \""${ASSEMBLYACC}"\" </dev/null | elink -target pubmed | esummary | xtract -pattern DocumentSummary -element Id -sep "," -tab "," | sed -n 'H;${x;s/\n/,/g;s/^,//;p;}')

	# Retrieve taxonomic ranks for assembly.
	# Change -element ScientificName to -element Rank,ScientificName if you want the output to specify
	# what Subrank1 and Subrank2 are (because it could be either a subspecies or a serotype).
	TAXONOMY=$(efetch -db taxonomy -id \""${TAXID}"\" -format xml </dev/null | xtract -pattern Taxon -element Taxon -def "NA" -block "*/Taxon" -if Rank -equals "superkingdom" -or Rank -equals "phylum" -or Rank -equals "class" -or Rank -equals "order" -or Rank -equals "family" -or Rank -equals "genus" -or Rank -equals "species" -or Rank -equals "subspecies" -or Rank -equals "serotype" -tab "____"  -sep ":" -element ScientificName)

	echo -e ""${SAFEID}"\t"${LINK}"\t"${ASSEMBLYNAME}"\t"${ORGNAME}"\t"${ASSEMBLYACC}"\t"${TAXID}"\t"${PMIDS}"\t"${TAXONOMY}"" | sed 's/____/\t/g' >> TAXON_NAME.tsv
	sleep 5s
done 
