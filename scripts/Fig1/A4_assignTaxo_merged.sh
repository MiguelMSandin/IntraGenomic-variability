#!/bin/bash

DB="../DB/PR2/pr2_version_4.14.0_SSU_taxo_long_V4.fasta"

for FILE in $(ls illumina/merged | grep \.fasta)
do
	echo "  Analyzing ${FILE}"
	vsearch --usearch_global "data/illumina/merged/${FILE}" --db $DB --id 0 --blast6out "data/illumina/merged/${FILE/.fasta/_id.tsv}" --log "data/illumina/merged/${FILE/.fasta/_id.log}"
	echo "----------------------------------------"
done

echo "Done"

