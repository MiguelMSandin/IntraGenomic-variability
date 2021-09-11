#!/bin/bash

DB="../DB/PR2/pr2_version_4.14.0_SSU_taxo_long_V4.fasta"
FASTA="data/illumina/dada2/default_pipeline/out/ASVs_taxo.fasta"
OUT=${FASTA/.fasta/_id.tsv}

vsearch --usearch_global $FASTA --db $DB --id 0 --blast6out $OUT --log ${FASTA/.fasta/_id.log}
