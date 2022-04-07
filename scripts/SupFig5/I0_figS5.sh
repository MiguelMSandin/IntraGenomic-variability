#!/bin/bash

mkdir data/illumina/lulu
cd data/illumina/lulu

cp ../dada2/default_pipeline/out/ASVs_taxo.fasta .

vsearch \
  --threads 4
  --usearch_global ASVs_taxo.fasta \
  --db ASVs_taxo.fasta \
  --self \
  --id .84 \
  --iddef 1 \
  --userout match_list.txt \
  --userfields query+target+id \
  --maxaccepts 0 \
  --query_cov .9 \
  --maxhits 10
