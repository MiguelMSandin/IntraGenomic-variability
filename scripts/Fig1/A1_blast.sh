#!/bin/bash

# General parameters for BLASTn and Output format
OUT_FMT="6 qseqid sseqid sacc stitle sscinames staxids sskingdoms sblastnames pident slen length mismatch gapopen qstart qend sstart send evalue bitscore"

DB="../DB/NCBI/" # Database with the taxonomy files
DBnt="$DB/nt" # Specify the header of the database file 'nt'

THREADS="2"


# For Sanger sequencing data -----------------------------------------------------------------------

# Create blast directory
[ ! -d "data/sanger/blast/" ] && mkdir -p "data/sanger/blast/"

# And set parameters
FILE="data/sanger/raw/raw.fasta"
BLAST_TSV="data/sanger/blast/raw_blast.tsv"

# Run blast
export BLASTDB=$DB
blastn -num_threads $THREADS -max_target_seqs 100 -evalue 1.00e-10 -query $FILE -out $BLAST_TSV -db $DBnt -outfmt "$OUT_FMT"


# For Sanger sequencing data concatenated ----------------------------------------------------------

# Create blast directory
[ ! -d "data/sanger/blast/" ] && mkdir -p "data/sanger/blast/"

# And set parameters
FILE="data/sanger/raw/raw_concatenated.fasta"
BLAST_TSV="data/sanger/blast/raw_concatenated_blast.tsv"

# Run blast
export BLASTDB=$DB
blastn -num_threads $THREADS -max_target_seqs 100 -evalue 1.00e-10 -query $FILE -out $BLAST_TSV -db $DBnt -outfmt "$OUT_FMT"


# For ONT data -------------------------------------------------------------------------------------

# Define variables and parameters
FILE="data/minion/raw/raw.fasta"
BLAST_TSV="data/minion/blast/raw_blast.tsv"

# Create blast 
[ ! -d "data/minion/blast/" ] && mkdir -p "data/minion/blast/"

# Run blast
export BLASTDB=$DB
blastn -num_threads $THREADS -max_target_seqs 100 -evalue 1.00e-10 -query $FILE -out $BLAST_TSV -db $DBnt -outfmt "$OUT_FMT"
