#!/bin/bash

# General parameters for BLASTn and Output format
OUT_FMT="6 qseqid sseqid sacc stitle sscinames staxids sskingdoms sblastnames pident slen length mismatch gapopen qstart qend sstart send evalue bitscore"

DB="../DB/NCBI/" # Database with the taxonomy files
DBnt="$DB/nt" # Specify the header of the database file 'nt'

THREADS="2"


# Sanger sequences ---------------------------------------------------------------------------------

# Define variables and parameters
FILE="data/minion/consensus/raw_concatenated_polycystines_consensus.fasta"
BLAST_TSV="data/minion/blast/raw_concatenated_polycystines_consensus_blast.tsv"

# Run blast
export BLASTDB=$DB
blastn -num_threads $THREADS -max_target_seqs 100 -evalue 1.00e-10 -query $FILE -out $BLAST_TSV -db $DBnt -outfmt "$OUT_FMT"


# MinION sequences ---------------------------------------------------------------------------------

# Define variables and parameters
FILE="data/minion/consensus/demultiplexed_polycystines_consensus.fasta"
BLAST_TSV="data/minion/blast/demultiplexed_polycystines_consensus_blast.tsv"

# Run blast
export BLASTDB=$DB
blastn -num_threads $THREADS -max_target_seqs 100 -evalue 1.00e-10 -query $FILE -out $BLAST_TSV -db $DBnt -outfmt "$OUT_FMT"
