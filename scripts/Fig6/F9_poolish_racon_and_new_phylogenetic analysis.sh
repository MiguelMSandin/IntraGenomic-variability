#!/bin/bash

# Poolishing ---------------------------------------------------------------------------------------
minimap2 "data/minion/consensus/demultiplexed_polycystines_consensus.fasta" "data/minion/raw/raw.fastq" > "data/minion/consensus/minimap_racon.paf"

racon "data/minion/raw/raw.fastq" "data/minion/consensus/minimap_racon.paf" "data/minion/consensus/demultiplexed_polycystines_consensus.fasta" > "data/minion/consensus/demultiplexed_polycystines_consensus_polished.fasta"

# BLAST --------------------------------------------------------------------------------------------
OUT_FMT="6 qseqid sseqid sacc stitle sscinames staxids sskingdoms sblastnames pident slen length mismatch gapopen qstart qend sstart send evalue bitscore"

DB="../DB/NCBI/" # Database with the taxonomy files
DBnt="$DB/nt" # Specify the header of the database file 'nt'

THREADS="2"

FILE="data/minion/consensus/demultiplexed_polycystines_consensus_polished.fasta"
BLAST_TSV="data/minion/blast/demultiplexed_polycystines_consensus_polished_blast.tsv"

# Run blast
export BLASTDB=$DB
blastn -num_threads $THREADS -max_target_seqs 100 -evalue 1.00e-10 -query $FILE -out $BLAST_TSV -db $DBnt -outfmt "$OUT_FMT"

# Remove old consensus sequence and realign new poolished sequence ---------------------------------
sequenceSelect.py -f "data/phylo/reference_consensus.fasta" -o "tmp/tmp" -p "minion" -a r -v

# Remove positions composed of only gaps
trimal -noallgaps -in "tmp/tmp" -out "tmp/tmp1"

# Align new poolished sequences --------------------------------------------------------------------
mafft --thread 2 --add "data/minion/consensus/demultiplexed_polycystines_consensus_polished.fasta" "tmp/tmp1" > "data/phylo/reference_consensus_poolished.fasta"
rm -f "tmp/tmp" "tmp/tmp1"

# And reconstruct again the phylogeny --------------------------------------------------------------
cd data/phylo/

ALIGNMENT="reference_consensus_poolished.fasta"
TRIM_THRESHOLD="30"

# Trimming alignment
TRIMMED=${ALIGNMENT/.fasta/_trimed$TRIM_THRESHOLD.fasta}
trimal -in $ALIGNMENT -out $TRIMMED -gt 0.$TRIM_THRESHOLD

# Launching the phylogenetic inference
BS="1000"
THREADS="2"

# In order to replicate a bit the analysis, we perform one inference under the model GTR+CAT
OUTPUT=${TRIMMED/.fasta/_cat_BS$BS}
raxmlHPC-PTHREADS-SSE3 -T $THREADS -m GTRCAT -c 25 -p $RANDOM -x $(date +%s) -d -f a -N $BS -n $OUTPUT -s $TRIMMED

# And another inference under the model GTR+GAMMA
OUTPUT=${TRIMMED/.fasta/_gamma_BS$BS}
raxmlHPC-PTHREADS-SSE3 -T $THREADS -m GTRGAMMA -p $RANDOM -x $(date +%s) -d -f a -N $BS -n $OUTPUT -s $TRIMMED
