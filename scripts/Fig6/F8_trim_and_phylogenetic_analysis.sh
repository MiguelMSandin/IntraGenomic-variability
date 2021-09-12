#!/bin/bash

cd data/phylo/

ALIGNMENT="reference_consensus_checked.fasta"
TRIM_THRESHOLD="30"


# Trimming alignment
TRIMMED=${ALIGNMENT/.fasta/_trimed$TRIM_THRESHOLD.fasta}
trimal -in $ALIGNMENT -out $TRIMMED -gt 0.$TRIM_THRESHOLD

# And launching the phylogenetic inference
BS="1000"
THREADS="2"

# In order to replicate a bit the analysis, we perform one inference under the model GTR+CAT
OUTPUT=${TRIMMED/.fasta/_cat_BS$BS}
raxmlHPC-PTHREADS-SSE3 -T $THREADS -m GTRCAT -c 25 -p $RANDOM -x $(date +%s) -d -f a -N $BS -n $OUTPUT -s $TRIMMED

# And another inference under the model GTR+GAMMA
OUTPUT=${TRIMMED/.fasta/_gamma_BS$BS}
raxmlHPC-PTHREADS-SSE3 -T $THREADS -m GTRGAMMA -p $RANDOM -x $(date +%s) -d -f a -N $BS -n $OUTPUT -s $TRIMMED
