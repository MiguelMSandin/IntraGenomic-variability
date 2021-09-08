#!/bin/bash

ALIGNMENT="data/phylo/rads_ref_align_consensus.fasta"
TRIM_THRESHOLD="05"


# Trimming checked alignment
TRIMMED=${FILE/.fasta/_trimed$TRIM_THRESHOLD.fasta}
trimal -in $ALIGNMENT -out $TRIMMED -gt $TRIM_THRESHOLD


# And launching the phylogenetic inference
BS="100"
THREADS="2"
SEED=$(date +%s)

# In order to replicate a bit the analysis, we perform one inference under the model GTR+CAT
OUTPUT=${FILE/.fasta/_cat_BS$BS}
raxmlHPC-PTHREADS-SSE3 -T $THREADS -m GTRCAT -c 25 -p $RANDOM -x $SEED -d -f a -N $BS -n $OUTPUT -s $FILE

# And another inference under the model GTR+GAMMA
OUTPUT=${FILE/.fasta/_gamma_BS$BS}
raxmlHPC-PTHREADS-SSE3 -T $THREADS -m GTRGAMMA -p $RANDOM -x $SEED -d -f a -N $BS -n $OUTPUT -s $FILE
