#!/bin/bash

# Get the consensus sequences from sanger
sed 's/^>/>sanger_/g' "data/sanger/consensus/raw_concatenated_polycystines_consensus.fasta" > tmp/tmp1
sed 's/^>/>minion_/g' "data/minion/consensus/demultiplexed_polycystines_consensus.fasta" > tmp/tmp2

cat tmp/tmp1 tmp/tmp2 > "data/phylo/consensus.fasta"

# And add the consensus sequence to the reference alignment
mafft --thread 2 --add "data/phylo/consensus.fasta" "data/phylo/reference.fasta" > "data/phylo/reference_consensus.fasta"

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# ------------------------------------ Check alignment manualy! ------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

rm tmp/*

