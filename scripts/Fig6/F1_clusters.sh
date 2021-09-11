#!/bin/bash

# MinION sequences ---------------------------------------------------------------------------------
[ ! -d "data/minion/consensus" ] && mkdir -p "data/minion/consensus"

FILE="data/minion/alignments/demultiplexed_polycystines_align.fasta"

# Creating the distance tables
mothur "#dist.seqs(fasta=$FILE, processors=2)"

#---------------------------------------------------------------------------------------------------
##--------------------- Check the output with the script 'F2_explore_output.R' ---------------------
#---------------------------------------------------------------------------------------------------

## Using a very big 30% threshold

# Clustering
mothur "#dist.seqs(fasta=$FILE, output=lt, processors=10)"
mothur "#cluster(phylip=${FILE/.fasta/.phylip.dist}, cutoff=0.30)"

# And move file to the consensus folder 
[ ! -d "data/minion/consensus/clusters" ] && mkdir -p "data/minion/consensus/clusters"
mv data/minion/alignments/*dist data/minion/consensus/clusters/
mv data/minion/alignments/*phylip* data/minion/consensus/clusters/

rm -f mothur*
