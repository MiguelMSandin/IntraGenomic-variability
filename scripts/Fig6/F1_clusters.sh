#!/bin/bash

# Sanger sequences ---------------------------------------------------------------------------------

FILE="data/sanger/consensus/raw_concatenated_polycystines_align.fasta"

# Creating the distance tables
mothur "#dist.seqs(fasta=$FILE, processors=2)"

#---------------------------------------------------------------------------------------------------
##--------------------- Check the output with the script 'F2_explore_output.R' ---------------------
#---------------------------------------------------------------------------------------------------

## Using a classic 97% threshold

# Clustering
mothur "#dist.seqs(fasta=$FILE, output=lt, processors=10)"
mothur "#cluster(phylip=${FILE/.fasta/.phylip.dist}, cutoff=0.03)"

# Moving to a new subdirectory
[ ! -d "data/sanger/consensus/clusters" ] && mkdir -p "data/sanger/consensus/clusters"
mv data/sanger/consensus/*dist data/sanger/consensus/clusters/
mv data/sanger/consensus/*phylip* data/sanger/consensus/clusters/

# MinION sequences ---------------------------------------------------------------------------------

FILE="data/minion/consensus/demultiplexed_polycystines_align.fasta"

# Creating the distance tables
mothur "#dist.seqs(fasta=$FILE, processors=2)"

#---------------------------------------------------------------------------------------------------
##--------------------- Check the output with the script 'F2_explore_output.R' ---------------------
#---------------------------------------------------------------------------------------------------

## Using a very big 30% threshold

# Clustering
mothur "#dist.seqs(fasta=$FILE, output=lt, processors=10)"
mothur "#cluster(phylip=${FILE/.fasta/.phylip.dist}, cutoff=0.30)"
## 
[ ! -d "data/minion/consensus/clusters" ] && mkdir -p "data/minion/consensus/clusters"
mv data/minion/consensus/*dist data/minion/consensus/clusters/
mv data/minion/consensus/*phylip* data/minion/consensus/clusters/

rm -f mothur*

