#!/bin/bash

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# -------------------- Check alignments manualy before estimating the entropy!! --------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

for FILE in $(ls data/entropy/*.fasta)
do
	alignmentEntropy.py -f $FILE -g
done


# And moving the fasta files to different folders to keep an structure -----------------------------

[ ! -d "data/entropy/minion" ] && mkdir -p "data/entropy/minion"
mv data/entropy/minion*fasta data/entropy/minion

[ ! -d "data/entropy/sanger" ] && mkdir -p "data/entropy/sanger"
mv data/entropy/sanger*fasta data/entropy/sanger


rm -f tmp/*

