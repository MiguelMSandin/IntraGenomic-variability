#!/bin/bash

# Aligning all files to be analyzed ----------------------------------------------------------------
for FILE in $(ls data/entropy/V4/*fasta)
do
	echo "  Working on $FILE"
	mafft --globalpair --maxiterate 1000 $FILE > ${FILE/.fasta/_align.fasta}
	echo "----------------------------------------"
done
echo "Done"

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# ----------------------------------- Check alignments manualy!! -----------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

