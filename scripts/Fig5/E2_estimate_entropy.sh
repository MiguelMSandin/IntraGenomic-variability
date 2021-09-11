#!/bin/bash

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# -------------------- Check alignments manualy before estimating the entropy!! --------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Separate again the files -------------------------------------------------------------------------

# For Nassellaria
sequenceSelect.py -f "data/entropy/V4/V4_nass_align.fasta" -o "data/entropy/V4/align_sanger_V4_nass_Vil325.fasta" -p "sanger_Vil325" -a k
sequenceSelect.py -f "data/entropy/V4/V4_nass_align.fasta" -o "data/entropy/V4/align_sanger_V4_nass_Vil496.fasta" -p "sanger_Vil496" -a k
sequenceSelect.py -f "data/entropy/V4/V4_nass_align.fasta" -o "data/entropy/V4/align_illumina_V4_nass.fasta" -p "illumina_" -a k
sequenceSelect.py -f "data/entropy/V4/V4_nass_align.fasta" -o "data/entropy/V4/align_minion_V4_nass.fasta" -p "minion_" -a k
sequenceSelect.py -f "data/entropy/V4/V4_nass_align.fasta" -o "data/entropy/V4/align_ref_V4_nass.fasta" -p "ref_" -a k

# For Spumellaria
sequenceSelect.py -f "data/entropy/V4/V4_spum_align.fasta" -o "data/entropy/V4/align_sanger_V4_spum_Mge17-81.fasta" -p "sanger_Mge17-81" -a k
sequenceSelect.py -f "data/entropy/V4/V4_spum_align.fasta" -o "data/entropy/V4/align_sanger_V4_spum_Mge17-82.fasta" -p "sanger_Mge17-82" -a k
sequenceSelect.py -f "data/entropy/V4/V4_spum_align.fasta" -o "data/entropy/V4/align_illumina_V4_spum.fasta" -p "illumina_" -a k
sequenceSelect.py -f "data/entropy/V4/V4_spum_align.fasta" -o "data/entropy/V4/align_minion_V4_spum.fasta" -p "minion_" -a k
sequenceSelect.py -f "data/entropy/V4/V4_spum_align.fasta" -o "data/entropy/V4/align_ref_V4_spum.fasta" -p "ref_" -a k


# And finally estimate the entropy -----------------------------------------------------------------
for FILE in $(ls data/entropy/V4/align_*)
do
	alignmentEntropy.py -f $FILE -g -r 10
done

rm -f data/entropy/V4/V4_nass.fasta data/entropy/V4/V4_spum.fasta data/entropy/V4/align*fasta

