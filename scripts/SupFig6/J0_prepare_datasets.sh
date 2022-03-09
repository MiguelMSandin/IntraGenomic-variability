#!/bin/bash

# Before running these scripts you should run first the script 'F0_extract_sequences_and_align.sh' from Figure 6


# Separate sanger files and remove columns composed of only gaps -----------------------------------
[ ! -d "data/entropy" ] && mkdir -p "data/entropy"

sequenceSelect.py -f "data/sanger/alignments/raw_concatenated_polycystines_align_clean.fasta" -o "tmp/tmp" -p "Vil325" -a k
trimal -noallgaps -in "tmp/tmp" -out "data/entropy/sanger_rDNA_nass_Vil325.fasta"

sequenceSelect.py -f "data/sanger/alignments/raw_concatenated_polycystines_align_clean.fasta" -o "tmp/tmp" -p "Vil496" -a k
trimal -noallgaps -in "tmp/tmp" -out "data/entropy/sanger_rDNA_nass_Vil496.fasta"

sequenceSelect.py -f "data/sanger/alignments/raw_concatenated_polycystines_align_clean.fasta" -o "tmp/tmp" -p "Mge17-81" -a k
trimal -noallgaps -in "tmp/tmp" -out "data/entropy/sanger_rDNA_spum_Mge17-81.fasta"

sequenceSelect.py -f "data/sanger/alignments/raw_concatenated_polycystines_align_clean.fasta" -o "tmp/tmp" -p "Mge17-82" -a k
trimal -noallgaps -in "tmp/tmp" -out "data/entropy/sanger_rDNA_spum_Mge17-82.fasta"


# Separate minION files and remove columns composed of only gaps -----------------------------------
sequenceSelect.py -f "data/minion/alignments/demultiplexed_polycystines_align.fasta" -o "tmp/tmp" -p "Vil325" -a k
trimal -noallgaps -in "tmp/tmp" -out "data/entropy/minion_rDNA_nass_Vil325.fasta"

sequenceSelect.py -f "data/minion/alignments/demultiplexed_polycystines_align.fasta" -o "tmp/tmp" -p "Vil496" -a k
trimal -noallgaps -in "tmp/tmp" -out "data/entropy/minion_rDNA_nass_Vil496.fasta"

sequenceSelect.py -f "data/minion/alignments/demultiplexed_polycystines_align.fasta" -o "tmp/tmp" -p "Mge17-81" -a k
trimal -noallgaps -in "tmp/tmp" -out "data/entropy/minion_rDNA_spum_Mge17-81.fasta"

sequenceSelect.py -f "data/minion/alignments/demultiplexed_polycystines_align.fasta" -o "tmp/tmp" -p "Mge17-82" -a k
trimal -noallgaps -in "tmp/tmp" -out "data/entropy/minion_rDNA_spum_Mge17-82.fasta"
