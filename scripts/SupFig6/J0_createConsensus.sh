#!/bin/bash

sequenceSelect.py -f "data/sanger/alignments/raw_concatenated_polycystines_align_clean.fasta" -p Mge17-82 -o "raw_concatenated_polycystines_align_clean_mge17-82.fasta"

alignmentConsensus.py -f "raw_concatenated_polycystines_align_clean_mge17-82.fasta" -o "raw_concatenated_polycystines_align_clean_mge17-82_consensus.fasta" -t 10 -b 10

# Now import the consensus sequence to https://rnacentral.org/r2dt and edit the figure including the plots obtained 
