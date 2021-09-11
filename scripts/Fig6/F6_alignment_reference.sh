#!/bin/bash

# Align reference sequences (this might take a while...)
mafft --thread 2 --maxiterate 1000 --globalpair "data/phylo/seqs_18S.fasta" > "data/phylo/seqs_18S_align.fasta"

mafft --thread 2 --maxiterate 1000 --globalpair "data/phylo/seqs_28S.fasta" > "data/phylo/seqs_28S_align.fasta"

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# ----------------------------------- Check alignments manualy!! -----------------------------------
# --------------------------------------------------------------------------------------------------
# --- Some 18S and 28S sequences retrieved are actually the full rDNA, please correct manually! ----
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Concatenate 18S and 28S
fastaConcat.py -f "data/phylo/seqs_18S_align_checked.fasta" "data/phylo/seqs_28S_align_checked.fasta" -o "data/phylo/reference.fasta" -a
