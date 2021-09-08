#!/bin/bash

[ ! -d "data/phylo" ] && mkdir -p "data/phylo"

# Select sequences of PR2 to build a reference phylogeny
sequenceSelect.py -f "../DB/PR2/pr2_version_4.14.0_SSU_taxo_long.fasta" -o "data/phylo/rads.fasta" -p "Radiolaria"
fastaClean.py -f "data/phylo/rads.fasta" -o "data/phylo/rads_clean.fasta" -l "1600"

# Select an outgroup
sequenceSelect.py -f "../DB/PR2/pr2_version_4.14.0_SSU_taxo_long.fasta" -o "tmp/tmp" -p "Phaeodaria"
fastaClean.py -f "tmp/tmp" -o "tmp/tmp2" -l "1600"

# create a reference file
cat tmp/tmp2 data/phylo/rads_clean.fasta > data/phylo/rads_ref.fasta

# And align it
mafft --thread 2 --retree 2 --maxiterate 1000 "data/phylo/rads_ref.fasta" > "data/phylo/rads_ref_align.fasta"


# Use now the consensus sequences
sed -e 's/^>/>minion_/g' "data/sanger/consensus/raw_concatenated_polycystines_consensus.fasta" > tmp/tmp1
sed -e 's/^>/>sanger_/g' "data/minion/consensus/demultiplexed_polycystines_consensus.fasta" > tmp/tmp2

cat tmp/tmp1 tmp/tmp2 > "data/phylo/consensus.fasta"

# And add the consensus sequence to the reference alignment
mafft --thread 2 --globalpair --maxiterate 100 --add "data/phylo/consensus.fasta" "data/phylo/rads_ref_align.fasta" > "data/phylo/rads_ref_align_consensus.fasta"
