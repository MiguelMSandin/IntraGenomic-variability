#!/bin/bash

# Sanger sequences ---------------------------------------------------------------------------------

# Create a directory
[ ! -d "data/sanger/consensus" ] && mkdir -p "data/sanger/consensus"

# Select only Polycystine (Nassellaria and Spumellaria) reads
sequenceSelect.py -f "data/sanger/raw/raw_concatenated.fasta" -o "data/sanger/consensus/raw_concatenated_polycystines.fasta" -l "data/sanger/blast/polycystines_reads.list" -a k

# Build a solid alignment (this might take a while...)
mafft --globalpair --maxiterate 1000 "data/sanger/consensus/raw_concatenated_polycystines.fasta" > "data/sanger/consensus/raw_concatenated_polycystines_align.fasta"


# MinION sequences ---------------------------------------------------------------------------------

# Create a directory
[ ! -d "data/minion/consensus" ] && mkdir -p "data/minion/consensus"

# Select only Polycystine (Nassellaria and Spumellaria) reads
sequenceSelect.py -f "data/minion/raw/demultiplexed.fasta" -o "data/minion/consensus/demultiplexed_polycystines.fasta" -l "data/minion/blast/polycystines_reads.list" -a k

# Build a solid alignment (this might take a while...)
mafft --globalpair --maxiterate 1000 "data/minion/consensus/demultiplexed_polycystines.fasta" > "data/minion/consensus/demultiplexed_polycystines_align.fasta"

