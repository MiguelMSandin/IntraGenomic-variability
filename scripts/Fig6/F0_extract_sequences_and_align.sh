#!/bin/bash

# Sanger sequences ---------------------------------------------------------------------------------

# Create a directory
[ ! -d "data/sanger/alignments" ] && mkdir -p "data/sanger/alignments"

# Select only Polycystine (Nassellaria and Spumellaria) reads
sequenceSelect.py -f "data/sanger/raw/raw_concatenated.fasta" -o "data/sanger/alignments/raw_concatenated_polycystines.fasta" -l "data/info/sanger_polycystines_reads.list" -a k

# Build a solid alignment (this might take a while...)
mafft --globalpair --maxiterate 1000 "data/sanger/alignments/raw_concatenated_polycystines.fasta" > "data/sanger/alignments/raw_concatenated_polycystines_align.fasta"

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# ------------------------------------ Check alignment manualy! ------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------  There are sequences asign to Nassellaria or Spumellaria --------------------
# ------------------- that do not share such assignment in the other replicates! -------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

echo -e "Mge17-82_BC11_G8\nVil325_BC02_E4\nVil496_BC15_D6" > tmp/tmp

sequenceSelect.py -f "data/sanger/alignments/raw_concatenated_polycystines_align.fasta" -o "data/sanger/alignments/raw_concatenated_polycystines_align_clean.fasta" -l tmp/tmp -a r


# MinION sequences ---------------------------------------------------------------------------------

# Create a directory
[ ! -d "data/minion/alignments" ] && mkdir -p "data/minion/alignments"

# Select only Polycystine (Nassellaria and Spumellaria) reads
sequenceSelect.py -f "data/minion/raw/demultiplexed.fasta" -o "data/minion/alignments/demultiplexed_polycystines.fasta" -l "data/info/minion_polycystines_reads.list" -a k

# Build a solid alignment (this might take a while...)
mafft --globalpair --maxiterate 1000 --add "data/minion/alignments/demultiplexed_polycystines.fasta" "data/sanger/alignments/raw_concatenated_polycystines_align_clean.fasta" > "tmp/tmp.fasta"

# And remove Sanger sequences
grep ">" "data/sanger/alignments/raw_concatenated_polycystines_align_clean.fasta" > tmp/tmp
sed -i 's/>//g' tmp/tmp

sequenceSelect.py -f tmp/tmp.fasta -o tmp/tmp2.fasta -l tmp/tmp -a r
trimal -noallgaps -in "tmp/tmp2.fasta" -out "data/minion/alignments/demultiplexed_polycystines_align.fasta"

rm -f tmp/*
