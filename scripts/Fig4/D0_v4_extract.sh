#!/bin/bash

PRIMER_F="CCAGCASCYGCGGTAATTCC"
PRIMER_R="ACTTTCGTTCTTGATYRATGA"
ANTI_PRIMER_R="TCATYRATCAAGAACGAAAGT"

MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))

# Sanger concatenated reads ------------------------------------------------------------------------
INPUT="data/sanger/raw/raw_concatenated.fasta"
OUTPUT="${INPUT/.fasta/_V4.fasta}"
LOG="${INPUT/.fasta/_V4.log}"

CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH} -e 0.10"

cat ${INPUT} | \
	${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
	${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" > "${OUTPUT}"

# Select only Polycystine (Nassellaria and Spumellaria) reads
awk -F '\t' '$24 == "Spumellaria" { print $1}' "data/files/fig2_sanger.tsv" > "data/info/sanger_polycystines_reads.list"
awk -F '\t' '$24 == "Nassellaria" { print $1}' "data/files/fig2_sanger.tsv" >> "data/info/sanger_polycystines_reads.list"

sequenceSelect.py -f ${OUTPUT} -o "data/sanger/raw/polycystines_V4.fasta" -l "data/info/sanger_polycystines_reads.list" -a k

# And moving V4 files to a new directory
[ ! -d "data/sanger/V4" ] && mkdir -p "data/sanger/V4"
mv data/sanger/raw/*V4* data/sanger/V4/


# MinION demultiplexed reads -----------------------------------------------------------------------
INPUT="data/minion/raw/demultiplexed.fasta"
OUTPUT="${INPUT/.fasta/_V4.fasta}"
LOG="${INPUT/.fasta/_V4.log}"

CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH} -e 0.20"

cat ${INPUT} | \
	${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
	${CUTADAPT} -a "${ANTI_PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" > "${OUTPUT}"

# Select only Polycystine (Nassellaria and Spumellaria) reads
awk -F '\t' '$24 == "Spumellaria" { print $21"_"$22"_"$1}' "data/files/fig2_minion.tsv" > "data/info/minion_polycystines_reads.list"
awk -F '\t' '$24 == "Nassellaria" { print $21"_"$22"_"$1}' "data/files/fig2_minion.tsv" >> "data/info/minion_polycystines_reads.list"

sequenceSelect.py -f ${OUTPUT} -o "data/minion/raw/polycystines_V4.fasta" -l "data/info/minion_polycystines_reads.list" -a k

# And moving all V4 files to a new directory
[ ! -d "data/minion/V4" ] && mkdir -p "data/minion/V4"
mv data/minion/raw/*V4* data/minion/V4/

