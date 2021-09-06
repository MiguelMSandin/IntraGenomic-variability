#!/bin/bash

VSEARCH=$(which vsearch)
FILE="data/illumina/raw/R1.fastq.gz"
OUTPUT="data/illumina/raw/encoding_R1.log"

# Check quality encoding (33 or 64?)
"${VSEARCH}" \
	--fastq_chars ${FILE} 2> ${OUTPUT}

grep -e "fastq_ascii" $OUTPUT
