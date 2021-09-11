#!/bin/bash

SAMPLES=$(cut -f1 illumina/info/info_tags_cutadapt.txt)

[ ! -d "data/illumina/merged/" ] && mkdir -p "data/illumina/merged/"

for FILE in $SAMPLES
do
	VSEARCH=$(which vsearch)
	THREADS=2
	ENCODING=$(grep -e "fastq_ascii" "data/illumina/raw/encoding_R1.log" | cut -d " " -f 7)
	FORWARD="data/illumina/demultiplexed/${FILE}_R1_trimmed.fastq.gz"
	REVERSE="data/illumina/demultiplexed/${FILE}_R2_trimmed.fastq.gz"
	OUTPUT="data/illumina/merged/${FILE}_merged.fastq"
	
	# Merge read pairs
	"${VSEARCH}" \
		--threads ${THREADS} \
		--fastq_mergepairs ${FORWARD} \
		--reverse ${REVERSE} \
		--fastq_ascii ${ENCODING} \
		--fastqout ${OUTPUT} \
		--fastq_allowmergestagger \
		--quiet 2>> ${OUTPUT/.fastq/.log}
	
	sed -n '1~4s/^@/>/p;2~4p' ${OUTPUT} > ${OUTPUT/.fastq/.fasta}
done
