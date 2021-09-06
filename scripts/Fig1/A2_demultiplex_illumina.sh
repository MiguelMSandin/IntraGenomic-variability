#!/bin/bash

PRIMER_F="CCAGCASCYGCGGTAATTCC"
PRIMER_R="ACTTTCGTTCTTGATYRATGA"
TAGS="data/illumina/info/info_tags_cutadapt.txt"

INPUT_R1="data/illumina/raw/R1.fastq.gz"
INPUT_R2="data/illumina/raw/R2.fastq.gz"

[ ! -d "data/illumina/demultiplexed/" ] && mkdir -p "data/illumina/demultiplexed/"

while read SPLE_NAME TAG_SEQ; do
	echo "  Working on ${SPLE_NAME}"
	
	TMP_TAG_F=$(mktemp)
	TMP_TAG_R=$(mktemp)
	TMP_CUT_F1=$(mktemp)
	TMP_CUT_F2=$(mktemp)
	TMP_CUT_R1=$(mktemp)
	TMP_CUT_R2=$(mktemp)
	
	OUT="data/illumina/demultiplexed/${SPLE_NAME}"
	LOG="data/illumina/demultiplexed/${SPLE_NAME}.log"
	
	cutadapt -g "${TAG_SEQ}" -O ${#TAG_SEQ} --discard-untrimmed --no-indels -o ${TMP_TAG_F} -p ${TMP_TAG_R} ${INPUT_R1} ${INPUT_R2}
	cutadapt -g "${PRIMER_F}" -G "${PRIMER_R}" --report=minimal --discard-untrimmed --minimum-length 100 --no-indels -o ${TMP_CUT_F1} -p ${TMP_CUT_R1} ${TMP_TAG_F} ${TMP_TAG_R}  1>> ${LOG}
	
	cutadapt -g "${TAG_SEQ}" -O ${#TAG_SEQ} --discard-untrimmed --no-indels -o ${TMP_TAG_F} -p ${TMP_TAG_R} ${INPUT_R2} ${INPUT_R1}
	cutadapt -g "${PRIMER_F}" -G "${PRIMER_R}" --report=minimal --discard-untrimmed --minimum-length 100 --no-indels -o ${TMP_CUT_F2} -p ${TMP_CUT_R2} ${TMP_TAG_F} ${TMP_TAG_R}  1>> ${LOG}
	
	cat ${TMP_CUT_F1} ${TMP_CUT_F2} | gzip > "${OUT}_R1_trimmed.fastq.gz"
	cat ${TMP_CUT_R1} ${TMP_CUT_R2} | gzip > "${OUT}_R2_trimmed.fastq.gz"
	
	rm -f "${TMP_TAG_F}" "${TMP_TAG_R}" "${TMP_CUT_F1}" "${TMP_CUT_F2}" "${TMP_CUT_R1}" "${TMP_CUT_R2}"
	
	echo "----------------------------------------"
done < "${TAGS}"

echo "Done"
