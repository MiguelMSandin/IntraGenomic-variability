#!/bin/bash


# Demultiplexing -----------------------------------------------------------------------------------

INPUT="data/minion/raw/raw.fasta"
TAGS="data/info/info_minion_demultiplex.txt"

CUTADAPT="$(which cutadapt) -e 0.1 --no-indels --discard-untrimmed "

[ ! -d "data/minion/demultiplexed/" ] && mkdir -p "data/minion/demultiplexed/"

while read TAG_NAME TAG_SEQ ; do
    OUT="data/minion/demultiplexed/${TAG_NAME}.fasta"
    LOG="data/minion/demultiplexed/${TAG_NAME}.log"
    cat "${INPUT}" | \
        ${CUTADAPT} -g "${TAG_SEQ}" -O "${#TAG_SEQ}" - 2> "${LOG}" > "${OUT}"
done < "${TAGS}"


# And BLAST ----------------------------------------------------------------------------------------

# Concatenating all demultiplexed
for FILE in $(ls data/minion/demultiplexed/*fasta)
do
	PREF=${FILE/.fasta/_}
	sed -e "s\>\>$PREF\g" $FILE >> data/minion/raw/demultiplexed.fasta
done
sed -i 's\data/minion/demultiplexed/\\g' data/minion/raw/demultiplexed.fasta

# setting BLASTn parameters
FILE="data/minion/raw/demultiplexed.fasta"
BLAST_TSV="data/minion/blast/demultiplexed_blast.tsv"
OUT_FMT="6 qseqid sseqid sacc stitle sscinames staxids sskingdoms sblastnames pident slen length mismatch gapopen qstart qend sstart send evalue bitscore"

DB="../DB/NCBI/"
DBnt="$DB/nt"

THREADS="2"

# Run blast
export BLASTDB=$DB
blastn -num_threads $THREADS -max_target_seqs 100 -evalue 1.00e-10 -query $FILE -out $BLAST_TSV -db $DBnt -outfmt "$OUT_FMT"
