#!/bin/bash

# Sanger V4
INPUT="data/sanger/V4/polycystines_V4.fasta"
OUTPUT=${INPUT/.fasta/_similarities.tsv}

vsearch \
	--threads 2 \
	--allpairs_global $INPUT \
	--self \
	--id 0.20 \
	--userfields query+target+id \
	--userout $OUTPUT


# MinION V4
INPUT="data/minion/V4/polycystines_V4.fasta"
OUTPUT=${INPUT/.fasta/_similarities.tsv}

vsearch \
	--threads 2 \
	--allpairs_global $INPUT \
	--self \
	--id 0.20 \
	--userfields query+target+id \
	--userout $OUTPUT

