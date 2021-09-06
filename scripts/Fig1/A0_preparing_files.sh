#!/bin/bash

# For ONT data -------------------------------------------------------------------------------------

# Define variables
FASTQ="data/minion/raw/raw.fastq.gz"
OUT="data/minion/raw/raw.fasta"

# Uncompress
gunzip -k $FASTQ

# Convert fastQ into fastQ
cat ${FASTQ/.gz/} | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > $OUT


# For Sanger sequencing data -----------------------------------------------------------------------

# Define variables
FASTA="data/sanger/raw/raw.fasta.gz"

# Uncompress
gunzip -k $FASTA

# Convert the multiLine fasta file into a single line fasta
multi2linefasta.py -f ${FASTA/.gz/} -o tmp
mv tmp ${FASTA/.gz/}

# Extract sequences for the different primers and cells
[ ! -d "tmp" ] && mkdir -p "tmp"
for PRIMER in Nass_S30 Nass_S69f Nass_d1r Spum_S30 Spum_S69 Spum_d1r
do
	grep -A 1 $PRIMER ${FASTA/.gz/} > tmp/$PRIMER
	sed -i "s/$PRIMER//g" tmp/$PRIMER
	sed -i 's/>_/>/g' tmp/$PRIMER
done

# Concatenate sequences belonging to the same replicate
fastaConcat.py -f tmp/Nass_S30 tmp/Nass_S69f tmp/Nass_d1r -o tmp/nass.fasta
fastaConcat.py -f tmp/Spum_S30 tmp/Spum_S69 tmp/Spum_d1r -o tmp/spum.fasta
cat tmp/nass.fasta tmp/spum.fasta > ${FASTA/.fasta.gz/_concatenated.fasta}

# Remove temporary files
rm -f tmp/*


# Download PR2 database and preparing it for use ---------------------------------------------------
[ ! -d "~/DB/PR2" ] && mkdir -p "~/DB/PR2"

cd ~/DB/PR2

wget https://github.com/pr2database/pr2database/releases/download/v4.14.0/pr2_version_4.14.0_SSU_taxo_long.fasta.gz

gunzip -k pr2_version_4.14.0_SSU_taxo_long.fasta.gz

# Extract V4 region
PRIMER_F="CCAGCASCYGCGGTAATTCC"
PRIMER_R="ACTTTCGTTCTTGATYRATGA"
ANTI_PRIMER_R="TCATYRATCAAGAACGAAAGT"

MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))

INPUT="pr2_version_4.14.0_SSU_taxo_long.fasta"
OUTPUT="${INPUT/.fasta/_V4.fasta}"
LOG="${INPUT/.fasta/_V4.log}"

CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${MIN_LENGTH} -e 0.10"

cat $INPUT | \
     $CUTADAPT -g $PRIMER_F -O $MIN_F - 2> $LOG | \
     $CUTADAPT -a $ANTI_PRIMER_R -O $MIN_R - 2>> $LOG > $OUTPUT

# And the correspondent version for DADA2
wget https://github.com/pr2database/pr2database/releases/download/v4.14.0/pr2_version_4.14.0_SSU_dada2.fasta.gz
gunzip -k pr2_version_4.14.0_SSU_dada2.fasta.gz

INPUT="pr2_version_4.14.0_SSU_dada2.fasta"
OUTPUT="${INPUT/.fasta/_V4.fasta}"
LOG="${INPUT/.fasta/_V4.log}"

cat $INPUT | \
     $CUTADAPT -g $PRIMER_F -O $MIN_F - 2> $LOG | \
     $CUTADAPT -a $ANTI_PRIMER_R -O $MIN_R - 2>> $LOG > $OUTPUT

rm -f pr2_version_4.14.0_SSU_dada2.fasta


# Download NCBI taxonomy database ------------------------------------------------------------------
[ ! -d "../NCBI" ] && mkdir -p "../NCBI"
cd ../NCBI

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*
rm -f nt-nucl-metadata.json
for FILE in $(ls nt*tar.gz)
do
	echo "  Extracting $FILE"
	tar -xvzf $FILE
	rm -f $FILE
done
rm nt*tar.gz

# And the taxonomy files
wget ftp://ftp.ncbi.nih.gov/blast/db/taxdb.tar.gz
tar -xvzf taxdb.tar.gz
rm -f taxdb.tar.gz


cd ../../IntraGenomic-variability
