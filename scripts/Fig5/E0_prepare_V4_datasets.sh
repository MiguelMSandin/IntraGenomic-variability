#!/bin/bash

[ ! -d "data/ref" ] && mkdir -p "data/ref"

# Creating the reference files ---------------------------------------------------------------------
# Creating a reference of Nassellaria_V4
sequenceSelect.py -f "../DB/PR2/pr2_version_4.14.0_SSU_taxo_long_V4.fasta" -o "data/ref/pr2_V4_nass_all.fasta" -p "Nassellaria" -a k
grep -E "Orosphaeridae|Collophidiidae|Collosphaeridae|Sphaerozoidae|Collodaria" "data/ref/pr2_V4_nass_all.fasta" > "tmp/list"
sed -i 's/>//g' "tmp/list"
sequenceSelect.py -f "data/ref/pr2_V4_nass_all.fasta" -o "data/ref/ref_V4_nass.fasta" -l "tmp/list" -a r
rm -f tmp/list

# Creating a reference of Spumellaria_V4
sequenceSelect.py -f "../DB/PR2/pr2_version_4.14.0_SSU_taxo_long_V4.fasta" -o "data/ref/ref_V4_spum.fasta" -p "Spumellaria" -a k

# Removing duplicates and moving to the entropy folder
[ ! -d "data/entropy" ] && mkdir -p "data/entropy"

fastaClean.py -f "data/ref/pr2_V4_nass.fasta" -o "data/entropy/pr2_V4_nass_noDup.fasta"
fastaClean.py -f "data/ref/pr2_V4_spum.fasta" -o "data/entropy/pr2_V4_spum_noDup.fasta"


# Separate sanger files ----------------------------------------------------------------------------
sequenceSelect.py -f "data/sanger/V4/polycystines_V4.fasta" -o "data/entropy/sanger_V4_nass_Vil325.fasta" -p "Vil325" -a k
sequenceSelect.py -f "data/sanger/V4/polycystines_V4.fasta" -o "data/entropy/sanger_V4_nass_Vil496.fasta" -p "Vil496" -a k
sequenceSelect.py -f "data/sanger/V4/polycystines_V4.fasta" -o "data/entropy/sanger_V4_spum_Mge17-81.fasta" -p "Mge17-81" -a k
sequenceSelect.py -f "data/sanger/V4/polycystines_V4.fasta" -o "data/entropy/sanger_V4_spum_Mge17-82.fasta" -p "Mge17-82" -a k


# Separate minION files ----------------------------------------------------------------------------
sequenceSelect.py -f "data/minion/V4/polycystines_V4.fasta" -o "data/entropy/minion_V4_nass.fasta" -p "Vil" -a k
sequenceSelect.py -f "data/minion/V4/polycystines_V4.fasta" -o "data/entropy/minion_V4_spum.fasta" -p "Mge17" -a k


# Separate illumina dada2 file ---------------------------------------------------------------------
sequenceSelect.py -f "data/illumina/dada2/default_pipeline/out/ASVs_taxo.fasta" -o "tmp/tmp" -p "Nassellaria" -a k
grep -E "Orosphaeridae|Collophidiidae|Collosphaeridae|Sphaerozoidae|Collodaria" "tmp/tmp" > "tmp/list"
sed -i 's/>//g' "tmp/list"
sequenceSelect.py -f "tmp/tmp" -o "data/entropy/illumina_V4_nass.fasta" -l "tmp/list" -a r
rm -f tmp/list tmp/tmp
sequenceSelect.py -f "data/illumina/dada2/default_pipeline/out/ASVs_taxo.fasta" -o "data/entropy/illumina_V4_spum.fasta" -p "Spumellaria" -a k


# Aligning all files to be analyzed ----------------------------------------------------------------
for FILE in $(ls data/entropy/*fasta)
do
	mafft --localpair --maxiterate 1000 $FILE > ${FILE/.fasta/_align.fasta}
done

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
# ----------------------------------- Check alignments manualy!! -----------------------------------
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

for FILE in $(ls data/entropy/*align.fasta)
do
	alignmentEntropy.py -f $FILE -g -r 10
done


# And moving the fasta files to different folders to keep an structure -----------------------------

for FOLDER in ref sanger illumina minion
do
	[ ! -d "data/entropy/$FOLDER" ] && mkdir -p "data/entropy/$FOLDER"
	mv data/entropy/${FOLDER}_*fasta "data/entropy/$FOLDER/"
done
mv data/entropy/pr2_*fasta data/entropy/ref/
