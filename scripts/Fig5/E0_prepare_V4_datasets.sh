#!/bin/bash

[ ! -d "data/entropy/V4" ] && mkdir -p "data/entropy/V4"

# Creating the reference files ---------------------------------------------------------------------
# Creating a reference of Nassellaria_V4
sequenceSelect.py -f "../DB/PR2/pr2_version_4.14.0_SSU_taxo_long_V4.fasta" -o "data/entropy/V4/pr2_V4_nass_all.fasta" -p "Nassellaria" -a k
grep -E "Orosphaeridae|Collophidiidae|Collosphaeridae|Sphaerozoidae|Collodaria" "data/entropy/V4/pr2_V4_nass_all.fasta" > "tmp/list"
sed -i 's/>//g' "tmp/list"
sequenceSelect.py -f "data/entropy/V4/pr2_V4_nass_all.fasta" -o "data/entropy/V4/pr2_V4_nass_noClean.fasta" -l "tmp/list" -a r
rm -f tmp/list

# Creating a reference of Spumellaria_V4
sequenceSelect.py -f "../DB/PR2/pr2_version_4.14.0_SSU_taxo_long_V4.fasta" -o "data/entropy/V4/pr2_V4_spum_noClean.fasta" -p "Spumellaria" -a k

# Removing duplicates and moving to the entropy folder

fastaClean.py -f "data/entropy/V4/pr2_V4_nass_noClean.fasta" -o "data/entropy/V4/ref_V4_nass.fasta"
sed -i 's/^>/>ref_/g' "data/entropy/V4/ref_V4_nass.fasta"
fastaClean.py -f "data/entropy/V4/pr2_V4_spum_noClean.fasta" -o "data/entropy/V4/ref_V4_spum.fasta"
sed -i 's/^>/>ref_/g' "data/entropy/V4/ref_V4_spum.fasta"

rm -f data/entropy/V4/pr2_V4_nass_all.fasta data/entropy/V4/pr2_V4_nass_noClean.fasta data/entropy/V4/pr2_V4_spum_noClean.fasta


# Separate sanger files ----------------------------------------------------------------------------

sequenceSelect.py -f "data/sanger/V4/polycystines_V4.fasta" -o "data/entropy/V4/sanger_V4_nass.fasta" -p "Vil" -a k
sed -i 's/^>/>sanger_/g' "data/entropy/V4/sanger_V4_nass.fasta"
sequenceSelect.py -f "data/sanger/V4/polycystines_V4.fasta" -o "data/entropy/V4/sanger_V4_spum.fasta" -p "Mge17" -a k
sed -i 's/^>/>sanger_/g' "data/entropy/V4/sanger_V4_spum.fasta"


# Separate minION files ----------------------------------------------------------------------------
sequenceSelect.py -f "data/minion/V4/polycystines_V4.fasta" -o "data/entropy/V4/minion_V4_nass.fasta" -p "Vil" -a k
sed -i 's/^>/>minion_/g' "data/entropy/V4/minion_V4_nass.fasta"
sequenceSelect.py -f "data/minion/V4/polycystines_V4.fasta" -o "data/entropy/V4/minion_V4_spum.fasta" -p "Mge17" -a k
sed -i 's/^>/>minion_/g' "data/entropy/V4/minion_V4_spum.fasta"


# Separate illumina dada2 file ---------------------------------------------------------------------
sequenceSelect.py -f "data/illumina/dada2/default_pipeline/out/ASVs_taxo.fasta" -o "tmp/tmp" -p "Nassellaria" -a k
grep -E "Orosphaeridae|Collophidiidae|Collosphaeridae|Sphaerozoidae|Collodaria" "tmp/tmp" > "tmp/list"
sed -i 's/>//g' "tmp/list"
sequenceSelect.py -f "tmp/tmp" -o "data/entropy/V4/illumina_V4_nass.fasta" -l "tmp/list" -a r
sed -i 's/^>/>illumina_/g' "data/entropy/V4/illumina_V4_nass.fasta"
rm -f tmp/list tmp/tmp
sequenceSelect.py -f "data/illumina/dada2/default_pipeline/out/ASVs_taxo.fasta" -o "data/entropy/V4/illumina_V4_spum.fasta" -p "Spumellaria" -a k
sed -i 's/^>/>illumina_/g' "data/entropy/V4/illumina_V4_spum.fasta"


# And finally merge all different datasets for Nassellaria and for Spumellaria

cat data/entropy/V4/*nass* > "data/entropy/V4/V4_nass.fasta"
cat data/entropy/V4/*spum* > "data/entropy/V4/V4_spum.fasta"

rm -f data/entropy/V4/illumina* data/entropy/V4/minion* data/entropy/V4/ref* data/entropy/V4/sanger*

