#!/bin/bash

# Sanger sequences ---------------------------------------------------------------------------------

ALIGNMENT="data/sanger/alignments/raw_concatenated_polycystines_align_clean.fasta"
OUTPUT="data/sanger/consensus/raw_concatenated_polycystines_consensus.fasta"

[ ! -d "data/sanger/consensus" ] && mkdir -p "data/sanger/consensus"

for CELL in Mge17-81 Mge17-82 Vil325 Vil496
do
	echo "Working on $CELL"
	sequenceSelect.py -f $ALIGNMENT -o tmp/$CELL -p $CELL -a k
	(( $(grep -c '>' tmp/$CELL) < 4 )) && echo "  Warning: There are less than 4 sequences matching $CELL"
	alignmentConsensus.py -f tmp/$CELL -o $OUTPUT -m -r -v
	rm -f tmp/$CELL
done
sed -i 's\tmp/\\g' $OUTPUT


# MinION sequences ---------------------------------------------------------------------------------

# Download script
wget https://microbiology.se/sw/consension.zip
unzip consension.zip
mv consension scripts/Fig6/
rm -f consension.zip

# Modify script to get the size/abundance of the OTUs
sed "194s/.*/\t@tmp = scalar(split(',', @otus[\$o]));\n\tpush(@save_names, \"@otu_names[\$o]_@tmp\");/" scripts/Fig6/consension > scripts/Fig6/F3_consension
chmod +x scripts/Fig6/F3_consension

# Start the clustering
ALIGNMENT="data/minion/alignments/demultiplexed_polycystines_align.fasta"
DIST_LIST="data/minion/consensus/clusters/demultiplexed_polycystines_align.phylip.opti_mcc.list"
MIN_CLUSTER_SIZE="3"
OUTPUT="data/minion/consensus/demultiplexed_polycystines_consensus.fasta"

scripts/Fig6/F3_consension -i $ALIGNMENT -t $DIST_LIST -o $OUTPUT -K $MIN_CLUSTER_SIZE -R T -T T

sed -i 's/-//g' $OUTPUT
sed -i 's/\.//g' $OUTPUT
