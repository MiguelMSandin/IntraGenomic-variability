#!/bin/bash

# Sanger sequences ---------------------------------------------------------------------------------

ALIGNMENT="data/sanger/consensus/raw_concatenated_polycystines_align.fasta"
DIST_LIST="data/sanger/consensus/clusters/raw_concatenated_polycystines_align.phylip.opti_mcc.list"
MIN_CLUSTER_SIZE="3"
OUTPUT="${ALIGNMENT/_align.fasta/_consensus.fasta}"

scripts/Fig6/F3_consension_new -i $ALIGNMENT -t $DIST_LIST -o $OUTPUT -K $MIN_CLUSTER_SIZE -R T -T T

sed -i 's/-//g' $OUTPUT

# MinION sequences ---------------------------------------------------------------------------------

ALIGNMENT="data/minion/consensus/demultiplexed_polycystines_align.fasta"
DIST_LIST="data/minion/consensus/clusters/demultiplexed_polycystines_align.phylip.opti_mcc.list"
MIN_CLUSTER_SIZE="3"
OUTPUT="${ALIGNMENT/_align.fasta/_consensus.fasta}"

scripts/Fig6/F3_consension_new -i $ALIGNMENT -t $DIST_LIST -o $OUTPUT -K $MIN_CLUSTER_SIZE -R T -T T

sed -i 's/-//g' $OUTPUT
