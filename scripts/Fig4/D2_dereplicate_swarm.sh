#!/bin/bash

# Sanger V4 ----------------------------------------------------------------------------------------
INPUT="data/sanger/V4/polycystines_V4.fasta"
DEREP=${INPUT/.fasta/_derep.fasta}
LOG=${INPUT/.fasta/_derep.log}

# Dereplicating
vsearch \
	--derep_fulllength $INPUT \
	--sizeout \
	--fasta_width 0 \
	--relabel_keep \
	--output $DEREP 2> $LOG

# swarm
for DIST in 1 2 3
do
	swarm -t 2 -z \
		-d $DIST \
		-s ${DEREP/_derep.fasta/_swarm_d$DIST.stats} \
		-w ${DEREP/_derep.fasta/_swarm_d$DIST.fasta} \
		-o ${DEREP/_derep.fasta/_swarm_d$DIST.log} \
		$DEREP
done

# Moving all swarm output files to a new directory
[ ! -d "data/sanger/V4/swarm" ] && mkdir -p "data/sanger/V4/swarm"
mv data/sanger/V4/*_swarm_* data/sanger/V4/swarm/

# Calculate similarities of the clusters
for FILE in $(ls data/sanger/V4/swarm/*fasta)
do
	vsearch --id 0.2 \
		--allpairs_global $FILE \
		--userfields query+target+id \
		--userout ${FILE/.fasta/.similarities}
	
	# And add self hits
	grep ">" $FILE > tmp/tmp
	sed -i 's/>//g' tmp/tmp
	while read SEQ; do echo -e "${SEQ/>/}\t${SEQ/>/}\t100" >> ${FILE/.fasta/.similarities}; done < "tmp/tmp"
done

rm -f tmp/tmp

# MinION V4 ----------------------------------------------------------------------------------------
INPUT="data/minion/V4/polycystines_V4.fasta"
DEREP=${INPUT/.fasta/_derep.fasta}
LOG=${INPUT/.fasta/_derep.log}

# Dereplicating
vsearch \
	--derep_fulllength $INPUT \
	--sizeout \
	--fasta_width 0 \
	--relabel_keep \
	--output $DEREP 2> $LOG

# swarm
for DIST in $(seq 1 10)
do
	swarm -t 2 -z \
		-d $DIST \
		-s ${DEREP/_derep.fasta/_swarm_d$DIST.stats} \
		-w ${DEREP/_derep.fasta/_swarm_d$DIST.fasta} \
		-o ${DEREP/_derep.fasta/_swarm_d$DIST.log} < $DEREP
done

# Moving all swarm output files to a new directory
[ ! -d "data/minion/V4/swarm" ] && mkdir -p "data/minion/V4/swarm"
mv data/minion/V4/*_swarm_* data/minion/V4/swarm

# And calculate similarities of the clusters
for FILE in $(ls data/minion/V4/swarm/*fasta)
do
	vsearch --id 0.2 \
		--allpairs_global $FILE \
		--userfields query+target+id \
		--userout ${FILE/.fasta/.similarities}
done
