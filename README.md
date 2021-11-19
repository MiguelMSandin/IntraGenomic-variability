# IntraGenomic-variability
  
This repository gathers all scripts, data and information needed to replicate the following study:  
Sandin MM, Romac S, Not F (2021) **Intra-genomic rDNA gene variability of Nassellaria and Spumellaria (Rhizaria, Radiolaria) assessed by Sanger, MinION and Illumina sequencing**. bioRxiv 2021.10.05.463214; doi: [doi.org/10.1101/2021.10.05.463214](https://doi.org/10.1101/2021.10.05.463214)
  
## Dependencies
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)  
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  
- [vsearch](https://github.com/torognes/vsearch)  
- [mafft](https://mafft.cbrc.jp/alignment/software/)  
- [trimAl](http://trimal.cgenomics.org/)  
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)  
- [mothur](https://mothur.org/)  
- [minimap2](https://github.com/lh3/minimap2)  
- [Racon](https://github.com/isovic/racon)  
- [consension](https://microbiology.se/software/consension/)  
- [python](https://www.python.org/)  
    -   **Required modules**: argparse, Bio, sys, re, statistics, numpy, pandas, math.  
- [R](https://www.r-project.org/)  
    -   **Required packages**: dada2, Biostrings, ape, magrittr, data.table, dplyr, tidyr, ggplot2, treemapify, seqinr, colorRamps, vegan, circlize, forcats, cowplot, lulu (or mumu) and colorspace.  
- An alignment editor software, such as [aliview](https://ormbunkar.se/aliview/) or [seaview](http://doua.prabi.fr/software/seaview)  
### in-house dependencies
Accessible in [fasta-functions](https://github.com/MiguelMSandin/fasta-functions)  
- [multi2linefasta.py](https://github.com/MiguelMSandin/fasta-functions/blob/main/scripts/multi2linefasta.py)  
- [fastaConcat.py](https://github.com/MiguelMSandin/fasta-functions/blob/main/scripts/fastaConcat.py)  
- [sequenceSelect.py](https://github.com/MiguelMSandin/fasta-functions/blob/main/scripts/sequenceSelect.py)  
- [fastaStats.py](https://github.com/MiguelMSandin/fasta-functions/blob/main/scripts/fastaStats.py)  
- [alignmentConsensus.py](https://github.com/MiguelMSandin/fasta-functions/blob/main/scripts/alignmentConsensus.py)  
- [alignmentEntropy.py](https://github.com/MiguelMSandin/fasta-functions/blob/main/scripts/alignmentEntropy.py)  
- [fastaClean.py](https://github.com/MiguelMSandin/fasta-functions/blob/main/scripts/fastaClean.py)  
### Optional (yet recomended)  
- [Rstudio](https://rstudio.com/products/rstudio/download/) 
  
## Structure of the repository
The folder **scripts** contains all the scripts used for the analysis and plotting of the figures presented in the [study](https://doi.org/10.1101/2021.10.05.463214) mentioned above, estructured with subdirectories for each figure.  
Folder **data/info** contains the contextualized information in order to perform and analyze the raw data (e.g.; tags for demultiplexing), that can be permanently downloaded from figshare using the following DOI: https://doi.org/10.6084/m9.figshare.16922764.  
  
All scripts should be self-explained and ready to be run. Tested on a laptop Ubuntu 20.04, with 4 processors and 16 GB of RAM.  

Raw figures were edited in [inkscape](https://inkscape.org/).  
  
