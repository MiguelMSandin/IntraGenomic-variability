[![Published in Environmental Microbiology](https://img.shields.io/badge/published%20in-Environmental%20Microbiology-blue.svg)](https://doi.org/10.1111/1462-2920.16081)

# IntraGenomic-variability
  
This repository gathers all scripts, data and information needed to replicate the following study:  
Sandin MM, Romac S, Not F (**2022**) Intra-genomic rRNA gene variability of Nassellaria and Spumellaria (Rhizaria, Radiolaria) assessed by Sanger, MinION and Illumina sequencing. *Environ Microbiol*. 24(7):2979-2993. doi: [10.1111/1462-2920.16081](https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1462-2920.16081)  
  
## Dependencies
- [python](https://www.python.org/)  
    -   **Required modules**: argparse, Bio, sys, re, statistics, numpy, pandas, math.  
- [R](https://www.r-project.org/)  
    -   **Required packages**: dada2, Biostrings, ape, magrittr, data.table, dplyr, tidyr, ggplot2, treemapify, seqinr, colorRamps, vegan, circlize, forcats, cowplot, lulu (or mumu) and colorspace.  
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)  
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  
- [vsearch](https://github.com/torognes/vsearch)  
- [mafft](https://mafft.cbrc.jp/alignment/software/)  
- [trimAl](http://trimal.cgenomics.org/)  
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)  
- [mothur](https://mothur.org/)  
- [consension](https://microbiology.se/software/consension/)  
- [minimap2](https://github.com/lh3/minimap2)  
- [Racon](https://github.com/isovic/racon)  
- An alignment editor software, such as [aliview](https://ormbunkar.se/aliview/) or [seaview](http://doua.prabi.fr/software/seaview)  
### in-house dependencies
Accessible in the [random](https://github.com/MiguelMSandin/random) repository:  
- [multi2linefasta.py](https://github.com/MiguelMSandin/random/blob/main/fasta/multi2linefasta.py)  
- [fastaConcat.py](https://github.com/MiguelMSandin/random/blob/main/fasta/fastaConcat.py)  
- [sequenceSelect.py](https://github.com/MiguelMSandin/random/blob/main/fasta/sequenceSelect.py)  
- [fastaStats.py](https://github.com/MiguelMSandin/random/blob/main/fasta/fastaStats.py)  
- [alignmentConsensus.py](https://github.com/MiguelMSandin/random/blob/main/fasta/alignmentConsensus.py)  
- [alignmentEntropy.py](https://github.com/MiguelMSandin/random/blob/main/fasta/alignmentEntropy.py)  
- [fastaClean.py](https://github.com/MiguelMSandin/random/blob/main/fasta/fastaClean.py)  
### Optional (yet recomended)  
- [Rstudio](https://rstudio.com/products/rstudio/download/) 
  
## Structure of the repository
The folder **scripts** contains all the scripts used for the analysis and plotting of the figures presented in the [study](https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1462-2920.16081) mentioned above, estructured with subdirectories for each figure.  
Folder **data/info** contains the contextualized information in order to perform and analyze the raw data (e.g.; tags for demultiplexing), that can be permanently downloaded from figshare using the following DOI: https://doi.org/10.6084/m9.figshare.16922764. Raw sequencing results have also been deposited in SRA under the accession number: [PRJNA816840](https://www.ncbi.nlm.nih.gov/sra/PRJNA816840).  
  
All scripts should be self-explained and ready to be run. Tested on a laptop Ubuntu 20.04, with 4 processors and 16 GB of RAM.  

Raw figures were edited in [inkscape](https://inkscape.org/).  
  
