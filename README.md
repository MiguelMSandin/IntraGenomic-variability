# IntraGenomic-variability
  
This repository gathers all scripts, data and information needed to replicate the following study: ADD_REFERENCE
  
## Dependencies
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)  
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)  
- [vsearch](https://github.com/torognes/vsearch)  
- [mafft](https://mafft.cbrc.jp/alignment/software/)  
- [trimAl](http://trimal.cgenomics.org/)  
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)  
- [mothur](https://mothur.org/)  
- [consension](https://microbiology.se/software/consension/)  
- [python](https://www.python.org/)  
    -   **Required modules**: argparse, Bio, sys.  
- [R](https://www.r-project.org/)  
    -   **Required packages**: dada2, Biostrings, ape, magrittr, data.table, dplyr, tidyr, ggplot2, treemapify, seqinr, colorRamps, vegan, circlize, forcats, cowplot, lulu (or mumu) and colorspace.  
### in-house dependencies
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
The folder **scripts** contains all the scripts used for the analysis and plotting of the figures, estructured with subdirectories for each figure.  
Folder **data** is subdivided into: **illumina**, **minion** and **sanger** (containing the raw fasta/q files used in this study: TO BE UPLOADED along with the permanent accession number at the repository ENA) in the subdirectories **raw**, along with contextualized data in the subdirectories **info**.  
  
All scripts should be self-explained and ready to be run. Tested on Ubuntu 20.04.
  
  
