 
library(lulu)
library(data.table)

setwd("data/illumina/lulu/")

system("cp ../dada2/default_pipeline/out/ASVs_taxo.fasta .")
system("vsearch --threads 2 --usearch_global ASVs_taxo.fasta --db ASVs_taxo.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10")

abun <- fread("ASVs_abun.tsv", sep="\t")
simm <- fread("match_list.txt", sep="\t")

simm$V1 <- gsub("\\|.*", "", simm$V1)
simm$V2 <- gsub("\\|.*", "", simm$V2)

tmp <- abun$ASV
abun$ASV <- gsub("\\|.*", "", abun$ASV)
taxo <- data.frame(asv=abun$ASV, taxo=tmp) 
rownames(abun)  <- abun$ASV; abun$ASV <- NULL

curated <- lulu(abun, simm)
