 
library(lulu)
library(data.table)

setwd("/home/mmsandin/Desktop/PhD/0_Thesis/3.1_Chapter/github/data/illumina/lulu/")

abun <- fread("ASVs_abun.tsv", sep="\t")
simm <- fread("match_list.txt", sep="\t")

simm$V1 <- gsub("\\|.*", "", simm$V1)
simm$V2 <- gsub("\\|.*", "", simm$V2)

tmp <- abun$ASV
abun$ASV <- gsub("\\|.*", "", abun$ASV)
taxo <- data.frame(asv=abun$ASV, taxo=tmp) 
rownames(abun)  <- abun$ASV; abun$ASV <- NULL

curated <- lulu(abun, simm)
