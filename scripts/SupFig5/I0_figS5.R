
library(lulu)
library(data.table)
library(vegan)
library(tidyr)
library(circlize)

system("mkdir data/illumina/lulu")
setwd("~/IntraGenomic-variability/data/illumina/lulu")
system("cp ../dada2/default_pipeline/out/ASVs_taxo.fasta .")
system("vsearch --usearch_global ASVs_taxo.fasta --db ASVs_taxo.fasta --self --id .84 --iddef 1 --userout match_list.txt -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10")

abun <- fread("ASVs_abun.tsv", sep="\t")
simm <- fread("match_list.txt", sep="\t")

simm$V1 <- gsub("\\|.*", "", simm$V1)
simm$V2 <- gsub("\\|.*", "", simm$V2)

tmp <- abun$ASV
abun$ASV <- gsub("\\|.*", "", abun$ASV)
taxo <- data.frame(parent_id=abun$ASV, taxo=tmp)
taxo$taxo <- gsub("asv\\d+\\|", "", taxo$taxo)
abun$ASV <- NULL
abun <- as.data.frame(abun)
row.names(abun)  <- taxo$parent_id

curated <- lulu(abun, simm)
system("rm -f lulu*")

out <- curated$otu_map
out$child <- rownames(out)

out <- merge(out, taxo, by="parent_id", all.x=TRUE)

# write.table(out, "table_lulu.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
# write.table(curated$curated_table, "ASVs_abun_lulu.tsv", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

out <- out %>% separate(taxo, c("domain", "supergroup", "phylum", "class", "order", "family", "genus", "species"), sep="\\|", remove=TRUE)

out$taxo <- fifelse(out$class=="Polycystinea", out$order, out$class)
out$taxo <- fifelse(out$taxo=="Bacillariophyta", "Diatoms", 
                      fifelse(out$taxo=="Dinophyceae", "Dinoflagellates", 
                              fifelse(out$taxo=="Labyrinthulomycetes", "Fungi", 
                                      fifelse(out$taxo=="Basidiomycota", "Fungi", 
                                              fifelse(grepl("Sphaerozoidae|Collosphaeridae|Collophidiidae", out$family), "Collodaria", out$taxo)))))


file <- subset(out, class=="Polycystinea" & spread > 2 & total >= summary(out$total)[[3]])
file <- subset(out, spread >= 2 & total >= 3)

file <- file %>% select(child, parent_id, total)
colnames(file) <- c("from", "to", "value")

file$value <- c(decostand(file$value, "log"))

file$to <- factor(file$to, levels=unique(file[order(file$value, decreasing=TRUE),]$to))
file$from <- factor(file$from, levels=unique(file[order(file$value, decreasing=TRUE),]$from))

# write.table(filec, "data/files/fig3.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

pdf(paste("../../../plots/FigS5/FigS5_lulu_circular.pdf", sep=""), width=11.69, height=8.27, paper='special')
circular_plot <- chordDiagram(file, 
                              transparency = 0.5, 
                              order=union(unique(file[order(file$value, decreasing=TRUE),]$from), unique(file$to)))
dev.off()
