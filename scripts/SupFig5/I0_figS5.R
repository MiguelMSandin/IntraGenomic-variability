
library(lulu)
library(data.table)
library(vegan)
library(tidyr)
library(ggplot2)
library(treemapify)
library(colorRamps)
library(colorspace)


system("mkdir data/illumina/lulu")
setwd("data/illumina/lulu/")
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

out <- curated$otu_map
out$child <- rownames(out)

out <- merge(out, taxo, by="parent_id", all.x=TRUE)

write.table(out, "table_lulu.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

abunl <- curated$curated_table
abunl <- decostand(abunl, "total", 2)
abunl$parent_id <- rownames(abunl)

abunl <- merge(abunl, taxo, by="parent_id")
abunl <- abunl %>% separate(taxo, c("domain", "supergroup", "phylum", "class", "order", "family", "genus", "species"), sep="\\|", remove=TRUE)

abunl$taxo <- fifelse(abunl$class=="Polycystinea", abunl$order, abunl$class)
abunl$taxo <- fifelse(abunl$taxo=="Bacillariophyta", "Diatoms", 
                      fifelse(abunl$taxo=="Dinophyceae", "Dinoflagellates", 
                              fifelse(abunl$taxo=="Labyrinthulomycetes", "Fungi", 
                                      fifelse(abunl$taxo=="Basidiomycota", "Fungi", 
                                              fifelse(grepl("Sphaerozoidae|Collosphaeridae|Collophidiidae", abunl$family), "Collodaria", abunl$taxo)))))

file <- melt(as.data.table(abunl), id.vars=grep("Mge17|Vil", names(abunl), value=TRUE, invert=TRUE))

file <- file[value>0 & !is.na(order)]

file <- file %>% separate(variable, c("cell", "replicate"), sep="_", remove=TRUE)
file$cell <- factor(file$cell, levels=c("Mge17-9", "Mge17-81", "Mge17-124", "Mge17-82", "Vil490", "Vil480", "Vil496", "Vil497"))

file$specimen <- fifelse(file$cell=="Mge17-9", "Extotoxon undulatum", 
                         fifelse(file$cell=="Mge17-124",  "Carpocanium obliqua", 
                                 fifelse(file$cell=="Vil490",   "Pterocorys zanclea", 
                                         fifelse(file$cell=="Vil496",   "Eucyrtidium acuminatum", 
                                                 fifelse(file$cell=="Mge17-81",  "Rhizosphaera trigonacantha",  
                                                         fifelse(file$cell=="Mge17-82",  "Spongosphaera streptacantha",  
                                                                 fifelse(file$cell=="Vil480",   "Tetrapyle octacantha", 
                                                                         fifelse(file$cell=="Vil497", "Arachnospongus varians", "NA"))))))))
file$specimen <- factor(file$specimen, levels=c("Extotoxon undulatum", "Rhizosphaera trigonacantha", "Carpocanium obliqua", "Spongosphaera streptacantha", 
                                                "Pterocorys zanclea", "Tetrapyle octacantha", "Eucyrtidium acuminatum", "Arachnospongus varians"))

file <- file %>% group_by(cell, specimen, parent_id, taxo) %>% summarise(value=sum(value))

color <- as.character(file$taxo)
for(i in unique(color)){
    color[which(as.character(file$taxo)==i)]=qualitative_hcl(length(unique(color)))[grep(i, unique(color))]
}; rm(i)
color[which(as.character(file$taxo)=="Nassellaria")]="springgreen3"
color[which(as.character(file$taxo)=="Acantharea")]="yellow3"
color[which(as.character(file$taxo)=="Spumellaria")]="steelblue3"
color[which(as.character(file$taxo)=="Collodaria")]="orangered3"

tmp <- unique(data.frame(taxo=as.character(file$taxo), color=color))
color <- as.character(tmp$color)
names(color) <- as.character(tmp$taxo); rm(tmp)

(treemap <- ggplot(file, aes(area=value, fill=taxo))+
        geom_treemap()+
        facet_wrap(~specimen, nrow=4)+
        geom_treemap_text(aes(label= paste(taxo, parent_id, sep="\n")), place = "centre")+
        scale_fill_manual(values=color)+
        theme(legend.position="none"))

pdf(paste("/home/mmsandin/Desktop/PhD/0_Thesis/3.1_Chapter/figures/figS5_lulu/FigS5_lulu.pdf", sep=""), width=8.27, height=11.69, paper='special')
plot(treemap)
dev.off()

