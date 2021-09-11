#----
#---- Loading packages  ----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)

library(seqinr)

library(ggplot2)
library(treemapify)
library(colorRamps)

#----
#---- Set working directory ------------------------------------------------------------------------
setwd("~/IntraGenomic-variability/")
#----
#---- Opening Sanger files -------------------------------------------------------------------------

# Sanger ___________________________________________________________________________________________
sanger <- fread("data/sanger/blast/raw_concatenated_blast.tsv", sep="\t")
colnames(sanger) <- c("qseqid", "sseqid", "sacc", "stitle", "sscinames", "staxids", "sskingdoms", "sblastnames", "pident", "slen", "lengthAlign", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extract the first blast match
sanger1st <- data.frame()
for(i in unique(sanger$qseqid)){
  p <- grep(paste0("^", i, "$"), unique(sanger$qseqid))
  l <- length(unique(sanger$qseqid))
  if((p %% (round(l*0.1,0))) == 0){message("  ", round(p/l*100, 0), "%")}
  df <- subset(sanger, qseqid==i)[1,]
  if(grepl("U|uncultured", df$sscinames)){
    group <- unique(subset(sanger, qseqid==i)$sscinames)
    df$group <- group[!grepl("U|uncultured", group)][1]
  }else{
    df$group <- df$sscinames
  }
  sanger1st <- rbind(sanger1st, df)
}; rm(i, p, l, df, group)

sanger1st <- sanger1st %>% separate(qseqid, c("cell", "replicate", "read"), sep="_", remove=FALSE)

# Annotating the different taxonomic groups 
unique(sanger1st$group)
sanger1st$taxo <- fifelse(grepl("Eucyrtidium", sanger1st$group), "Nassellaria", 
                          fifelse(grepl("Triastrum|Heliodiscus|Styptosphaera|Spongodiscidae|Euchitonia|Rhizosphaera|Dictyocoryne", sanger1st$group), "Spumellaria", 
                                  fifelse(sanger1st$group=="Oikopleura dioica","Tunicates", 
                                          fifelse(sanger1st$group=="Oikopleura longicauda","Tunicates", 
                                                  fifelse(sanger1st$group=="Cyclotrichium cyclokaryon","Ciliophora", 
                                                          fifelse(sanger1st$group=="Paraphysomonas imperforata","Chrysophyceae", 
                                                                  fifelse(sanger1st$group=="Chaetoceros curvisetus","Diatoms", 
                                                                          fifelse(sanger1st$group=="Cryptosporidium canis","Apicomplexa", 
                                                                                  fifelse(sanger1st$group=="Oblongichytrium sp. HK9","Fungi", 
                                                                                          fifelse(sanger1st$group=="Schizochytrium minutum","Fungi", sanger1st$group))))))))))
unique(sanger1st$taxo)

# "Mge17-82_BC11_H7" couldn't retrieve any name
# Exploring in the full blast file for other matches
tmp <- sanger1st[is.na(sanger1st$taxo),]$qseqid
tmp <- subset(sanger, qseqid==tmp)
unique(tmp$sscinames)
unique(tmp$sblastnames)
# Other matches are related to Syndiniales, Dinoflagellates and alveolates
# Therefore assigning manually to Dinoflagellates
sanger1st$taxo[grepl("Mge17-82_BC11_H7", sanger1st$qseqid)] <- "Dinoflagellates"

table(sanger1st$taxo)

# Merging with the sequences 
tmp <-  seqinr::read.fasta("data/sanger/raw/raw_concatenated.fasta", seqtype= "AA", as.string = T)
tmp <- data.frame(qseqid=gsub(">", "", paste(getAnnot(tmp))), sequence=paste0(tmp))

all(tmp$qseqid %in% sanger1st$qseqid)

sanger1st <- merge(sanger1st, tmp, by="qseqid"); rm(tmp)

# And export the table
write.table(sanger1st, "data/files/fig2_sanger.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

#___________________________________________________________________________________________________
#_________ Remember to check manually the taxonomic assignation in (e.g.) libreOffice Calc _________
#___________________________________________________________________________________________________

#---- Opening MinION files -------------------------------------------------------------------------

# MinION ___________________________________________________________________________________________
minion <- fread("data/minion/blast/raw_blast.tsv", sep="\t")
colnames(minion) <- c("qseqid", "sseqid", "sacc", "stitle", "sscinames", "staxids", "sskingdoms", "sblastnames", "pident", "slen", "lengthAlign", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Take only demultiplexed reads
demultiplexed <- seqinr::read.fasta("data/minion/raw/demultiplexed.fasta", seqtype= "AA", as.string = T)
demultiplexed <- data.frame(name=gsub(">", "", paste(getAnnot(demultiplexed))), sequence=paste0(demultiplexed))
demultiplexed$name <- gsub(" .*", "", demultiplexed$name)
demultiplexed <- demultiplexed %>% separate(name, c("cell", "replicate", "qseqid"), sep="_", remove=TRUE)

minion <- minion[minion$qseqid %in% demultiplexed$qseqid]

# Extract the first blast match
minion1st <- data.frame()
for(i in unique(minion$qseqid)){
  p <- grep(paste0("^", i, "$"), unique(minion$qseqid))
  l <- length(unique(minion$qseqid))
  if((p %% (round(l*0.1,0))) == 0){message("  ", round(p/l*100, 0), "%")}
  df <- subset(minion, qseqid==i)[1,]
  if(grepl("U|uncultured", df$sscinames)){
    group <- unique(subset(minion, qseqid==i)$sscinames)
    df$group <- group[!grepl("U|uncultured", group)][1]
  }else{
    df$group <- df$sscinames
  }
  minion1st <- rbind(minion1st, df)
}; rm(i, p, l, df, group)

minion1st <- merge(minion1st, demultiplexed, by="qseqid", all=TRUE, sort=FALSE)

# Annotating the different taxonomic groups 
unique(minion1st$group)
minion1st$taxo <- fifelse(grepl("Eucyrtidium", minion1st$group), "Nassellaria", 
                          fifelse(grepl("Triastrum|Styptosphaera|Euchitonia|Rhizosphaera", minion1st$group), "Spumellaria", (
                            fifelse(minion1st$group=="Aureobasidium pullulans","Fungi", 
                                    fifelse(minion1st$group=="Gonyaulax elongata","Dinoflagellates", 
                                            fifelse(minion1st$group=="Paraphysomonas imperforata","Chrysophyceae", 
                                                    fifelse(minion1st$group=="Chaetoceros curvisetus","Diatoms", 
                                                            fifelse(minion1st$group=="Eimeria gaimardi","Alveolates", 
                                                                    fifelse(minion1st$group=="Alexandrium margalefii","Dinoflagellates", 
                                                                            fifelse(minion1st$group=="Hepatozoon tenuis","Fungi", 
                                                                                    fifelse(minion1st$group=="Ceratocorys sp. NY002","Dinoflagellates", minion1st$group)))))))))))
table(minion1st$taxo)

# Write table
write.table(minion1st, "data/files/fig2_minion.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

#___________________________________________________________________________________________________
#_________ Remember to check manually the taxonomic assignation in (e.g.) libreOffice Calc _________
#___________________________________________________________________________________________________

#---- Opening Illumina files -----------------------------------------------------------------------

illumina <- fread("data/illumina/dada2/default_pipeline/out/ASVs_final_tax_abun.tsv")

illumina <- melt(illumina, id.vars=grep("Mge17|Vil", names(illumina), value=TRUE, invert=TRUE))

illumina <- illumina %>% separate(variable, c("cell", "replicate"), sep="_", remove=TRUE)

illuOut <- select(illumina, c("cell", "replicate", "ASV", "sequence", "Supergroup", "Division", "Class", "Order", "Family", "Species", "min_BS", "tpres", "tabun", "value"))

illuOut <- illuOut[value>0 & !is.na(Order)]

illuOut$taxo <- fifelse(illuOut$Class=="Polycystinea", illuOut$Order, illuOut$Class)
illuOut$taxo <- fifelse(illuOut$taxo=="Bacillariophyta", "Diatoms", 
                        fifelse(illuOut$taxo=="Dinophyceae", "Dinoflagellates", 
                                fifelse(illuOut$taxo=="Labyrinthulomycetes", "Fungi", 
                                        fifelse(illuOut$taxo=="Basidiomycota", "Fungi", 
                                                fifelse(grepl("Sphaerozoidae|Collosphaeridae|Collophidiidae", illuOut$Family), "Collodaria", illuOut$taxo)))))


sort(table(illuOut$taxo), decreasing=TRUE)

write.table(illuOut, "data/files/fig2_illumina.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)


#----
#---- Creating a single file with all the information for plotting ---------------------------------

tmps <- fread("data/files/fig2_sanger.tsv", sep="\t")
tmpm <- fread("data/files/fig2_minion.tsv", sep="\t")
tmpi <- fread("data/files/fig2_illumina.tsv", sep="\t")

tmps <- select(tmps, c("cell", "replicate", "qseqid", "sequence", "taxo")); tmps$method <- "Sanger"
tmpm <- select(tmpm, c("cell", "replicate", "qseqid", "sequence", "taxo")); tmpm$method <- "MinION"
tmpi <- select(tmpi, c("cell", "replicate", "ASV", "sequence", "taxo")); tmpi$method <- "Illumina"
colnames(tmpi) <- c("cell", "replicate", "qseqid", "sequence", "taxo", "method")

file <- rbind(tmps, tmpm, tmpi); rm(tmps, tmpm, tmpi)

file <- file %>% group_by(method, cell, taxo) %>% summarise(reads=length(unique(sequence)))

file$method <- factor(file$method, levels=c("Sanger", "MinION", "Illumina"))

file$cell <- factor(file$cell, levels=c("Mge17-9", "Mge17-124", "Vil325", "Vil490", "Vil496", "Mge17-81", "Mge17-82", "Vil480", "Vil497"))

file$species <- fifelse(file$cell=="Mge17-9", "Extotoxon undulatum", 
                        fifelse(file$cell=="Mge17-124",  "Carpocanium obliqua", 
                                fifelse(file$cell=="Vil325",   "Eucyrtidium cienkowski", 
                                        fifelse(file$cell=="Vil490",   "Pterocorys zanclea", 
                                                fifelse(file$cell=="Vil496",   "Eucyrtidium acuminatum", 
                                                        fifelse(file$cell=="Mge17-81",  "Rhizosphaera trigonacantha",  
                                                                fifelse(file$cell=="Mge17-82",  "Spongosphaera streptacantha",  
                                                                        fifelse(file$cell=="Vil480",   "Tetrapyle octacantha", 
                                                                                fifelse(file$cell=="Vil497", "Arachnospongus varians", "NA")))))))))

file$species <- factor(file$species, levels=c("Extotoxon undulatum", "Carpocanium obliqua", "Eucyrtidium cienkowski", "Pterocorys zanclea", "Eucyrtidium acuminatum",
                                              "Rhizosphaera trigonacantha", "Spongosphaera streptacantha", "Tetrapyle octacantha", "Arachnospongus varians"))

file <- file[!is.na(file$taxo),]

write.table(file, "data/files/fig2.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

#----
#---- Figure 2: TreeMaps ---------------------------------------------------------------------------

library(colorspace)

color <- as.character(file$taxo)
for(i in unique(color)){
  #color[which(as.character(file$taxo)==i)]=primary.colors(length(unique(color)))[grep(i, unique(color))]
  color[which(as.character(file$taxo)==i)]=qualitative_hcl(length(unique(color)))[grep(i, unique(color))]
}; rm(i)
color[which(as.character(file$taxo)=="Nassellaria")]="springgreen3"
color[which(as.character(file$taxo)=="Acantharea")]="yellow3"
color[which(as.character(file$taxo)=="Spumellaria")]="steelblue3"
color[which(as.character(file$taxo)=="Collodaria")]="orangered3"

tmp <- unique(data.frame(taxo=as.character(file$taxo), color=color))
color <- as.character(tmp$color)
names(color) <- as.character(tmp$taxo); rm(tmp)

(treemap <- ggplot(file, aes(area=reads, fill=taxo))+
  geom_treemap()+
  facet_grid(species~method, space="free")+
  geom_treemap_text(aes(label= paste(taxo, reads, sep="\n")), place = "centre")+
  scale_fill_manual(values=color)+
  theme(legend.position="none"))

pdf(paste("plots/Fig2.pdf", sep=""), width=8.27, height=11.69, paper='special')
plot(treemap)
dev.off()


#---- 
#---- Exploring average similarity of illumina reads  ----

sim <- fread("data/illumina/dada2/default_pipeline/out/ASVs_taxo_id.tsv", sep="\t")
colnames(sim) <- c("query", "target", "id", "alnlen", "mism", "opens", "qlo", "qhi", "tlo", "thi", "evalue", "bits")

summary(sim$id)
sd(sim$id)

#---- 
  