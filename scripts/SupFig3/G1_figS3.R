#----
#---- Loading packages  ----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)

library(vegan)

library(ggplot2)
library(treemapify)
library(colorRamps)
library(colorspace)

#----
#---- Set working directory ------------------------------------------------------------------------
setwd("~/IntraGenomic-variability/")
#----
#---- Opening Illumina file ------------------------------------------------------------------------

file <- fread("data/illumina/dada2/default_pipeline/out/ASVs_final_tax_abun.tsv")

# Normalize data pero sample, to be able to compare each sample.
abun <- select(file, grep("Mge|Vil", names(file), value=TRUE))
abun <- decostand(abun, "total", 2)

# Extract only interested columns for taxonomy and further cleaning to merge with the transformed abundances data
taxo <- select(file, c("ASV", "sequence", "Class", "Order", "Family", "Species", "min_BS", "tpres", "tabun"))

taxo$taxo <- fifelse(taxo$Class=="Polycystinea", taxo$Order, taxo$Class)
taxo$taxo <- fifelse(taxo$taxo=="Bacillariophyta", "Diatoms", 
                        fifelse(taxo$taxo=="Dinophyceae", "Dinoflagellates", 
                                fifelse(taxo$taxo=="Labyrinthulomycetes", "Fungi", 
                                        fifelse(taxo$taxo=="Basidiomycota", "Fungi", 
                                                fifelse(grepl("Sphaerozoidae|Collosphaeridae|Collophidiidae", taxo$Family), "Collodaria", taxo$taxo)))))

# Merging files 
file <- cbind(taxo, abun); rm(taxo, abun)

# Preparing the table to be plotted
file <- melt(file, id.vars=grep("Mge17|Vil", names(file), value=TRUE, invert=TRUE))

file <- file[value>0 & !is.na(Order)]

file <- file %>% separate(variable, c("cell", "replicate"), sep="_", remove=TRUE)
file$cell <- factor(file$cell, levels=c("Mge17-9", "Mge17-81", "Mge17-124", "Mge17-82", "Vil490", "Vil480", "Vil496", "Vil497"))

file$species <- fifelse(file$cell=="Mge17-9", "Extotoxon undulatum", 
                        fifelse(file$cell=="Mge17-124",  "Carpocanium obliqua", 
                                        fifelse(file$cell=="Vil490",   "Pterocorys zanclea", 
                                                fifelse(file$cell=="Vil496",   "Eucyrtidium acuminatum", 
                                                        fifelse(file$cell=="Mge17-81",  "Rhizosphaera trigonacantha",  
                                                                fifelse(file$cell=="Mge17-82",  "Spongosphaera streptacantha",  
                                                                        fifelse(file$cell=="Vil480",   "Tetrapyle octacantha", 
                                                                                fifelse(file$cell=="Vil497", "Arachnospongus varians", "NA"))))))))
file$species <- factor(file$species, levels=c("Extotoxon undulatum", "Rhizosphaera trigonacantha", "Carpocanium obliqua", "Spongosphaera streptacantha", 
                                              "Pterocorys zanclea", "Tetrapyle octacantha", "Eucyrtidium acuminatum", "Arachnospongus varians"))

file <- file %>% group_by(cell, species, ASV, taxo) %>% summarise(value=sum(value))

write.table(file, "data/files/figS3.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

#----
#---- Figure 2: TreeMaps ---------------------------------------------------------------------------

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
  facet_wrap(~species, nrow=4)+
  geom_treemap_text(aes(label= paste(taxo, ASV, sep="\n")), place = "centre")+
  scale_fill_manual(values=color)+
  theme(legend.position="none"))

pdf(paste("plots/FigS3.pdf", sep=""), width=8.27, height=11.69, paper='special')
plot(treemap)
dev.off()


#---- 