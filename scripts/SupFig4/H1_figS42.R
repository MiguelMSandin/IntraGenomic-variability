#----
#---- Loading packages  ----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)

library(ggplot2)
library(treemapify)
library(colorRamps)

#----
#---- Set working directory ------------------------------------------------------------------------
setwd("~/Desktop/PhD/PhD_Data/_SingleCell/data/")
#----
#---- Opening Illumina file ------------------------------------------------------------------------

file <- fread("illumina/dada2/200508/vsearch/ASVs_200508_abun_checked.tsv")

# Normalize data pero sample, to be able to compare each sample.
abun <- select(file, grep("Mge|Vil", names(file), value=TRUE))
abun <- decostand(abun, "total", 2)

# Extract only interested columns for taxonomy and further cleaning to merge with the transformed abundances data
taxo <- select(file, c("query", "AcNu", "Class", "Order", "Species", "id", "tpres", "tabun"))
taxo$taxo <- fifelse(taxo$Class=="Polycystinea", taxo$Order, taxo$Class)
taxo$taxo <- with(taxo, fifelse(taxo=="Spumellarida", "Spumellaria", 
                                fifelse(taxo=="Bacillariophyta", "Diatoms", 
                                        fifelse(taxo=="Dinophyceae", "Dinoflagellates", 
                                                fifelse(taxo=="Labyrinthulomycetes", "fungi", 
                                                        fifelse(taxo=="Basidiomycota", "fungi", taxo))))))

# Saving the median value of the reads for future filtering 
summary(taxo$tabun)
tmp <- median(taxo$tabun)

# Merging files 
file <- cbind(taxo, abun); rm(taxo, abun)

# Preparing the table to be plotted
file <- melt(file, id.vars=grep("Mge17|Vil", names(file), value=TRUE, invert=TRUE))

# Subsetting 
file <- file[value>0 & tabun >= tmp & tpres >= 3]

file <- file %>% separate(variable, c("cell", "replicate"), sep="_", remove=TRUE)

# Summarising the table

file <- file %>% group_by(cell, taxo) %>% summarise(abundance=sum(value), presence=length(unique(query)))

file$query <- gsub("asv","ASV ", file$query)

file <- melt(file, id.vars=c("cell", "taxo"))

file$cell <- factor(file$cell, levels=c("Mge17-9", "Mge17-124", "Vil490", "Vil496", "Mge17-81", "Mge17-82", "Vil480", "Vil497"))
file$variable <- factor(file$variable, levels=c("presence", "abundance"))

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

treemap <- ggplot(file, aes(area=value, fill=taxo))+
  geom_treemap()+
  facet_grid(cell~variable)+
  #geom_treemap_text(aes(label=paste(taxo, value, sep="\n")), place = "centre")+
  geom_treemap_text(aes(label=ifelse(variable=="abundance", paste(taxo), paste(taxo, value, sep="\n"))), place = "centre")+
  scale_fill_manual(values=color)+
  theme(legend.position="none")
treemap

pdf(paste("../plots/FigS4_Treemaps_clean_AbunASVs.pdf", sep=""), width=8.27, height=11.69, paper='special')
plot(treemap)
dev.off()


#---- 