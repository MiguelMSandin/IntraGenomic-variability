#----
#---- Loading packages  ----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)

library(seqinr)
library(vegan)

library(ggplot2)

#----
#---- Set working directory ------------------------------------------------------------------------
setwd("~/IntraGenomic-variability/")
#----
#---- Setting variables and factors ----------------------------------------------------------------

entropy <- c()

dir <- "data/entropy/V4/"
for(i in grep(".*_positionInfo\\.tsv", dir(dir), value=TRUE)){
  tmp <- fread(paste0(dir, i))
  tmp$file <- gsub("_align_positionInfo\\.tsv", "", i)
  entropy <- rbind(entropy, tmp)
}; rm(i, tmp)

unique(entropy$file)
entropy$method <- fifelse(grepl("ref", entropy$file), "Reference", 
                        fifelse(grepl("sanger", entropy$file), "Sanger", 
                                fifelse(grepl("minion", entropy$file), "MinION", 
                                        fifelse(grepl("illumina", entropy$file), "Illumina", "NA"))))
entropy$method <- factor(entropy$method, levels=c("Reference", "Sanger", "Illumina", "MinION"))
unique(entropy$method)

entropy$taxo <- fifelse(grepl("nass", entropy$file), "Nassellaria", 
                          fifelse(grepl("spum", entropy$file), "Spumellaria", "NA"))
entropy$taxo <- factor(entropy$taxo, levels=c("Nassellaria", "Spumellaria"))
unique(entropy$taxo)

entropy$cell <- fifelse(grepl("Vil325", entropy$file), "Vil325", 
                        fifelse(grepl("Vil496", entropy$file), "Vil496", 
                                fifelse(grepl("Mge17-81", entropy$file), "Mge17-81", 
                                        fifelse(grepl("Mge17-82", entropy$file), "Mge17-82", as.character(entropy$method)))))

entropy$color <- paste(entropy$taxo, entropy$method, entropy$cell, sep="_")
entropy$color <- entropy$color %>% gsub("_Illumina$|_MinION$|_Reference$", "", .)
entropy$color <- factor(entropy$color, levels=c("Nassellaria_Reference", 
                                                "Nassellaria_Sanger_Vil325", "Nassellaria_Sanger_Vil496", 
                                                "Nassellaria_MinION", "Nassellaria_Illumina", 
                                                "Spumellaria_Reference", 
                                                "Spumellaria_Sanger_Mge17-81", "Spumellaria_Sanger_Mge17-82",
                                                "Spumellaria_Illumina", 
                                                "Spumellaria_MinION"))
unique(entropy$color)

write.table(entropy, "data/files/fig5.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

#----
#---- Plotting -------------------------------------------------------------------------------------

(dive <- ggplot(entropy, aes(x=position, y=shannon_clean, colour=color))+
  geom_point(aes())+
  geom_smooth(level=.99)+
  facet_grid(rows=vars(method), cols=vars(taxo), scales="free_x")+
  scale_colour_manual(values=c("grey60", "springgreen1", "springgreen3", "springgreen4", "springgreen4",
                               "grey60", "steelblue1", "steelblue3", "steelblue4", "steelblue4"))+
  theme_bw()+
  theme(legend.position="none"))

# Exporting the PDF of the plot
pdf("plots/Fig5.pdf", width=11.69, height=8.27, paper='special')
plot(dive)
dev.off()
# But this graph is too heavy too edit in inkscape, so let's remove all 0 points (without affecting geom_smoot !!)

(divep <- ggplot()+
  geom_point(subset(entropy, shannon_clean>0), mapping=aes(x=position, y=shannon_clean, colour=color))+
  geom_smooth(entropy, mapping=aes(x=position, y=shannon_clean, colour=color), level=.95)+
  facet_grid(rows=vars(method), cols=vars(taxo), scales="free_x")+
  scale_colour_manual(values=c("grey60", "springgreen1", "springgreen3", "springgreen4", "springgreen4",
                               "grey60", "steelblue1", "steelblue3", "steelblue4", "steelblue4"))+
  theme_bw()+
  theme(legend.position="none"))

# Exporting the PDF of the plot
pdf("plots/Fig5_clean.pdf", width=11.69, height=8.27, paper='special')
plot(divep)
dev.off()

#----
