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

dir <- "data/entropy/"
for(i in grep("rDNA.*_positionInfo\\.tsv", dir(dir), value=TRUE)){
  tmp <- fread(paste0(dir, i))
  tmp$file <- gsub("_positionInfo\\.tsv", "", i)
  entropy <- rbind(entropy, tmp)
}; rm(i, tmp)

unique(entropy$file)
entropy$method <- fifelse(grepl("sanger", entropy$file), "Sanger", 
                          fifelse(grepl("minion", entropy$file), "MinION", "NA"))
entropy$method <- factor(entropy$method, levels=c("Sanger", "MinION"))
unique(entropy$method)

entropy$taxo <- fifelse(grepl("nass", entropy$file), "Nassellaria", 
                          fifelse(grepl("spum", entropy$file), "Spumellaria", "NA"))
entropy$taxo <- factor(entropy$taxo, levels=c("Nassellaria", "Spumellaria"))
unique(entropy$taxo)

entropy$cell <- fifelse(grepl("Vil325", entropy$file), "Vil325", 
                        fifelse(grepl("Vil496", entropy$file), "Vil496", 
                                fifelse(grepl("Mge17-81", entropy$file), "Mge17-81", 
                                        fifelse(grepl("Mge17-82", entropy$file), "Mge17-82", as.character(entropy$method)))))
unique(entropy$cell)

entropy$color <- paste(entropy$method, entropy$taxo, entropy$cell, sep="_")
entropy$color <- factor(entropy$color, levels=c("Sanger_Nassellaria_Vil325", "Sanger_Nassellaria_Vil496", 
                                                "Sanger_Spumellaria_Mge17-81", "Sanger_Spumellaria_Mge17-82",
                                                "MinION_Nassellaria_Vil325", "MinION_Nassellaria_Vil496", 
                                                "MinION_Spumellaria_Mge17-81", "MinION_Spumellaria_Mge17-82"))
unique(entropy$color)

write.table(entropy, "data/files/figS5.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

#----
#---- Plotting -------------------------------------------------------------------------------------

(dive <- ggplot(entropy, aes(x=position, y=shannon_clean, colour=color))+
   geom_point(aes())+
   geom_smooth(level=.99)+
   facet_wrap(method~taxo, ncol=1)+
   scale_colour_manual(values=c("springgreen1", "springgreen3", "steelblue1", "steelblue3", 
                                "springgreen1", "springgreen3", "steelblue1", "steelblue3"))+
   theme_bw()+
   theme(legend.position="none"))

# Exporting the PDF of the plot
pdf("plots/FigS6.pdf", width=11.69, height=8.27, paper='special')
plot(dive)
dev.off()
# But this graph is too heavy too edit in inkscape, so let's remove all 0 points (without affecting geom_smoot !!)

(divep <- ggplot()+
    geom_point(subset(entropy, shannon_clean>0), mapping=aes(x=position, y=shannon_clean, colour=color))+
    geom_smooth(entropy, mapping=aes(x=position, y=shannon_clean, colour=color), level=.95)+
    facet_wrap(method~taxo, ncol=1, scales="free_x")+
    scale_colour_manual(values=c("springgreen1", "springgreen3", "steelblue1", "steelblue3", 
                                 "springgreen1", "springgreen3", "steelblue1", "steelblue3"))+
    theme_bw()+
    theme(legend.position="none"))

# Exporting the PDF of the plot
pdf("plots/FigS6_clean.pdf", width=11.69, height=8.27, paper='special')
plot(divep)
dev.off()

#----
