#----
#---- Loading packages  ----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)

library(seqinr)
library(ape)

library(ggplot2)
library(forcats)
library(cowplot)

#----
#---- Set working directory ------------------------------------------------------------------------
setwd("~/IntraGenomic-variability/")
#---- 
#---- Preparing the data files for Sanger ----------------------------------------------------------

rm(list=ls())

# Abundances _______________________________________________________________________________________
sabun <- data.frame(query=labels(read.FASTA("data/sanger/V4/polycystines_V4_derep.fasta", type="DNA")))
sabun <- separate(sabun, query, c("query", "abun"), sep=";")
sabun$file <- "raw"
sabun$abun <- gsub("size=", "", sabun$abun)

# Vil325 and Vil496 share the exact same V4 (Vil325-BC01-A2 = 101): 
# Checking manually corresponds to 54 reads correspond to Vil325 and 47 to Vil496
sabun$abun[grepl("Vil325_BC01_A2", sabun$query)&grepl("raw", sabun$file)] <- 54
sabun <- rbind(sabun, c("Vil496_BC_manualRead", 47, "raw"))

dir <-"data/sanger/V4/swarm/"
for(i in grep("fasta", dir(dir), value=TRUE)){
  tmp <- data.frame(query=labels(read.FASTA(paste0(dir, i), type="DNA")))
  tmp$query <- gsub(";$", "", tmp$query)
  tmp <- separate(tmp, query, c("query", "abun"), sep=";")
  tmp$file <- gsub(".*swarm_", "", i) %>% gsub("\\.fasta*", "", .)
  tmp$abun <- gsub("size=", "", tmp$abun)
  sabun <- rbind(sabun, tmp)
}; rm(i, tmp)

sabun$cell <- gsub("_BC.*", "", sabun$query)

sabun$file <- factor(sabun$file, levels=c("raw", "d1", "d2", "d3"))
sabun$cell <- ifelse(sabun$file!="raw" & sabun$cell=="Vil325", "Vil325_&_Vil496", as.character(sabun$cell))
sabun$cell <- factor(sabun$cell, levels=c("Vil325", "Vil496", "Vil325_&_Vil496", "Mge17-81", "Mge17-82"))
sabun$abun <- as.numeric(sabun$abun)

# Now ordering the amplicons by file, cell and abundance

setorderv(sabun, cols=c("file", "cell", "abun"), c(1, 1, -1))
sabun$order <- 1:nrow(sabun)


# Similarities _____________________________________________________________________________________

ssim <- fread("data/sanger/V4/polycystines_V4_similarities.tsv")
colnames(ssim) <- c("query", "target", "id")
ssim$file <- "raw"

dir <- "data/sanger/V4/swarm/"
for(i in grep("similarities", dir(dir), value=TRUE)){
  tmp <- fread(paste0(dir, i), sep="\t")
  colnames(tmp) <- c("query", "target", "id")
  tmp$file <- gsub(".*swarm_", "", i) %>% gsub("\\.similarities*", "", .)
  ssim <- rbind(ssim, tmp)
}; rm(i, tmp)

ssim$cell <- gsub("_BC.*", "", ssim$query)
ssim$cellt <- gsub("_BC.*", "", ssim$target)
  
ssim <- ssim[ssim$cell == ssim$cellt]

ssim$file <- factor(ssim$file, levels=c("raw", "d1", "d2", "d3"))
ssim$cell <- factor(ssim$cell, levels=c("Vil325", "Vil496", "Mge17-81", "Mge17-82"))


#----
#---- Figure 4: Intracellular diversity Sanger -----------------------------------------------------

(plotSabun <- ggplot(sabun, aes(x=order, y=abun, colour=cell))+
  geom_point()+
  facet_wrap(~file, nrow=1, scales="free_x")+
  scale_colour_manual(values=c("springgreen1", "springgreen3", "springgreen2", "steelblue1", "steelblue3"))+
  theme_bw()+
  theme(legend.position="none"))


(plotSsim <- ggplot(ssim, aes(x=cell, y=id, fill=cell))+
  geom_jitter(color="grey80")+
  geom_boxplot(alpha=0.8)+
  facet_wrap(~file, nrow=1, scales="free_x")+
  scale_fill_manual(values=c("springgreen1", "springgreen3", "steelblue1", "steelblue3"))+
  theme_bw()+
  theme(legend.position="none"))

write.table(sabun, "data/fig4_abundance.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
write.table(ssim, "data/fig4_similarity.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

(Sanger <- plot_grid(plotSabun, plotSsim, labels=c("Abundance", "Intracellular simmilarity"), ncol=1, nrow=2))

pdf("plots/Figure4_Intracell_var.pdf", width=11.69, height=8.27, paper='special')
plot(Sanger)
dev.off()


#----
#----
