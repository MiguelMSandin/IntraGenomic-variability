#----
#---- Loading packages  ----------------------------------------------------------------------------

library(data.table)
library(ggplot2)

#---- Set working directory ------------------------------------------------------------------------
setwd("~/IntraGenomic-variability/data/")
#----   
#----  Sanger --------------------------------------------------------------------------------------
dist <- fread("sanger/consensus/raw_concatenated_polycystines_align.dist")
colnames(dist) <- c("query", "target", "dist")

dist$qcell <- gsub("_.*", "", dist$query)
dist$tcell <- gsub("_.*", "", dist$target)

dist$sim <- fifelse(dist$qcell==dist$tcell, "Intra-cell", "Between-cells")

summary(dist$dist)
hist(dist$dist)

ggplot(dist, aes(x=sim, y=dist))+
    geom_jitter(color="grey80")+
    geom_boxplot(alpha=0.6)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=60, hjust=0.9, vjust=0.9))+
    theme(legend.position="none")

ggplot(dist, aes(x=qcell, y=dist))+
    geom_jitter(color="grey80")+
    geom_boxplot(alpha=0.6)+
    facet_grid(~sim, scales="free")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=60, hjust=0.9, vjust=0.9))+
    theme(legend.position="none")

# Exploring only intracellular differences

intra <- subset(dist, sim=="Intra-cell" & dist!=1)

summary(intra$dist)
hist(intra$dist)

ggplot(intra, aes(x=qcell, y=dist))+
    geom_jitter(color="grey80")+
    geom_hline(yintercept=0.03, color="orangered1")+
    geom_hline(yintercept=0.05, color="orangered2")+
    geom_boxplot(alpha=0.6)+
    theme_bw()+
    theme(axis.text.x=element_text(angle=60, hjust=0.9, vjust=0.9))+
    theme(legend.position="none")

for(cell in unique(intra$qcell)){
    cat("For", cell, "the summary statistics are:\n")
    print(summary(subset(intra, qcell==cell)$dist))
    cat("\n")
}

#----   
#----  MinION --------------------------------------------------------------------------------------

dist <- fread("minion/consensus/demultiplexed_polycystines_align.dist")
colnames(dist) <- c("query", "target", "dist")

dist$qcell <- gsub("_.*", "", dist$query)
dist$tcell <- gsub("_.*", "", dist$target)

dist$sim <- fifelse(dist$qcell==dist$tcell, "Intra-cell", "Between-cells")

summary(dist$dist)
hist(dist$dist)

ggplot(dist, aes(x=sim, y=dist))+
    geom_jitter(color="grey80")+
    geom_boxplot(alpha=0.6)+
    geom_hline(yintercept=0.30, color="orangered1")+
    geom_hline(yintercept=0.25, color="orangered2")+
    geom_hline(yintercept=0.03, color="orangered2")+
    theme_bw()+
    theme(axis.text.x=element_text(angle=60, hjust=0.9, vjust=0.9))+
    theme(legend.position="none")

#----   