#----
#---- Loading packages  ----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)

library(vegan)

library(circlize)

#----
#---- Set working directory ------------------------------------------------------------------------
setwd("~/IntraGenomic-variability/")
#----
#---- Figure 3: CirclesPlots -----------------------------------------------------------------------

file <- fread("data/illumina/dada2/default_pipeline/out/ASVs_final_tax_abun.tsv")

(med <- summary(file$tabun)[[3]])

file <- file[Class=="Polycystinea" & tpres>2 & tabun > med]

file <- melt(file, id.vars=grep("Mge17|Vil", names(file), value=TRUE, invert=TRUE))

file <- select(file, c("ASV", "variable", "value"))
colnames(file) <- c("from", "to", "value")

file$value <- c(decostand(file$value, "log"))


file$to <- factor(file$to, levels=c("Mge17-124_M8R25", "Mge17-124_M8R26", "Mge17-124_M8R41", 
                                        "Mge17-81_M8R10", "Mge17-81_M8R12", "Mge17-81_M8R9", 
                                        "Mge17-82_M8R15", "Mge17-82_M8R16", "Mge17-82_M8R17", 
                                        "Mge17-9_M8R36", "Mge17-9_M8R38", "Mge17-9_M8R40", 
                                        "Vil480_M8R31", "Vil480_M8R32", "Vil480_M8R35", 
                                        "Vil490_M8R19", "Vil490_M8R20", "Vil490_M8R23", 
                                        "Vil496_M8R6", "Vil496_M8R7", "Vil496_M8R8", 
                                        "Vil497_M8R28", "Vil497_M8R29", "Vil497_M8R30"))

tmp <- c("Mge17-124_M8R25"="springgreen3", "Mge17-124_M8R26"="springgreen3", "Mge17-124_M8R41"="springgreen3", 
         "Mge17-81_M8R10"="steelblue3", "Mge17-81_M8R12"="steelblue3", "Mge17-81_M8R9"="steelblue3", 
         "Mge17-82_M8R15"="steelblue3", "Mge17-82_M8R16"="steelblue3", "Mge17-82_M8R17"="steelblue3", 
         "Mge17-9_M8R36"="springgreen3", "Mge17-9_M8R38"="springgreen3", "Mge17-9_M8R40"="springgreen3", 
         "Vil480_M8R31"="steelblue3", "Vil480_M8R32"="steelblue3", "Vil480_M8R35"="steelblue3", 
         "Vil490_M8R19"="springgreen3", "Vil490_M8R20"="springgreen3", "Vil490_M8R23"="springgreen3", 
         "Vil496_M8R6"="springgreen3", "Vil496_M8R7"="springgreen3", "Vil496_M8R8"="springgreen3", 
         "Vil497_M8R28"="steelblue3", "Vil497_M8R29"="steelblue3", "Vil497_M8R30"="steelblue3")


pdf(paste("plots/Fig3_CircularPlot.pdf", sep=""), width=11.69, height=8.27, paper='special')
circular_plot <- chordDiagram(file, 
                              transparency = 0.5, 
                              order=union(unique(file$from), unique(file$to)), 
                              grid.col=tmp)
dev.off()

