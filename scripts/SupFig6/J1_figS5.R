
library(seqinr)
library(ggplot2)
library(ggseqlogo)

setwd("/home/mmsandin/Desktop/PhD/0_Thesis/3.1_Chapter/github/data/sanger/alignments/")

data <-  read.fasta("raw_concatenated_polycystines_align_clean_mge17-82.fasta", seqtype= "AA", as.string = T)
data <- data.frame(name=paste(getAnnot(data)), sequence=paste(data, sep=""))

positions <- list(c(607, 626), c(674, 688), c(1007, 1034), c(1526, 1542),  c(1630, 1644))

seqs <- c()
for(j in positions){
    for(i in 1:nrow(data)){
        seqs[i] <- toupper(substr(data[i,2], unlist(j)[1], unlist(j)[2]))
        seqs[i] <- ifelse(grepl("^-*$", seqs[i]), seqs[i], gsub("-", "x", seqs[i]))
        seqs[i] <- gsub("T", "U", seqs[i])
    }
    (logo <- ggplot() + 
            geom_logo(seqs, method='p', namespace='ACUGx') + theme_logo() +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                  axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()))
    pdf(paste0("logo_", unlist(j)[1], "-", unlist(j)[2], ".pdf"), width=1.25, height=0.4, paper='special')
    plot(logo)
    dev.off()
}; rm(j, i)
