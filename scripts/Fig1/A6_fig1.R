#----
#---- Loading packages  ----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)

library(ggplot2)

library(ape)
library(seqinr)

#----
#---- Set working directory ------------------------------------------------------------------------
setwd("~/IntraGenomic-variability/")
#----
#---- Open Sanger files ----------------------------------------------------------------------------

# Sanger ___________________________________________________________________________________________
sanger <- fread("data/sanger/blast/raw_blast.tsv", sep="\t")
colnames(sanger) <- c("qseqid", "sseqid", "sacc", "stitle", "sscinames", "staxids", "sskingdoms", "sblastnames", "pident", "slen", "lengthAlign", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Extract the first blast match
sanger1st <- data.frame()
for(i in unique(sanger$qseqid)){
  p <- grep(paste0("^", i, "$"), unique(sanger$qseqid))
  l <- length(unique(sanger$qseqid))
  if((p %% (round(l*0.1,0))) == 0){message("  ", round(p/l*100, 0), "%")}
  df <- subset(sanger, qseqid==i)[1,]
  sanger1st <- rbind(sanger1st, df)
}; rm(i, p, l, df)

sanger1st <- sanger1st %>% separate(qseqid, c("order", "primer", "cell", "replicate", "sequence"), sep="_", remove=FALSE)

# Getting length of the sequences 
tmp <-  seqinr::read.fasta("data/sanger/raw/raw.fasta", seqtype= "AA", as.string=TRUE)
tmp <- data.frame(qseqid=gsub(">", "", paste(seqinr::getAnnot(tmp))), sequence=paste0(tmp), length=nchar(as.character(paste0(tmp))))

all(tmp$qseqid %in% sanger1st$qseqid)

# And merging with the final dataset
sanger1st <- merge(sanger1st, tmp, by="qseqid", all=TRUE, sort=FALSE); rm(tmp)

#---- Open Minion files ----------------------------------------------------------------------------

# MinION ___________________________________________________________________________________________
minion <- fread("data/minion/blast/raw_blast.tsv", sep="\t")
colnames(minion) <- c("qseqid", "sseqid", "sacc", "stitle", "sscinames", "staxids", "sskingdoms", "sblastnames", "pident", "slen", "lengthAlign", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

demultiplexed <- data.frame()
for(i in grep("\\.fasta$", dir("data/minion/demultiplexed/"), value=TRUE)){
  if(file.size(paste0("data/minion/demultiplexed/", i))==0){
    demultiplexed <- rbind(demultiplexed, data.frame(file=gsub("\\.fasta$", "", i), qseqid="NA"))
  }else{
    lab <- labels(read.FASTA(paste0("data/minion/demultiplexed/", i), type="DNA"))
    lab <- gsub(" .*", "", lab)
    demultiplexed <- rbind(demultiplexed, data.frame(file=rep(gsub("\\.fasta$", "", i), length(lab)), qseqid=lab))
  }
}; (rm(i, lab))

#write.table(demultiplexed, "files/minion_demultiplexed.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

minion <- minion[minion$qseqid %in% demultiplexed$qseqid]

# Extract the first blast match
minion1st <- data.frame()
for(i in unique(minion$qseqid)){
  p <- grep(i, unique(minion$qseqid))
  l <- length(unique(minion$qseqid))
  if((p %% (round(l*0.1, 0))) == 0){message("  ", round(p/l*100, 0), "%")}
  df <- subset(minion, qseqid==i)[1,]
  minion1st <- rbind(minion1st, df)
}; rm(i, p, l, df)

minion1st <- merge(minion1st, demultiplexed, by="qseqid", all=TRUE, sort=FALSE)

minion1st <- minion1st %>% separate(file, c("cell", "replicate"), sep="_", remove=TRUE)


# Getting length of the sequences 
tmp <-  seqinr::read.fasta("data/minion/raw/raw.fasta", seqtype= "AA", as.string=TRUE)
tmp <- data.frame(qseqid=paste(seqinr::getAnnot(tmp)), sequence=paste0(tmp), length=nchar(as.character(paste0(tmp))))
tmp$qseqid <- gsub(">", "", tmp$qseqid) %>% gsub(" .*", "", .)

tmp <- tmp[tmp$qseqid %in% demultiplexed$qseqid,]

all(tmp$qseqid %in% minion1st$qseqid)

# And merging with the final dataset
minion1st <- merge(minion1st, tmp, by="qseqid", all=TRUE, sort=FALSE); rm(tmp, demultiplexed)

#---- Open Illumina files --------------------------------------------------------------------------

# Illumina _________________________________________________________________________________________

illumina <- data.frame()
for(i in grep("_id\\.tsv$", dir("data/illumina/merged/"), value=TRUE)){
  p <- grep(i, grep("_id\\.tsv$", dir("data/illumina/merged/"), value=TRUE))
  l <- length(grep("_id\\.tsv$", dir("data/illumina/merged/"), value=TRUE))
  if((p %% (round(l*0.1,0))) == 0){message("  ", round(p/l*100, 0), "%")}

  df <- fread(paste0("data/illumina/merged/", i), sep="\t")
  df$cell <- gsub("_.*", "", i)
  df$replicate <- gsub("_merged_id.tsv", "", i) %>% gsub(".*_", "", .)
  
  dff <- seqinr::read.fasta(paste0("data/illumina/merged/", gsub("_id\\.tsv", "\\.fasta", i)), seqtype= "AA", as.string=TRUE)
  dff <- data.frame(V1=paste(seqinr::getAnnot(dff)), sequence=paste0(dff), length=nchar(as.character(paste0(dff))))
  dff$V1 <- gsub(">", "", dff$V1) %>% gsub(" .*", "", .)
  dff <- dff[dff$V1 %in% df$V1,]
  
  df <- merge(df, dff, by="V1", all=TRUE, sort=FALSE)
  
  illumina <- rbind(illumina, df)
}; rm(i, p, l, df, dff)

colnames(illumina) <- c("query", "target", "id", "alnlen", "mism", "opens", "qlo", "qhi", "tlo", "thi", "evalue", "bits", "cell", "replicate", "sequence", "length")


#----
#---- Creating a single file with all the information for plotting ---------------------------------

tmps <- select(sanger1st, c("cell", "replicate", "qseqid", "sseqid", "pident", "length")); tmps$method <- "Sanger"
tmpm <- select(minion1st, c("cell", "replicate", "qseqid", "sseqid", "pident", "length")); tmpm$method <- "MinION"
tmpi <- select(illumina, c("cell", "replicate", "query", "target", "id", "length")); tmpi$method <- "Illumina"
colnames(tmpi) <- c("cell", "replicate", "qseqid", "sseqid", "pident", "length", "method")

file <- rbind(tmps, tmpm, tmpi); rm(tmps, tmpm, tmpi)

file <- melt(file, id.vars=c("cell", "replicate", "qseqid", "sseqid", "method"))

file$method <- factor(file$method, levels=c("Sanger", "MinION", "Illumina"))

file$cell <- factor(file$cell, levels=c("Mge17-9", "Mge17-124", "Vil325", "Vil490", "Vil496", "Mge17-81", "Mge17-82", "Vil480", "Vil497", "NA"))

file$variable <- factor(file$variable, levels=c("length", "pident"))

file$order <- ifelse(file$cell=="Mge17-9"|file$cell=="Mge17-124"|file$cell=="Vil325"|file$cell=="Vil490"|file$cell=="Vil496", 
                     "Nassellaria", "Spumellaria")

#---- Saving file with the data --------------------------------------------------------------------

if(!dir.exists("data")){dir.create("data")}

write.table(file, "data/fig1.tsv", quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)


#----
#----
#---- Figure 1: Boxplots ---------------------------------------------------------------------------


ggplot(file, aes(x=method, y=value))+
  geom_boxplot()+
  #facet_wrap(~variable, scales="free")+
  facet_grid(variable~., scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))



(boxplots <- ggplot(file, aes(x=cell, y=value, fill=order))+
  geom_boxplot()+
  facet_grid(variable~method, scales="free")+
  scale_fill_manual(values=c("springgreen3", "steelblue3", "grey60"))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60, hjust=0.9, vjust=0.9))+
  theme(legend.position="none"))

if(!dir.exists("plots")){dir.create("plots")}

pdf(paste0("plots/Fig1_Boxplot_allVariables_perCell.pdf", sep=""), width=11.69, height=8.27, paper='special')
plot(boxplots)
dev.off()



#----


