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
setwd("~/Desktop/PhD/PhD_Data/_SingleCell/data/")
#----
#----
#---- Select MinION sequences ----------------------------------------------------------------------

file <-"minion/consensus/demultiplexed/minion_demultiplexed.fasta"

minion <-  seqinr::read.fasta(file, seqtype= "AA", as.string = T)
minion <- data.frame(name=gsub(">", "", paste(getAnnot(minion))), sequence=paste0(minion))
minion$file <- "MinION"
minion$read <- gsub(" .*", "", minion$name) %>% gsub(".*_", "", .)


# Removing MinION sequences that BLASTed other than Nassellaria and Spumellaria
tmp <- fread("minion/blast/raw_blast.tsv")
colnames(tmp) <- c("qseqid", "sseqid", "sacc", "stitle", "sscinames", "staxids", "sskingdoms", "sblastnames", "pident", "slen", "lengthAlign", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

all(minion$read %in% tmp$qseqid)
tmp <- tmp[tmp$qseqid %in% minion$read]

tmps <-data.frame()
for(i in unique(tmp$qseqid)){
  ss <- subset(tmp, qseqid==i)[1,]
  if(grepl("U|uncultured", ss$sscinames)){
    group <- unique(subset(tmp, qseqid==i)$sscinames)
    ss$group <- group[!grepl("U|uncultured", group)][1]
  }else{
    ss$group <- ss$sscinames
  }
  tmps <- rbind(tmps, ss)
}; rm(i, ss, group)

unique(tmps$group)
tmp <- unlist(c(tmps[grepl("Eucyrtidium|Styptosphaera|Rhizosphaera|Triastrum|Euchitonia", tmps$group),1]))

minion <- minion[grepl(paste(tmp, collapse="|"), minion$name) ,]; rm(tmp, tmps)

minion$nname <- gsub(" .*", "", minion$name) %>% paste("MinION", ., sep="-")


mspum <- minion[grepl("Mge", minion$name),]
mnass <- minion[grepl("Vil", minion$name),]

seqinr::write.fasta(sequences=as.list(mspum$sequence), names=mspum$nname, nbchar=80, file.out=paste0("files/entropy_rDNA/minion_demultiplexed_spum.fasta")) 
seqinr::write.fasta(sequences=as.list(mnass$sequence), names=mnass$nname, nbchar=80, file.out=paste0("files/entropy_rDNA/minion_demultiplexed_nass.fasta")) 

# And now the reference sanger alignment to align minion sequences
sanger <-  seqinr::read.fasta("sanger/raw/clones_concatenated_aligned.fasta", seqtype= "AA", as.string = T)
sanger <- data.frame(name=gsub(">", "", paste(getAnnot(sanger))), sequence=paste0(sanger))

sspum <- sanger[grepl("Mge", sanger$name),]
sspum$nname <- paste0("Sanger-", sspum$name)
snass <- sanger[grepl("Vil", sanger$name),]
snass$nname <- paste0("Sanger-", snass$name)

seqinr::write.fasta(sequences=as.list(sspum$sequence), names=sspum$nname, nbchar=80, file.out=paste0("files/entropy_rDNA/clones_concatenated_aligned_spum.fasta")) 
seqinr::write.fasta(sequences=as.list(snass$sequence), names=snass$nname, nbchar=80, file.out=paste0("files/entropy_rDNA/clones_concatenated_aligned_nass.fasta")) 




#----
#---- Align files ----------------------------------------------------------------------------------

setwd("files/entropy_rDNA/")

# But first remove all position that contains only gaps 

for(i in grep("clones", dir(), value=TRUE)){
  system(paste0("trimal -in ", i, " -out ", gsub("\\.fasta", "_trimmed.fasta", i), " -noallgaps"))
}

# This step might take a while ... 
# Also remember to check the number of threads !!!!

# For Nassellaria ____________________________________________________________________________________
new <- "minion_demultiplexed_nass.fasta"
reference <- "clones_concatenated_aligned_nass_trimmed.fasta"
out <- "entropy_nass_rDNA.fasta"

system(paste0("mafft",
              " --thread 4 --threadit 2 --inputorder --anysymbol --op 5 --ep 0.0 --leavegappyregion",
              " --add ./",     new, 
              " --localpair --maxiterate 100 ./", reference,
              " > ./",         out))


# For Spumellaira ____________________________________________________________________________________
new <- "minion_demultiplexed_spum.fasta"
reference <- "clones_concatenated_aligned_spum_trimmed.fasta"
out <- "entropy_spum_rDNA.fasta"

system(paste0("mafft",
              " --thread 4 --threadit 2 --inputorder --anysymbol --op 5 --ep 0.0 --leavegappyregion",
              " --add ./",     new, 
              " --localpair --maxiterate 100 ./", reference,
              " > ./",         out))

#________________________________________________________________________________________________#
#______________________________ Remember to heck manually !!!!!!!! ______________________________#
#________________ MinION sequences (or other positions) might screw the alignemnt _______________#
#________________________________________________________________________________________________#

setwd("../../")


#----
#---- Calculate Entropy ----------------------------------------------------------------------------

setwd("files/entropy_rDNA/")

# For clarity, move all files we are not interested in to another folder
dir.create("raw")
system(paste0("mv clones* minion* -t raw/"))

# And clear the environment
rm(list=ls())

# Load the function for the entropy ________________________________________________________________
# for further information check: https://github.com/MiguelMSandin/DNA-alignment-entropy ____________

positionInfo <- function(data, verbose=TRUE){
  info <- list()
  if(verbose){cat("Analyzing positions \n")}
  for (i in 1:(nchar(as.character(data$sequence[1])))){
    ss <- substr(data[,2],i,i)
    ss <- toupper(ss)
    df <- data.frame(base=unique(ss))
    tmp <- c()
    for ( j in 1:length(unique(ss))){
      tmp[j] <- length(grep(df[j,1], ss))
    }
    df$rep <- tmp
    rownames(df) <- as.character(df[,1])
    df[,1] <- NULL
    
    dff <- subset(df, rownames(df)!="-")
    
    info$shan[i] <- vegan::diversity(t(df))                                    # Shannon entropy of position
    info$shanc[i] <- ifelse(nrow(dff)==0, NA, vegan::diversity(t(dff)))        # Shannon entropy of position removing gaps ("-")
    info$rich[i] <- length(unique(ss))                                  # Position richness
    info$richc[i] <- ifelse(nrow(dff)==0, NA, length(rownames(dff)))    # Position richness removing gaps ("-")
    info$uniq[i] <- list(unique(ss))                                    # Unique bases in position
    info$repe[i] <- list(tmp)                                           # Repetitions of the unique bases in position
    
    if(verbose){
      if((i %% (round(nchar(as.character(data$sequence[1]))*0.1,0))) == 0){message("  ", round(i/nchar(as.character(data$sequence[1]))*100, 0), "%")}
    }
  }
  
  if(verbose){cat("Creating the data frame with the information \n")}
  
  df <- data.frame(posi=c(1:length(info$shan)),
                   shan=info$shan,
                   shanc=info$shanc,
                   rich=info$rich,
                   richc=info$richc,
                   uniq=unlist(lapply(info$uniq, function(x) paste(x, collapse="|"))),
                   repe=unlist(lapply(info$repe, function(x) paste(x, collapse="|"))))
  if(verbose){cat("Finished \n")}
  return(df)
  
}


# Choose the aligned fasta file
files <- c("entropy_nass_rDNA.fasta", "entropy_spum_rDNA.fasta")

entropy <- data.frame()
for(f in files){
  cat("        Analyzing", f, " \n")
  s <-  seqinr::read.fasta(f, seqtype= "AA", as.string = T)
  s <- data.frame(name=paste(getAnnot(s)), sequence=paste0(s))
  s$name <- gsub(">", "", s$name)
  
  (sample <- unique(gsub("_.*", "", s$name)))
  
  df <- data.frame()
  for(i in sample){
    cat("Reading sequences starting by '", i, "'. ", grep(paste0("^", i, "$"), sample), "/", length(sample), "\n", sep="")
    ss <- s[grepl(paste0("^", i, "_"), s$name),]
    dfs <- positionInfo(ss)
    dfs$sample <- i
    df <- rbind(df, dfs)
    cat("_____________________________", "\n")
  }; rm(i, ss, dfs)
  
  df$file <- gsub("_.*", "", f)
  entropy <- rbind(entropy, df)
  
  cat("Finished \n")
}; rm(f, s, sample, df)


write.table(entropy, "entropy.tsv", quote=FALSE, row.names=FALSE, sep="\t")


#----
#---- Setting variables and factors ----------------------------------------------------------------

setwd("../../")
entropy <- fread("files/entropy_rDNA/entropy.tsv")

entropy$sample <- factor(entropy$sample, levels=c("Sanger-Vil325", "Sanger-Vil496", "MinION-Vil325", "MinION-Vil496", 
                                                  "Sanger-Mge17-81", "Sanger-Mge17-82", "MinION-Mge17-81", "MinION-Mge17-82"))

entropy$method <- gsub("-.*", "", entropy$sample)
entropy$method <- factor(entropy$method, levels=c("Sanger", "MinION"))

entropy$order <- fifelse(grepl("Vil", entropy$sample), "Nassellaria", "Spumellaria")


#----
#---- Plotting -------------------------------------------------------------------------------------

dive <- ggplot(entropy, aes(x=posi, y=shanc, colour=sample))+
  geom_point(aes())+
  geom_smooth(level=.99)+
  facet_wrap(method~order, ncol=1)+
  scale_colour_manual(values=c("springgreen1", "springgreen4", "springgreen1", "springgreen4", 
                               "steelblue1", "steelblue4", "steelblue1", "steelblue4"))+
  theme_bw()+
  theme(legend.position="none")
dive

# Exporting the PDF of the plot
pdf("../plots/FigS5_alignment_entropy.pdf", width=11.69, height=8.27, paper='special')
plot(dive)
dev.off()
# But this graph is too heavy too edit in inkscape, so let's remove all 0 points (without affecting geom_smoot !!)

divep <- ggplot()+
  geom_point(subset(entropy, shanc>0), mapping=aes(x=posi, y=shanc, colour=sample))+
  geom_smooth(entropy, mapping=aes(x=posi, y=shanc, colour=sample), level=.95)+
  facet_wrap(method~order, ncol=1)+
  scale_colour_manual(values=c("springgreen1", "springgreen4", "springgreen1", "springgreen4", 
                               "steelblue1", "steelblue4", "steelblue1", "steelblue4"))+
  theme_bw()+
  theme(legend.position="none")
divep

# Exporting the PDF of the plot
pdf("../plots/FigS5_alignment_entropy_clean.pdf", width=11.69, height=8.27, paper='special')
plot(divep)
dev.off()


#----
#----





