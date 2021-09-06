#----
#---- Load libraries -------------------------------------------------------------------------------

library(dada2); packageVersion("dada2")

library(Biostrings)
library(ape)
library(magrittr)
library(data.table)
library(dplyr)
library(tidyr)

library(ggplot2)
library(treemapify)

#----
######### Set Working directory and file names -----------------------------------------------------

setwd("~/IntraGenomic-variability/data/illumina/dada2/")

project_name <- "default_pipeline/"

# Create directories that will be used to stored the files at the different stage of the processing
fastq_dir <-    "../demultiplexed/" 
filtered_dir <- "/filtered/"  # fastq filtered
qual_dir <- "/quality/"
dada2_dir <-    "/out/"           # dada2 results
database_dir <- "~/DB/PR2/pr2_version_4.14.0_SSU_dada2_V4.fasta" # databases
tax_levels <- c("Kingdom", "Supergroup","Division", "Class", "Order", "Family", "Genus", "Species")

dir.create(project_name)
dir.create(paste0(project_name, filtered_dir))
dir.create(paste0(project_name, dada2_dir))
dir.create(paste0(project_name, qual_dir))

# Forward and reverse fastq filenames
fnFs <- sort(list.files(fastq_dir, pattern="_R1_trimmed.fastq", full.names=TRUE))
fnRs <- sort(list.files(fastq_dir, pattern="_R2_trimmed.fastq", full.names=TRUE))

# Extract sample names
sample_names <- gsub("_R1_.*", "", basename(fnFs))

# Place filtered files in filtered/ subdirectory
filtFs <- paste0(project_name, filtered_dir, paste0(sample_names, "_R1_filt.fastq.gz"))
filtRs <- paste0(project_name, filtered_dir, paste0(sample_names, "_R2_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# Saving workspace
dir.create(paste0(project_name, "/workSpace"))
save(project_name, fastq_dir, filtered_dir, qual_dir, dada2_dir, database_dir, tax_levels, fnFs, fnRs, sample_names, filtFs, filtRs,
     file=paste0(project_name, "workSpace/01_directories.RData"))
#load(paste0(project_name, "workSpace/01_directories.RData"))


# Loading all previously saved workingSpaces if you want to run again some functions 
if(FALSE){
  for(i in dir(paste0(project_name, "/workSpace"))){
    cat("Loading file ", grep(i, dir(paste0(project_name, "/workSpace"))), "/", length(dir(paste0(project_name, "/workSpace"))), ": ", i, "\n", sep="")
    load(paste0(project_name, "/workSpace/", i))
  }; rm(i)
}


#----
#---- Check if files look OK after demultiplexing and trimming -------------------------------------


if(!length(fnFs)==length(fnRs)){stop(message("Different number of files in R1 and R2 file paths."))}


### Compute number of paired reads

df <- data.frame()  
# loop through all the R1 files and R2 files to check the number of reads
for(i in 1:length(fnFs)) { 
  geomF <- fastq.geometry(fnFs[i])    # use the dada2 function fastq.geometry
  geomR <- fastq.geometry(fnRs[i])    
  df_one_row <- data.frame(F_n_seq=geomF[1], F_file_name=basename(fnFs[i]), 
                           R_n_seq=geomR[1], R_file_name=basename(fnRs[i]))    # extract the information on number of sequences and file name 
  df <- rbind(df, df_one_row)    # add one line to data frame
}; rm(i, geomF, geomR, df_one_row)

write.table(df, file=paste0(project_name, "n_seq.tsv"), sep="\t", row.names = FALSE, na="", quote=FALSE)

if(!all(df$F_n_seq==df$R_n_seq)){stop(message("Different number of reads in R1 and R2 files."))}


# plot the histogram with number of sequences

p <- ggplot(df, aes(x=F_n_seq)) + 
  geom_histogram(alpha = 0.5, position="identity", binwidth = 10)+
  theme_bw()

pdf(paste0(project_name, "n_seq_histogram.pdf"), width=11.69, height=8.27/2, paper='special')
plot(p)
dev.off()

## Check if R1 and R2 fastq files contain reads in matched order

match <- data.frame()  
for(i in 1:length(fnFs)) { 
  tmpF <- read.fastq(fnFs[i])
  tmpR <- read.fastq(fnRs[i])
  lF <- gsub(" (1|2):", " d:", labels(tmpF))
  lR <- gsub(" (1|2):", " d:", labels(tmpR))
  match <- rbind(match, data.frame(file_F=fnFs[i], file_R=fnRs[i], match=all(lF==lR)))
}; rm(i, tmpF, tmpR, lF, lR)

if(!all(match$match)){
  stop("R1 and R2 fastq files contain reads in UN-MATCHING order.")
} else {
  cat("\n ## R1 and R2 files look good. Let's proceed with DADA2 pipeline. :) \n")
}

rm(df, p, match)

#----
#---- Plot quality ----------------------------------------------------------------------------------


pdf(paste0(project_name, qual_dir, "indiv_F_Qplots.pdf"), width=11.69, height=8.27, paper='special')
for(i in file.path(fnFs)){
  cat(i)
  print(plotQualityProfile(i))
}
dev.off()

pdf(paste0(project_name, qual_dir, "indiv_R_Qplots.pdf"), width=11.69, height=8.27, paper='special')
for(i in file.path(fnRs)){
  print(plotQualityProfile(i))
}; rm(i)
dev.off()

pdf(paste0(project_name, qual_dir, "global_F_Qplots.pdf"), width=11.69, height=8.27, paper='special')
plotQualityProfile(file.path(fnFs),aggregate=T)
dev.off()

pdf(paste0(project_name, qual_dir, "global_R_Qplots.pdf"), width=11.69, height=8.27, paper='special')
plotQualityProfile(file.path(fnRs),aggregate=T)
dev.off()

#----
############################################## DADA-2 ##############################################
#### Filtering reads                   ######################################################

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=0,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=2)
head(out)

# Saving workspace
save(out, file=paste0(project_name, "workSpace/02_filterAndTrim.RData"))
#load(paste0(project_name, "workSpace/02_filterAndTrim.RData"))


#### Learn the error rates             ######################################################

# Learning error rates for R1 files
errF <- learnErrors(filtFs, multithread=2)

pdf(paste0(project_name, "Errors_F.pdf"), width=11.69, height=8.27, paper='special')
plotErrors(errF, nominalQ=TRUE)
dev.off()

save(errF, file=paste0(project_name, "workSpace/03_ErrorF.RData"))
#load(paste0(project_name, "workSpace/03_ErrorF.RData"))


# Learning error rates for R2 files
errR <- learnErrors(filtRs, multithread=2)

pdf(paste0(project_name, "Errors_R.pdf"), width=11.69, height=8.27, paper='special')
plotErrors(errR, nominalQ=TRUE)
dev.off()

save(errR, file=paste0(project_name, "workSpace/03_ErrorR.RData"))
#load(paste0(project_name, "workSpace/03_ErrorR.RData"))


#### Running DADA2 algorithm           ######################################################

# But first dereplicate fastq files    ______________________________________________________

# Dereplicating R1 files
derepFs <- derepFastq(filtFs)

# Dereplicating R2 files
derepRs <- derepFastq(filtRs)

# Name the derep-class objects by the sample names
names(derepFs) <- sample_names
names(derepRs) <- sample_names

save(derepFs, derepRs, file=paste0(project_name, "workSpace/04_derep.RData"))
#load(paste0(project_name, "workSpace/04_derep.RData"))


# Running DADA2 algorithm for R1 files
dadaFs <- dada(derepFs, err=errF, multithread=2)

save(dadaFs, file=paste0(project_name, "workSpace/04_dadaF.RData"))
#load(paste0(project_name, "workSpace/04_dadaF.RData"))

dadaFs[[1]]


# Running DADA2 algorithm for R2 files
dadaRs <- dada(derepRs, err=errR, multithread=2)

save(dadaRs, file=paste0(project_name, "workSpace/04_dadaR.RData"))
#load(paste0(project_name, "workSpace/04_dadaR.RData"))

dadaRs[[1]]


#### Merge paired reads                ######################################################

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

save(mergers, file=paste0(project_name, "workSpace/05_mergers.RData"))
#load(paste0(project_name, "workSpace/05_mergers.RData"))

# Inspect the merger data.frame from the first sample
head(mergers[[1]])


#### Construct sequence table          ######################################################

seqtab <- makeSequenceTable(mergers)

save(seqtab, file=paste0(project_name, "workSpace/06_seqtab.RData"))
#load(paste0(project_name, "workSpace/06_seqtab.RData"))

dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#### Remove chimeras                   ######################################################

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=2, verbose=TRUE)

save(seqtab.nochim, file=paste0(project_name, "workSpace/07_seqtabNoChim.RData"))
#load(paste0(project_name, "workSpace/07_seqtabNoChim.RData"))

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
cat("    Chimeras accounted for the ", 100-round(sum(seqtab.nochim)/sum(seqtab)*100, 2), "% of the merged sequence reads. \n", sep="")
cat("    Total number of sequences: ", sum(seqtab.nochim), ". \n", sep="")


#### Track reads through the pipeline  ######################################################

getN <- function(x){sum(getUniques(x))}
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
track <- as.data.frame(track)
track$outPerCent <- track$nonchim/track$input*100
rownames(track) <- sample_names
head(track)

write.table(track, file=paste0(project_name, "track_sequences.tsv"), sep="\t", na="", quote=FALSE)


#### Assign taxonomy                   ######################################################

taxa <- assignTaxonomy(seqtab.nochim, database_dir, taxLevels=tax_levels, outputBootstraps=TRUE, multithread=2, verbose=TRUE)

# if problems with the memory:
# write.table(seqtab.nochim, paste0(project_name, dada2_dir, "seqtab.tsv"), sep="\t", quote=FALSE)
# And send to cluster. Check files "B2.1_dada_assignTaxo.R" and "B2.1_dada_assignTaxo.sh"


save(taxa, file=paste0(project_name, "workSpace/08_taxa.RData"))
#load(paste0(project_name, "workSpace/08_taxa.RData"))

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


#----
#### Exporting output                  ######################################################    

write.table(taxa$tax, paste0(project_name, dada2_dir, "taxa.tsv"), sep="\t", quote=FALSE)  
write.table(taxa$boot, paste0(project_name, dada2_dir, "taxa_boot.tsv"), sep="\t", quote=FALSE)

write.table(seqtab.nochim, paste0(project_name, dada2_dir, "seqtab_raw.tsv"), sep="\t", quote=FALSE)
write.table(t(seqtab.nochim), paste0(project_name, dada2_dir, "seqtab.tsv"), sep="\t", quote=FALSE)


fasta <- data.frame(name=paste0("asv", seq(1:ncol(seqtab.nochim))), sequence=paste0(colnames(seqtab.nochim)))
seqinr::write.fasta(sequences=as.list(fasta$sequence), names=fasta$name, nbchar=80, file.out=paste0(project_name, dada2_dir, "ASVs.fasta"))

##### Exporting ASVs + taxonomy + abundance table ############################################

file <- data.frame(ASV=paste0("asv", seq(1:ncol(seqtab.nochim))), 
                   sequence=paste0(colnames(seqtab.nochim)),
                   as.data.frame(taxa$tax),
                   min_BS=apply(taxa$boot, 1, min),
                   tpres=apply(t(seqtab.nochim), 1, function(x) sum(x > 0)),
                   tabun=apply(t(seqtab.nochim), 1, sum),
                   as.data.frame(t(seqtab.nochim)))
rownames(file) <- seq(1:nrow(file))

colnames(file) <- gsub("\\.", "-", colnames(file))

write.table(file, paste0(project_name, dada2_dir, "ASVs_final_tax_abun.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#----



