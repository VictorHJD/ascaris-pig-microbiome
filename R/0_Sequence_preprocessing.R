##Project: Ascaris - Pig Microbiome
##Aim: Sequencing pre-processing and denoising.
##Author: Víctor Hugo Jarquín-Díaz

##Load libraries 

library(ggplot2)
library(dada2)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ShortRead)

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doQualEval <- TRUE

doFilter <- FALSE

doTax <- FALSE

doPS <- FALSE

###################dada2 pipeline#######################
# Setting seed for reproducibility
set.seed(022019)

#######################################################
path <- c(
  ## Ascaris-Pig (16S) Full sequencing Run (Good run)
  "2021_16_Ascaris_Main", 
  ## Ascaris-Pig (16S) Preliminary test sequencing Run
  "2021_16_Ascaris_Quanti") 

fullpath <- paste0("/fast/AG_Forslund/Victor/data/Ascaris/", path)

names(fullpath) <- path

fastqList <- lapply(fullpath, function (path) { 
  fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) 
  fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE)
  fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE)
  list(fastqF=fastqF, fastqR=fastqR)
})

if(doQualEval){
  readlenght <- lapply(fastqList, function (x) {
    con <- file(x[["fastqF"]][1],"r")
    ## first line
    secondLine <- readLines(con, n=2)[[2]]
    ### simple check whether it's V2 or V3 data
    nchar(secondLine)
  })
  
  allFastqF <- lapply(fastqList, function (x) {
    readFastq(x[["fastqF"]])
  })
  
  allFastqR <- lapply(fastqList, function (x) {
    readFastq(x[["fastqR"]])
  })
  
  
  sampleQual <- function (x) {
    ## sample quality scores of 100,000 sequences 
    qmat <- as(quality(x)[sample(100000)], "matrix")
    cols <- seq(1, ncol(qmat), by=10)
    sapply(cols, function (i) {
      mean(qmat[, i], na.rm=TRUE)
    })
  }
  
  qualityF <- lapply(allFastqF, sampleQual)
  qualityR <- lapply(allFastqR, sampleQual)
  
  shouldL <- max(unlist(lapply(qualityF, length)))
  
  qualityFilledF <- lapply(qualityF, function (x) {
    c(x, rep(NA, times=shouldL - length(x)))
  })
  
  qualityFilledR <- lapply(qualityR, function (x) {
    c(x, rep(NA, times=(shouldL - length(x))))
  })
  
  qualityDFF <- Reduce("cbind",  qualityFilledF)
  qualityDFR <- Reduce("cbind",  qualityFilledR)
  
  colnames(qualityDFF) <- path
  colnames(qualityDFR) <- path
  
  qualityDFFL <- reshape2::melt(qualityDFF)
  qualityDFFL$direction <- "forward"
  
  qualityDFRL <- reshape2::melt(qualityDFR)
  qualityDFRL$direction <- "reverse"
  
  qualityDFL <- rbind(qualityDFFL, qualityDFRL)
  
  qualityDFL$position <- qualityDFL$Var1*10 -10
  
  # New facet label names for direction variable
  dir.labs <- c("Forward", "Reverse")
  names(dir.labs) <- c("forward", "reverse")
  
  ggplot(qualityDFL, aes(position, value, color=Var2)) +
    geom_line() +
    facet_wrap(~direction, labeller = labeller(direction= dir.labs))+
    theme_bw()+
    xlab("Position")+
    ylab("Quality score")+
    geom_vline(xintercept = 240, colour= 'black', linetype="dashed") +                                                                                                     
    geom_hline(yintercept = 30, colour= 'black', linetype="dashed") + 
    scale_color_discrete(name = "Run name", labels = c("Main", "Test"))-> A

}

rm(qualityDFF, qualityDFFL, qualityDFL, qualityDFR, qualityDFRL, qualityF, qualityR, qualityFilledR, qualityFilledF, dir.labs)

###############################################
## concluding from this that we can truncate: 
## at 240 for Rev # at 240 for Fwd for Main run
## at 200 fro Fwd # at 170 for Rev for Quanti run 

samplesList <- lapply (fastqList, function (x){
  samples<- gsub("-", "_", basename(x[["fastqF"]]))
  samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", samples)
  samples<- gsub("S\\d+_", "\\1", samples)
  paste(basename(dirname(x[["fastqF"]])), samples, sep="_")
})

fastqFall <- unlist(lapply(fastqList, "[[", "fastqF"))
fastqRall <- unlist(lapply(fastqList, "[[", "fastqR"))

samplesAll <- unlist(samplesList)

samplesAll<- gsub("2021_16_Ascaris_Main_", "P2_", samplesAll)
samplesAll<- gsub("2021_16_Ascaris_Quanti_", "P1_", samplesAll)

samples<- gsub("P2_", "", samplesAll)
samples<- gsub("P1_", "", samples)
samples<- unique(samples)

# Reading in sequences, quality filtering and trimming.
## Run 1: Quantification data (P1)

path1<- fullpath[[2]]

# File parsing
fastqF1 <- fastqList$`2021_16_Ascaris_Quanti`$fastqF
fastqR1 <- fastqList$`2021_16_Ascaris_Quanti`$fastqR
if(length(fastqF1) != length(fastqR1)) stop("Forward and reverse files do not match.")

##Individual quality check
plotQualityProfile(fastqF1[1:12])
plotQualityProfile(fastqR1[1:12])

#Creation of a folder for filtrated reads 
filt_path1<- "/fast/AG_Forslund/Victor/data/Ascaris/filtered_Ascaris_Quanti"
if(!file_test("-d", filt_path1)) dir.create(filt_path1)

filtFs1 <- file.path(filt_path1, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs1) <- samples
filtRs1 <- file.path(filt_path1, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs1) <- samples

# Filtering: The parameters are run specific based on the quality assessment

if(doFilter){
  filter_track_1 <- filterAndTrim(fwd=file.path(fastqF1), filt=file.path(filtFs1),
                rev=file.path(fastqR1), filt.rev=file.path(filtRs1),
                truncLen=c(200,170), 
                maxEE=c(2,2), truncQ=2, trimLeft = c(17, 21), ##Remove primers
                rm.phix=TRUE,
                compress=TRUE, verbose=TRUE, multithread=TRUE, matchIDs=TRUE)  ## forward and reverse not matching otherwise 
}

##Check the proportion of reads that passed the filtering 
sum(filter_track_1[,"reads.out"])/sum(filter_track_1[,"reads.in"])

### Over 82% passed for Run 1
## Check which one passed
filtFiles1 <- list.files(filt_path1, pattern=".fastq.gz$", full.names=TRUE) 
filtFs1 <- grep("_F_filt.fastq.gz", filtFiles1, value = TRUE)
filtRs1 <- grep("_R_filt.fastq.gz", filtFiles1, value = TRUE)
#Quality after trimming
if(length(filtFs1) != length(filtRs1)) stop("Forward and reverse files do not match.")

plotQualityProfile(filtFs1[1:12])
plotQualityProfile(filtRs1[1:12])

## Run 2: Main run data (P1)

path2<- fullpath[[1]]

# File parsing
fastqF2 <- fastqList$`2021_16_Ascaris_Main`$fastqF
fastqR2 <- fastqList$`2021_16_Ascaris_Main`$fastqR
if(length(fastqF2) != length(fastqR2)) stop("Forward and reverse files do not match.")

##Individual quality check
plotQualityProfile(fastqF2[1:12])
plotQualityProfile(fastqR2[1:12])

#Creation of a folder for filtrated reads 
filt_path2<- "/fast/AG_Forslund/Victor/data/Ascaris/filtered_Ascaris_Main"
if(!file_test("-d", filt_path2)) dir.create(filt_path2)

filtFs2 <- file.path(filt_path2, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs2) <- samples
filtRs2 <- file.path(filt_path2, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs2) <- samples

if(doFilter){
  filter_track_2 <- filterAndTrim(fwd=file.path(fastqF2), filt=file.path(filtFs2),
                                  rev=file.path(fastqR2), filt.rev=file.path(filtRs2),
                                  truncLen=c(245,245), 
                                  maxEE=c(2,2), truncQ=2, trimLeft = c(17, 21), ##Remove primers
                                  rm.phix=TRUE,
                                  compress=TRUE, verbose=TRUE, multithread=TRUE, matchIDs=TRUE)  ## forward and reverse not matching otherwise 
}

##Check the proportion of reads that passed the filtering 
sum(filter_track_2[,"reads.out"])/sum(filter_track_2[,"reads.in"])

### Over 78% passed for Run 2
#Quality after trimming
## Check which one passed
filtFiles2 <- list.files(filt_path2, pattern=".fastq.gz$", full.names=TRUE) 
filtFs2 <- grep("_F_filt.fastq.gz", filtFiles2, value = TRUE)
filtRs2 <- grep("_R_filt.fastq.gz", filtFiles2, value = TRUE)
#Quality after trimming
if(length(filtFs2) != length(filtRs2)) stop("Forward and reverse files do not match.")

plotQualityProfile(filtFs2[1:12])
plotQualityProfile(filtRs2[1:12])

# Infer Sequence Variants
#This should be run on a run-by-run basis as not all runs will have the same error profiles
## Run 1
set.seed(100)
# Learn forward error rates
errF_1 <- learnErrors(filtFs1, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR_1 <- learnErrors(filtFs1, nbases=1e8, multithread=TRUE)

## Run 2
# Learn forward error rates
errF_2 <- learnErrors(filtFs2, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR_2 <- learnErrors(filtFs2, nbases=1e8, multithread=TRUE)

#Let's look at the error profiles for each of the dada2 runs
plotErrors(errF_1, nominalQ=TRUE)
plotErrors(errF_2, nominalQ=TRUE)
plotErrors(errR_1, nominalQ=TRUE)
plotErrors(errR_2, nominalQ=TRUE)

# Sample Inference
#apply the core inference algorithm to the filtered and trimmed sequence data
## Run1
dadaF1 <- dada(filtFs1, err=errF_1, multithread=TRUE)
dadaR1 <- dada(filtRs1, err=errR_1, multithread=TRUE)

## Run2
dadaF2 <- dada(filtFs2, err=errF_2, multithread=TRUE)
dadaR2 <- dada(filtRs2, err=errR_2, multithread=TRUE)

# Merge sequences and make tables
## Run 1
mergers_1 <- mergePairs(dadaF1, filtFs1, dadaR1, filtRs1, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers_1[[1]])
seqtab_1 <- makeSequenceTable(mergers_1)

##Check size of fragments
table(nchar(getSequences(seqtab_1))) ##--> Amplicon size ranges between 183 to 320, ~426bp is expected as amplicon
##Weird but let's have a look at the final outcome and then decide whether is useful to include this data.

###Remove of the chimeras,
##1) Per-sample: The samples in a sequence table are independently checked for bimeras,
#and sequence variants are removed (zeroed-out) from samples independently
seqtab_nochim_1 <- removeBimeraDenovo(seqtab_1, method="per-sample", multithread=TRUE)

##2) Pooled: The samples in the sequence table are all pooled together for bimera
#identification
seqtab_nochim_1 <- removeBimeraDenovo(seqtab_nochim_1, method="pooled", multithread=TRUE)

#Look at fraction of chimeras. 
dim(seqtab_1)
## 762 ASVs before chimera removal
dim(seqtab_nochim_1)
## 718 ASVs after chimera removal in two steps
sum(seqtab_nochim_1)/sum(seqtab_1)
#Here, chimeras made up about 5.8% of the ASVs, but that was only about 1% of total sequence reads

## Run 2
mergers_2 <- mergePairs(dadaF2, filtFs2, dadaR2, filtRs2, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers_2[[1]])
seqtab_2 <- makeSequenceTable(mergers_2)
##Check size of fragments
table(nchar(getSequences(seqtab_2))) ##--> Amplicon size ranges between 228 to 440, ~426bp is expected as amplicon

###Remove of the chimeras,
##1) Per-sample
#seqtab_nochim_2 <- removeBimeraDenovo(seqtab_2, method="per-sample", multithread=TRUE)
# To computational demanding :(
#Error in mcfork() : 
#unable to fork, possible reason: Cannot allocate memory

##1) Pooled
seqtab_nochim_2 <- removeBimeraDenovo(seqtab_2, method="pooled", multithread=TRUE)
##2) Consensus (default dada2)
seqtab_nochim_2 <- removeBimeraDenovo(seqtab_nochim_2, method="consensus", multithread=TRUE)

#Look at fraction of chimeras. 
dim(seqtab_2)
## 75369 ASVs before chimera removal
dim(seqtab_nochim_2)
## 5521 ASVs after chimera removal in two steps
sum(seqtab_nochim_2)/sum(seqtab_2)
#Here, a lot of chimeras!!! They made up about 92.67% of the ASVs, and that was only about 33% of total sequence reads!!! 

#Track Reads through pipeline (come back to this. Have to make for each run)
getN <- function(x) sum(getUniques(x))
#Run1
colnames(filter_track_1)<- c("input", "filtered")
sampleNames_1<- rownames(filter_track_1)
sampleNames_1 <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", sampleNames_1)
sampleNames_1<- gsub("-", "_", sampleNames_1)
sampleNames_1<- gsub("S\\d+_", "\\1", sampleNames_1)
rownames(filter_track_1) <- sampleNames_1

as.data.frame(filter_track_1)%>%
  tibble::rownames_to_column(var = "Barcode_name")-> tmp

track_1 <- cbind(sapply(dadaF1, getN), sapply(dadaR1, getN), sapply(mergers_1, getN), rowSums(seqtab_nochim_1))
colnames(track_1) <- c("denoisedF", "denoisedR", "merged", "nochim")
sampleNames_1<- rownames(track_1)
sampleNames_1 <- gsub("_F_filt.fastq\\.gz", "\\1", sampleNames_1)
rownames(track_1) <- sampleNames_1

as.data.frame(track_1)%>%
  tibble::rownames_to_column(var = "Barcode_name")%>%
  right_join(tmp, by = "Barcode_name", all = TRUE)%>%
  relocate(filtered)%>%
  relocate(input)%>%
  relocate(Barcode_name)->track_1

#Run2
colnames(filter_track_2)<- c("input", "filtered")
sampleNames_2<- rownames(filter_track_2)
sampleNames_2 <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", sampleNames_2)
sampleNames_2<- gsub("-", "_", sampleNames_2)
sampleNames_2<- gsub("S\\d+_", "\\1", sampleNames_2)
rownames(filter_track_2) <- sampleNames_2

as.data.frame(filter_track_2)%>%
  tibble::rownames_to_column(var = "Barcode_name")-> tmp

track_2 <- cbind(sapply(dadaF2, getN), sapply(dadaR2, getN), sapply(mergers_2, getN), rowSums(seqtab_nochim_2))
colnames(track_2) <- c("denoisedF", "denoisedR", "merged", "nochim")
sampleNames_2<- rownames(track_2)
sampleNames_2 <- gsub("_F_filt.fastq\\.gz", "\\1", sampleNames_2)
rownames(track_2) <- sampleNames_2

as.data.frame(track_2)%>%
  tibble::rownames_to_column(var = "Barcode_name")%>%
  right_join(tmp, by = "Barcode_name", all = TRUE)%>%
  relocate(filtered)%>%
  relocate(input)%>%
  relocate(Barcode_name)%>% 
  mutate(perc_original_sequences = nochim/input*100)->track_2

##Save track files

write_csv(track_1, "Tables/Ascaris_Quanti_Tracking.csv")
write_csv(track_2, "Tables/Ascaris_Main_Tracking.csv")

#Save sequence tables
saveRDS(seqtab_nochim_1, "Data/seqnochim_Ascaris_Quanti.rds")
saveRDS(seqtab_nochim_2, "Data/seqnochim_Ascaris_Quanti.rds")

# Assign Taxonomy 
## SILVA train set
if(doFilter){
  
# Assign taxonomy SILVA Train Set
tax_silva_1 <- assignTaxonomy(seqtab_nochim_1, "~/vjarqui/SILVA_db/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
colnames(tax_silva_1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

tax_silva_2 <- assignTaxonomy(seqtab_nochim_2, "~/vjarqui/SILVA_db/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
colnames(tax_silva_2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

## Add species assignment to taxonomy table
tax_species_silva_1 <- addSpecies(tax_silva_1, "~/vjarqui/SILVA_db/silva_species_assignment_v138.fa.gz", verbose=TRUE, allowMultiple = FALSE)
colnames(tax_species_silva_1)  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unname(head(tax_species_silva_1))

tax_species_silva_2 <- addSpecies(tax_silva_2, "~/vjarqui/SILVA_db/silva_species_assignment_v138.fa.gz", verbose=TRUE, allowMultiple = FALSE)
colnames(tax_species_silva_2)  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unname(head(tax_species_silva_2))

# Write to disk
saveRDS(tax_species_silva_1, "Data/tax_species_final_Ascaris_Quanti.rds")
saveRDS(tax_species_silva_2, "Data/tax_species_final_Ascaris_Main.rds")

}