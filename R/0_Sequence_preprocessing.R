##Project: Ascaris - Pig Microbiome
##Aim: Sequencing pre-processing and denoising.
##Author: Víctor Hugo Jarquín-Díaz
##Root repo setwd("~/Ascaris/ascaris/")
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

rm(qualityDFF, qualityDFFL, qualityDFL, qualityDFR, qualityDFRL, 
   qualityF, qualityR, qualityFilledR, qualityFilledF, readlenght, dir.labs)

###############################################
## concluding from this that we can truncate: 
## at 240 for Rev # at 240 for Fwd for Main run
## at 240 fro Fwd # at 200 for Rev for Quanti run 

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
if(length(fastqF1) != length(fastqR1)) stop("Forward and reverse files do not match")

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
if(length(filtFs1) != length(filtRs1)) stop("Forward and reverse files do not match")

# Filtering: The parameters are run specific based on the quality assessment

if(doFilter){
  filter_track_1 <- filterAndTrim(fwd=file.path(fastqF1), filt=file.path(filtFs1),
                rev=file.path(fastqR1), filt.rev=file.path(filtRs1),
                truncLen=c(240,240), maxN=0,
                maxEE=c(2,2), truncQ=2, trimLeft = c(17, 21), ##Remove primers
                rm.phix=TRUE,
                compress=TRUE, verbose=TRUE, multithread=TRUE, matchIDs=TRUE)  ## forward and reverse not matching otherwise 
}

##Check the proportion of reads that passed the filtering 
sum(filter_track_1[,"reads.out"])/sum(filter_track_1[,"reads.in"])

### Over 36% passed for Run 1
## Check which one passed
filtFiles1 <- list.files(filt_path1, pattern=".fastq.gz$", full.names=TRUE) 
filtFs1 <- grep("_F_filt.fastq.gz", filtFiles1, value = TRUE)
filtRs1 <- grep("_R_filt.fastq.gz", filtFiles1, value = TRUE)
#Quality after trimming
if(length(filtFs1) != length(filtRs1)) stop("Forward and reverse files do not match")

plotQualityProfile(filtFs1[1:12])
plotQualityProfile(filtRs1[1:12])

## Run 2: Main run data (P1)

path2<- fullpath[[1]]

# File parsing
fastqF2 <- fastqList$`2021_16_Ascaris_Main`$fastqF
fastqR2 <- fastqList$`2021_16_Ascaris_Main`$fastqR
if(length(fastqF2) != length(fastqR2)) stop("Forward and reverse files do not match")

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
                                  truncLen=c(240,240),  maxN=0,
                                  maxEE=c(2,2), truncQ=2, trimLeft = c(17, 21), ##Remove primers
                                  rm.phix=TRUE,
                                  compress=TRUE, verbose=TRUE, multithread=TRUE, matchIDs=TRUE)  ## forward and reverse not matching otherwise 
}

##Check the proportion of reads that passed the filtering 
sum(filter_track_2[,"reads.out"])/sum(filter_track_2[,"reads.in"])

### Over 79% passed for Run 2
#Quality after trimming
## Check which one passed
filtFiles2 <- list.files(filt_path2, pattern=".fastq.gz$", full.names=TRUE) 
filtFs2 <- grep("_F_filt.fastq.gz", filtFiles2, value = TRUE)
filtRs2 <- grep("_R_filt.fastq.gz", filtFiles2, value = TRUE)
#Quality after trimming
if(length(filtFs2) != length(filtRs2)) stop("Forward and reverse files do not match")

plotQualityProfile(filtFs2[1:12])
plotQualityProfile(filtRs2[1:12])

# Infer Sequence Variants
#This should be run on a run-by-run basis as not all runs will have the same error profiles
## Run 1
set.seed(100)
# Learn forward error rates
errF_1 <- learnErrors(filtFs1, nbases=1e8, multithread=TRUE) ##Estimation based on 244 samples
# Learn reverse error rates
errR_1 <- learnErrors(filtFs1, nbases=1e8, multithread=TRUE)

## Run 2
# Learn forward error rates
errF_2 <- learnErrors(filtFs2, nbases=1e9, multithread=TRUE) ##Estimation based on 195 samples
# Learn reverse error rates
errR_2 <- learnErrors(filtFs2, nbases=1e9, multithread=TRUE)

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
table(nchar(getSequences(seqtab_1))) ##--> Amplicon size ranges between 223 to 429, ~426bp is expected as amplicon

###Remove of the chimeras,
##1) Per-sample: The samples in a sequence table are independently checked for bimeras,
#and sequence variants are removed (zeroed-out) from samples independently --> To computational demanding for large dataset (not for now)
#seqtab_nochim_1 <- removeBimeraDenovo(seqtab_1, method="per-sample", multithread=TRUE)

##1) Pooled: The samples in the sequence table are all pooled together for bimera
#identification
seqtab_nochim_1 <- removeBimeraDenovo(seqtab_1, method="pooled", multithread=TRUE)

##2) Consensus (default dada2)
seqtab_nochim_1 <- removeBimeraDenovo(seqtab_nochim_1, method="consensus", multithread=TRUE)

#Look at fraction of chimeras. 
dim(seqtab_1)
## 4164 ASVs before chimera removal
dim(seqtab_nochim_1)
## 1141 ASVs after chimera removal in two steps
sum(seqtab_nochim_1)/sum(seqtab_1)
#Here, chimeras made up about 72.6% of the ASVs, but that was only about ~13% of total sequence reads

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
## 108383 ASVs before chimera removal
dim(seqtab_nochim_2)
## 8181 ASVs after chimera removal in two steps
sum(seqtab_nochim_2)/sum(seqtab_2)
#Here, a lot of chimeras!!! They made up about 92.45% of the ASVs, and that was only about 49% of total sequence reads!!! 

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

##Create an ASV matrix with ASV as rows and samples as columns
##Run 1
asvmat1 <- t(seqtab_nochim_1) #Removing sequence rownames for display only
asv.names1 <- as.data.frame(rownames(asvmat1))
rownames(asv.names1) <- paste0("ASV", 1:nrow(asv.names1))
rownames(asvmat1)<-NULL
rownames(asvmat1) <- paste0("ASV", 1:nrow(asvmat1))
colnames(asvmat1) <- gsub("_F_filt.fastq\\.gz", "\\1", colnames(asvmat1))
head(asvmat1)
write.csv(asvmat1, "/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Quanti_ASV_matrix.csv")

##Run 2 
asvmat2 <- t(seqtab_nochim_2) #Removing sequence row names for display only
asv.names2 <- as.data.frame(rownames(asvmat2))
rownames(asv.names2) <- paste0("ASV", 1:nrow(asv.names2))
rownames(asvmat2)<-NULL
rownames(asvmat2) <- paste0("ASV", 1:nrow(asvmat2))
colnames(asvmat2) <- gsub("_F_filt.fastq\\.gz", "\\1", colnames(asvmat2))
head(asvmat2)
write.csv(asvmat2, "/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Main_ASV_matrix.csv")

##Get count of ASVs detected by sample
##Run 1
asv.sample<- as.data.frame(asvmat1)
test<- data.frame()
for (i in 1:ncol(asv.sample)) {
  asv<- data.frame()
  asv[1,1]<- sum(asv.sample[,i]!=0)
  rownames(asv)<- paste0("Sample", i)
  test <- rbind(test, asv) ### Join all the "individual" data frames into the final data frame 
}
asv.sample<- as.matrix(test)
colnames(asv.sample)<- "ASVs_dada2"
rownames(asv.sample)<- sampleNames_1

as.data.frame(asv.sample)%>%
  tibble::rownames_to_column(var = "Barcode_name")%>%
  right_join(track_1, by = "Barcode_name", all = TRUE)%>%
  mutate(perc_original_sequences = nochim/input*100)%>%
  mutate(perc_ASV_sequences = ASVs_dada2/nochim*100)->track_1

##Run 2
asv.sample<- as.data.frame(asvmat2)
test<- data.frame()
for (i in 1:ncol(asv.sample)) {
  asv<- data.frame()
  asv[1,1]<- sum(asv.sample[,i]!=0)
  rownames(asv)<- paste0("Sample", i)
  test <- rbind(test, asv) ### Join all the "individual" data frames into the final data frame 
}
asv.sample<- as.matrix(test)
colnames(asv.sample)<- "ASVs_dada2"
rownames(asv.sample)<- sampleNames_2

as.data.frame(asv.sample)%>%
  tibble::rownames_to_column(var = "Barcode_name")%>%
  right_join(track_2, by = "Barcode_name", all = TRUE)%>%
  mutate(perc_original_sequences = nochim/input*100)%>%
  mutate(perc_ASV_sequences = ASVs_dada2/nochim*100)->track_2

rm(test,asv.sample, asv, asvmat1, asvmat2,i)

##Merge both track objects 
#track_all <- rbind(track_1, track_2)

##Save track files
write.csv(track_1, "Tables/Ascaris_Quanti_Tracking.csv")
write.csv(track_2, "Tables/Ascaris_Main_Tracking.csv")

#Save sequence tables
saveRDS(seqtab_nochim_1, "/fast/AG_Forslund/Victor/data/Ascaris/tmp/seqnochim_Ascaris_Quanti.rds")
saveRDS(seqtab_nochim_2, "/fast/AG_Forslund/Victor/data/Ascaris/tmp/seqnochim_Ascaris_Main.rds")

# Assign Taxonomy 
## SILVA train set
if(doTax){
  
# Assign taxonomy SILVA Train Set
tax_silva_1 <- assignTaxonomy(seqtab_nochim_1, "/home/vjarqui/SILVA_db/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
colnames(tax_silva_1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

tax_silva_2 <- assignTaxonomy(seqtab_nochim_2, "/home/vjarqui/SILVA_db/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
colnames(tax_silva_2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

## Add species assignment to taxonomy table
tax_species_silva_1 <- addSpecies(tax_silva_1, "/home/vjarqui/SILVA_db/silva_species_assignment_v138.fa.gz", verbose=TRUE, allowMultiple = FALSE)
colnames(tax_species_silva_1)  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unname(head(tax_species_silva_1))

tax_species_silva_2 <- addSpecies(tax_silva_2, "/home/vjarqui/SILVA_db/silva_species_assignment_v138.fa.gz", verbose=TRUE, allowMultiple = FALSE)
colnames(tax_species_silva_2)  <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unname(head(tax_species_silva_2))

##Extract ASV sequences
library(DECIPHER); packageVersion("DECIPHER")

##Run 1
dna1 <- DNAStringSet(getSequences(seqtab_nochim_1)) # Create a DNAStringSet from the ASVs
names(dna1)<- paste0("ASV", 1:length(dna1)) ##Give short names to each sequence
writeXStringSet(dna1, "/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Quanti_ASV.fasta") ##Export fasta seq

##Run 2
dna2 <- DNAStringSet(getSequences(seqtab_nochim_2)) # Create a DNAStringSet from the ASVs
names(dna2)<- paste0("ASV", 1:length(dna2)) ##Give short names to each sequence
writeXStringSet(dna2, "/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Main_ASV.fasta") ##Export fasta seq

##Alternative taxa assignment
##Using decipher
load("/home/vjarqui/SILVA_db/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
##Run 1
ids <- IdTaxa(dna1, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid1 <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid1) <- ranks; rownames(taxid1) <- getSequences(seqtab_nochim_1)

##Run 2
ids <- IdTaxa(dna2, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid2 <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid2) <- ranks; rownames(taxid2) <- getSequences(seqtab_nochim_2)

# Write to disk
saveRDS(tax_silva_1, "Data/tax_final_Ascaris_Quanti.rds")
saveRDS(tax_species_silva_1, "Data/tax_species_final_Ascaris_Quanti.rds")
saveRDS(taxid1, "Data/tax_species_decipher_Ascaris_Quanti.rds")

saveRDS(tax_silva_2, "Data/tax_final_Ascaris_Main.rds")
saveRDS(tax_species_silva_2, "Data/tax_species_final_Ascaris_Main.rds")
saveRDS(taxid2, "Data/tax_species_decipher_Ascaris_Main.rds")

##Create an Taxa matrix with ASV as rows and taxonomic level as columns (for Alessio)
taxamat1 <- tax_species_silva_1 # Removing sequence rownames for display only
rownames(taxamat1) <- NULL
rownames(taxamat1) <- paste0("ASV", 1:nrow(taxamat1))
head(taxamat1)
write.csv(taxamat1, "/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Quanti_taxa.csv")

taxamat2 <- tax_species_silva_2 # Removing sequence rown ames for display only
rownames(taxamat2) <- NULL
rownames(taxamat2) <- paste0("ASV", 1:nrow(taxamat2))
head(taxamat2)
write.csv(taxamat2, "/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Main_taxa.csv")

}

##Create a Phylogenetic tree
##Run 1
Align16S<- AlignSeqs(dna1, anchor= NA, verbose= FALSE,  processors= 30) ##Alignment

require(phangorn)
phangAlign16S <- phyDat(as(Align16S, "matrix"), type="DNA")
dm16S <- dist.ml(phangAlign16S) ## Distance matrix
treeNJ16S <- NJ(dm16S) # Note, tip order != sequence order
plot(treeNJ16S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= TRUE) ##Neighbour-Joining tree
fit1 <- pml(treeNJ16S, data=phangAlign16S)
fitGTR16S <- update(fit1, k=4, inv=0.2)
fitGTR16S <- optim.pml(fitGTR16S, model="GTR", optInv=TRUE, optGamma=TRUE,
                       rearrangement = "stochastic", control = pml.control(trace = 0))
plot(fitGTR16S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)

tree1<- fitGTR16S$tree

##Run 2
Align16S<- AlignSeqs(dna2, anchor= NA, verbose= FALSE,  processors= 30) ##Alignment

phangAlign16S <- phyDat(as(Align16S, "matrix"), type="DNA")
dm16S <- dist.ml(phangAlign16S) ## Distance matrix
treeNJ16S <- NJ(dm16S) # Note, tip order != sequence order
plot(treeNJ16S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= TRUE) ##Neighbour-Joining tree
fit1 <- pml(treeNJ16S, data=phangAlign16S)
fitGTR16S <- update(fit1, k=4, inv=0.2)
fitGTR16S <- optim.pml(fitGTR16S, model="GTR", optInv=TRUE, optGamma=TRUE, multicore=TRUE,
                       rearrangement = "stochastic", control = pml.control(trace = 0))
plot(fitGTR16S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)

tree2<- fitGTR16S$tree

##Compile phyloseq
if(Phylobj){
  ##To phyloseq
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  
  ##Load matrices Run 1
  asvmat<- read.csv("/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Quanti_ASV_matrix.csv")
  rownames(asvmat)<-asvmat$X
  asvmat$X<- NULL
  taxamat<- read.csv("/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Quanti_taxa.csv")
  rownames(taxamat)<-taxamat$X
  taxamat$X<- NULL
  dna<- readDNAStringSet("/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Quanti_ASV.fasta")
  ##Load sample data
  sample <- read.csv("Data/Pig_Ascaris_16S_Samples_P2.csv", dec=",", stringsAsFactors=FALSE)
  ##Change in System column. Pigs now are numbered from Pig1 to Pig14 
  ##Pig1 to Pig9 correspond to Experiment 1
  ##Pig10 to Pig14 correspond to Experiment 2 (Original names were Pig1, Pig2, Pig3, Pig4 and Pig5 but changed to avoid confusion)
  
  ##Add sample names 
  row.names(sample)<- sample$Barcode_name
  
  ##To make Phyloseq object
  ##1) Use the ASV matrix and transform it to "OTU table" format
  asv<- otu_table(asvmat, taxa_are_rows = T)
  sample_names(asv)
  ##2) Use sample dataframe and transform it to "sample data" format
  sample<- sample_data(sample)
  sample_names(sample)
  ##3) Use taxa matrix and transform it to "tax table" format
  tax<-tax_table(as.matrix(taxamat))
  sample_names(tax)
  
  PS.1 <- merge_phyloseq(asv, tax)
  
  ###Add phylogenetic tree
  PS.1 <- merge_phyloseq(asv, sample, tax, tree1)
  
  ##Load matrices Run 2
  asvmat<- read.csv("/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Main_ASV_matrix.csv")
  rownames(asvmat)<-asvmat$X
  asvmat$X<- NULL
  taxamat<- read.csv("/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Main_taxa.csv")
  rownames(taxamat)<-taxamat$X
  taxamat$X<- NULL
  dna<- readDNAStringSet("/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Main_ASV.fasta")
  
  ##To make Phyloseq object
  ##1) Use the ASV matrix and transform it to "OTU table" format
  asv<- otu_table(asvmat, taxa_are_rows = T)
  sample_names(asv)
  ##2) Use sample dataframe and transform it to "sample data" format
  sample<- sample_data(sample)
  sample_names(sample)
  ##3) Use taxa matrix and transform it to "tax table" format
  tax<-tax_table(as.matrix(taxamat))
  sample_names(tax)
  
  ###Add phylogenetic tree
  PS.2 <- merge_phyloseq(asv, sample, tax, tree2)
  
  table(sample$System, sample$Compartment, sample$Origin) ## ---> README sample overview (previous filtering)
  
  ##Merge both runs 
  #PS <- merge_phyloseq(PS.1, PS.2) 
  
  #keep<- rownames(seqtab.nochim)
  ### Sample data includes those that didn't worked, so let's eliminate them 
  #samdata <- samdata[samdata$BeGenDiv_Name %in% keep, ]
  #rownames(samdata) <- samdata$sample_names
  
    
  saveRDS(PS.1, file="/fast/AG_Forslund/Victor/data/Ascaris/tmp/PhyloSeq_Ascaris_Quanti.Rds")
  saveRDS(PS.2, file="/fast/AG_Forslund/Victor/data/Ascaris/tmp/PhyloSeq_Ascaris_Main.Rds")
  #saveRDS(sample, file="/SAN/Victors_playground/Ascaris_Microbiome/output/sample.Rds")
}

rm(list = ls())