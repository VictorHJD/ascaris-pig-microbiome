##Project: Ascaris - Pig Microbiome
##Aim: Sample and taxa preprocessing
##Author: Víctor Hugo Jarquín-Díaz
##Root repo setwd("../Ascaris/ascaris/")

##Load libraries 
library(iNEXT)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(data.table)
library(decontam); packageVersion("decontam")

reRun <- FALSE

##Load data 
if(!exists("PS.2")){
  if(isTRUE(reRun)){
    source("R/0_Sequence_preprocessing.R") ## Run the script at base directory of repository!   
  } else {
    PS<- readRDS(file = "/fast/AG_Forslund/Victor/data/Ascaris/tmp/PhyloSeq_Ascaris_Main.Rds") ##Work just with full run for now
    #sample <- read.csv("Data/Pig_Ascaris_16S_Samples_P2.csv", dec=",", stringsAsFactors=FALSE)
  }
}

##Check phyloseq 
summarize_phyloseq(PS)

## HOW many ASVs for off-target eukaryotes and archaea
table(tax_table(PS)[, "Kingdom"], exclude = NULL) ## ---> README results summary

## HOW many reads for off-target eukaryotes and archaea
by((otu_table(PS)), tax_table(PS)[, "Kingdom"], sum) ## --->  README results summary

###Check how many reads have every superkingdom
###Raw counts 
rawcounts <- data.frame(colSums(otu_table(PS)))
rawcounts[,2] <- rownames(rawcounts)
###Bacterial and Eukaryotic counts
rawcounts[,3] <- as.data.frame(colSums(otu_table(subset_taxa(PS, Kingdom%in%"Bacteria"))))
rawcounts[,4] <- as.data.frame(colSums(otu_table(subset_taxa(PS, Kingdom%in%"Eukaryota"))))
rawcounts[,5] <- as.data.frame(colSums(otu_table(subset_taxa(PS, Kingdom%in%"Archaea"))))
rawcounts[,6] <- as.data.frame(colSums(otu_table(subset_taxa(PS, Kingdom%in%NA))))
colnames(rawcounts) <- c("Raw_counts", "Barcode_name", "Bacteria_reads", "Eukaryota_reads", "Archaea_reads", "Unassigned_reads")
#rownames(rawcounts) <- c(1:nrow(rawcounts))
rawcounts <- data.frame(Barcode_name = rawcounts$Barcode_name, 
                             Raw_counts = rawcounts$Raw_counts,
                             Bacteria_reads= rawcounts$Bacteria_reads,
                             Eukaryota_reads= rawcounts$Eukaryota_reads,
                        Archaea_reads= rawcounts$Archaea_reads,
                        Unassigned_reads= rawcounts$Unassigned_reads) 

summary(rawcounts$Raw_counts)
sum(rawcounts$Raw_counts)

as.data.frame(PS@sam_data)->tmp

rawcounts<- left_join(rawcounts, tmp, by= "Barcode_name")

rawcounts <- rawcounts[order(rawcounts$Raw_counts),]
rawcounts$Index <- seq(nrow(rawcounts))

ggplot(data=rawcounts, aes(x=Index, y=Raw_counts, color=Origin)) + 
  geom_point()

###Before further filtering first check rarefaction curves
##species richness (q = 0), Shannon diversity (q = 1, the exponential of Shannon entropy)
##iNEXT uses the observed sample of abundance or incidence data 
##to compute diversity estimates and the associated 95% confidence intervals for the following two types of rarefaction and extrapolation (R/E)
#test<- asvmat[1:5,1:5]
#rare<- iNEXT(test,  q=c(0, 1), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)

#ggiNEXT(rare, type=1, se=TRUE, facet.var="order", color.var="site", grey=FALSE)  

##Filtering 
## 1) Eliminate "empty" samples 
PS <- prune_samples(sample_sums(PS)>0, PS)
summarize_phyloseq(PS)
##2) Sample filtering: Filtering samples with low counts  
PS <- prune_samples(sample_sums(PS)>=2000, PS)
summarize_phyloseq(PS)

##2) Taxa filtering: Remove "uncharachterized" ASVs
PS<- subset_taxa(PS, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) 

##2.1) Check for contaminants
#contamdf.freq <- isContaminant(PS, method="frequency", conc="quant_reading")
#head(contamdf.freq)

##Quick check: After this filtering are the negative controls gone?
plot_richness(PS, x= "Compartment", color = "Compartment" , measures = c("Observed","Chao1", "Shannon")) +
  geom_jitter(alpha= 0.005)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

##3) Eliminate the only negative control out of 6 that was "contaminated" 
PS<- subset_samples(PS, !(Compartment%in%c("Water"))) 

PS.alpha<- PS ##--> for alpha diversity
##3) Supervised Prevalence Filtering: Remove low prevalent taxa
##Create a prevalence dataframe 
Prevdf<- apply(X = otu_table(PS),
               MARGIN = 1,
               FUN = function(x){sum(x > 0)})

##Add taxonomy and total read counts to this data.frame
Prevdf<- data.frame(Prevalence = Prevdf,
                    TotalAbundance = taxa_sums(PS),
                    tax_table(PS))

plyr::ddply(Prevdf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

##Filter the Archaea phylum
phyla2Filter<- c("Euryarchaeota")
PS<- subset_taxa(PS, !Phylum %in% phyla2Filter)

Prevdf<- apply(X = otu_table(PS),
               MARGIN = 1,
               FUN = function(x){sum(x > 0)})
Prevdf<- data.frame(Prevalence = Prevdf,
                    TotalAbundance = taxa_sums(PS),
                    tax_table(PS))

ggplot(Prevdf, aes(TotalAbundance, Prevalence / nsamples(PS),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Log 10 Total Reads") + ylab("Prevalence [Prop. of Samples]") +
  theme_bw()+
  facet_wrap(~Phylum) + theme(legend.position="none")

#4) Remove low prevalent ASVs (This data is used for alpha diversity and differential abundance analysis)
##Remove ASVs that do not show appear more than 5 times in more than 10% of the samples
#wh0 <- genefilter_sample(PS, filterfun_sample(function(x) x > 5), A=0.01*nsamples(PS))
#PS<- prune_taxa(wh0, PS)


##5) Transform to even sampling depth
##5.1) Normalization of proprtions
##Normalization transformation to an even sample size
PS.Norm<- transform_sample_counts(PS, function(x) 1E6 * x/sum(x)) ##--> For beta diversity analysis

##OR 

##5.2) Transform to even sampling depth
## Rarefy without replacement
vegan::rarecurve(t(otu_table(PS)), step=50, cex=0.5)

## Rarefy without replacement to the min sequencing depth 
PS.Rare<- rarefy_even_depth(PS, rngseed=2020, sample.size=min(sample_sums(PS)), replace=F)

##Check how many samples ended after filtering 
table(PS@sam_data$System, PS@sam_data$Compartment) #---> README sample overview (after filtering)

## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case Phylum, Genus and Family)
PS.Fam<-  tax_glom(PS, "Family", NArm = F)

PS.Gen<-  tax_glom(PS, "Genus", NArm = T)

PS.Phy<-  tax_glom(PS, "Phylum", NArm = F)

plot_bar(PS.Phy, fill="Phylum") + 
  facet_wrap(~Compartment, scales= "free_x", nrow=1) +
  theme(legend.position = "none")

##Alpha diversity (rarefied)
plot_richness(PS.Rare, x= "Compartment", color = "Compartment" , measures = c("Observed","Chao1", "Shannon")) +
  geom_jitter(alpha= 0.005)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

##Estimate alpha diversity for individual samples (before merging by Pig)
alphadiv<- estimate_richness(PS)

alphadiv.rare<- estimate_richness(PS.Rare)
as.data.frame(PS@sam_data)->tmp
alphadiv.rare<-cbind(alphadiv.rare, tmp)

##Add sample data into a single data frame 
as.data.frame(PS@sam_data)->tmp

alphadiv<-cbind(alphadiv, tmp)

table(alphadiv$System, alphadiv$Compartment) ## ---> README sample overview (post filtering)

##Prune samples for different questions
##Pig samples (Just GI compartment)
PS.pig<- subset_samples(PS, !(Compartment%in%c("Mock_community", "Ascaris")))
sdt.pig <- data.table(as(sample_data(PS.Pig), "data.frame"), keep.rownames = T)

##Merge compartment data
##Pig samples 
sample_data(PS.Pig)$Replicates<- paste(sdt.pig$System, sdt.pig$Compartment, sep = ".")
sdt.pig$Replicate<- paste(sdt.pig$System, sdt.pig$Compartment, sep = ".")

sdt.pig%>%
  select(Origin, InfectionStatus, AnimalSpecies, WormSex, Compartment,System, Replicate)%>%
  distinct()->sdt.pig
sdt.pig<- sample_data(sdt.pig)
sample_names(sdt.pig) <- sdt.pig$Replicate

PS.pig<-merge_samples(PS.Pig, "Replicates")
sample_data(PS.pig)<- sdt.pig

alphadiv.pig<- estimate_richness(PS.pig) ###Estimate alpha diversity values
alphadiv.pig$Library_Size<- sample_sums(PS.pig)
vegan::rarecurve(otu_table(PS.pig), step=50, cex=0.5)

##Add sample data into a single data frame 
as.data.frame(PS.pig@sam_data)->tmp

alphadiv.pig<-cbind(alphadiv.pig, tmp)

table(alphadiv.pig$System, alphadiv.pig$Compartment) 

###Rarefy and estimate alpha diversity 
PS.pig.rare<- rarefy_even_depth(PS.pig, rngseed=2020, sample.size=min(sample_sums(PS.pig)), replace=T)
vegan::rarecurve(otu_table(PS.pig.rare), step=50, cex=0.5)

alphadiv.pig.rare<- estimate_richness(PS.pig.rare) ###Estimate alpha diversity values
##Add sample data into a single data frame 
as.data.frame(PS.pig.rare@sam_data)->tmp
alphadiv.pig.rare<-cbind(alphadiv.pig.rare, tmp)
alphadiv.pig.rare$Library_Size<- sample_sums(PS.pig)

##Using method from metagenomeSeq
#library("metagenomeSeq")
##Creating a MRexperiment
#asv<- as.matrix(t(PS.pig@otu_table))
#tax<- as.data.frame(PS.pig@tax_table)
#tax<- AnnotatedDataFrame(tax)
#sam<- as.data.frame(PS.pig@sam_data)
#sam <- AnnotatedDataFrame(sam)
#MRex.pig<- newMRexperiment(asv, phenoData=sam, featureData=tax)

##Normalize using Wrench
#cond<- PS.pig@sam_data$InfectionStatus
#MRex.pig <- wrenchNorm(MRex.pig, condition = cond)

#alphanorm <- MRcounts(MRex.pig, norm = TRUE, log = TRUE)
#write.table(alphanorm, "/fast/AG_Forslund/Victor/data/Ascaris/qiime2_run/alphadata_pig_normalised.txt", sep = "\t")
#First colum should be called "#OTU ID"

##Decontam pipeline 
##Inspect library size
#df <- as.data.frame(sample_data(PS.pig)) # Put sample_data into a ggplot-friendly data.frame
#df$LibrarySize <- sample_sums(PS.pig)
#df <- df[order(df$LibrarySize),]
#df$Index <- seq(nrow(df))
#ggplot(data=df, aes(x=Index, y=LibrarySize, color=Compartment)) + geom_point()

##Identify contaminants 
#contamdf.freq <- isContaminant(PS.pig, method="frequency", conc="Compartment")
#head(contamdf.freq)

####Normalization transformation to an even sample size
PS.pig.Norm<- transform_sample_counts(PS.pig, function(x) 1E6 * x/sum(x)) ##--> For beta diversity analysis

##Pig samples compartment and Ascaris (not SH)
PS.PA<- subset_samples(PS, !(System%in%c("SH", "Positive")))
sdt.PA <- data.table(as(sample_data(PS.PA), "data.frame"), keep.rownames = T)
sdt.PA$Replicate<- paste(sdt.PA$System, sdt.PA$Compartment, sep = ".")

##Get the worm ID to indicate the number of worm per individual to avoid loosing samples as "replicates"
wormID<- read.csv("Data/Worm_ID.csv")
##Select samples in sdt.PA
keep <- sdt.PA$rn

wormID%>%
  filter(Barcode_name%in%keep)%>%
  select(Barcode_name, Worm_ID)%>%
  right_join(sdt.PA, by = "Barcode_name")%>%
  mutate(Replicate= paste(Replicate, Worm_ID, sep = "."))%>%
  mutate(Replicate= gsub(".NA", "", Replicate))%>%
  arrange(Barcode_name)-> sdt.PA ##Make sure that is in the right order
  
##Merge compartment and Ascaris data
##Add replicate to sample data 
sample_data(PS.PA)$Replicates<- sdt.PA$Replicate

sdt.PA%>%
  select(Origin, InfectionStatus, AnimalSpecies, WormSex, Compartment,System, Replicate)%>%
  distinct()->sdt.PA

sdt.PA<- sample_data(sdt.PA)
sample_names(sdt.PA) <- sdt.PA$Replicate

PS.PA<-merge_samples(PS.PA, "Replicates")
sample_data(PS.PA)<- sdt.PA

table(sample_data(PS.PA)$System, sample_data(PS.PA)$Compartment)

alphadiv.PA<- estimate_richness(PS.PA) ###Estimate alpha diversity values

##Add sample data into a single data frame 
as.data.frame(PS.PA@sam_data)->tmp

alphadiv.PA<-cbind(alphadiv.PA, tmp)

table(alphadiv.PA$System, alphadiv.PA$Compartment) 

alphadiv.PA$Library_Size<- sample_sums(PS.PA)
vegan::rarecurve(otu_table(PS.PA), step=50, cex=0.5)

###Rarefy and estimate alpha diversity 
PS.PA.rare<- rarefy_even_depth(PS.PA, rngseed=2020, sample.size=min(sample_sums(PS.PA)), replace=T)
vegan::rarecurve(otu_table(PS.PA.rare), step=50, cex=0.5)

alphadiv.PA.rare<- estimate_richness(PS.PA.rare) ###Estimate alpha diversity values
##Add sample data into a single data frame 
as.data.frame(PS.PA.rare@sam_data)->tmp
alphadiv.PA.rare<-cbind(alphadiv.PA.rare, tmp)

####Normalization transformation to an even sample size
PS.PA.Norm<- transform_sample_counts(PS.PA, function(x) 1E6 * x/sum(x)) ##--> For beta diversity analysis

##Ascaris samples
PS.Asc<- subset_samples(PS, Compartment%in%c("Ascaris"))
sdt.Asc <- data.table(as(sample_data(PS.Asc), "data.frame"), keep.rownames = T)
sdt.Asc$Replicate<- paste(sdt.Asc$System, sdt.Asc$Compartment, sep = ".")

##Select samples in sdt.PA
keep <- sdt.Asc$rn

wormID%>%
  filter(Barcode_name%in%keep)%>%
  select(Barcode_name, Worm_ID)%>%
  right_join(sdt.Asc, by = "Barcode_name")%>%
  mutate(Replicate= paste(Replicate, Worm_ID, sep = "."))%>%
  mutate(Replicate= gsub(".NA", "", Replicate))%>%
  arrange(Barcode_name)-> sdt.Asc ##Make sure that is in the right order

##Add replicate to sample data 
sample_data(PS.Asc)$Replicates<- sdt.Asc$Replicate

sdt.Asc%>%
  select(Origin, InfectionStatus, AnimalSpecies, WormSex, Compartment,System, Replicate)%>%
  distinct()->sdt.Asc

sdt.Asc<- sample_data(sdt.Asc)
sample_names(sdt.Asc) <- sdt.Asc$Replicate

PS.Asc<-merge_samples(PS.Asc, "Replicates")
sample_data(PS.Asc)<- sdt.Asc

table(sample_data(PS.Asc)$System, sample_data(PS.Asc)$Compartment)

alphadiv.Asc<- estimate_richness(PS.Asc) ###Estimate alpha diversity values

##Add sample data into a single data frame 
as.data.frame(PS.Asc@sam_data)->tmp

alphadiv.Asc<-cbind(alphadiv.Asc, tmp)

table(alphadiv.Asc$System, alphadiv.Asc$Compartment) 

####Normalization transformation to an even sample size
PS.Asc.Norm<- transform_sample_counts(PS.Asc, function(x) 1E6 * x/sum(x)) ##--> For beta diversity analysis

##Store relevant files
##Phyloseq with all separated samples filtrated 
saveRDS(PS.alpha, "/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.alpha.Rds") ##Row alpha diversity 
saveRDS(PS.Norm, "/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Norm.Rds") ##All not merged samples but normalized for beta diversity
saveRDS(PS.pig, "/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.pig.Rds") ## Data just merged pigs not normalized for alpha diversity plots  
saveRDS(PS.pig.Norm, "/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.pig.Norm.Rds") ## Data just merged pigs normalized for beta diversity plots 
saveRDS(PS.Asc, "/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Asc.Rds") ## Data all Ascaris not normalized for alpha diversity plots 
saveRDS(PS.Asc.Norm, "/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Asc.Norm.Rds") ## Data all Ascaris normalized for beta diversity plots 
saveRDS(PS.PA, "/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Rds") ## Data merged pigs and Ascaris (not SH) not normalized for alpha diversity plots 
saveRDS(PS.PA.Norm, "/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Norm.Rds") ## Data merged pigs and Ascaris (not SH) normalized for beta diversity plots 

##Alpha diverisity tables
saveRDS(alphadiv, "Data/alphadiv.rds")
saveRDS(alphadiv.pig, "Data/alphadiv.pig.rds")
saveRDS(alphadiv.Asc, "Data/alphadiv.Asc.rds")
saveRDS(alphadiv.PA, "Data/alphadiv.PA.rds")

##Create biom format object for PICRUSt2
require("biomformat")
asvmat<- t(as.matrix(PS.PA@otu_table))
biom.PA<- make_biom(asvmat, matrix_element_type = "int")
write_biom(biom.PA,"/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST//biom_PA.biom") #-> For Picrust2

##Select sequences from the ASV in PS.PA
require("DECIPHER")
dna<- readDNAStringSet("/fast/AG_Forslund/Victor/data/Ascaris/tmp/Ascaris_Main_ASV.fasta")
keep <- data.frame(name = rownames(asvmat))
names(dna)
dna.PA<- dna[keep$name]
writeXStringSet(dna.PA, "/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/ASV_PA.fasta") #-> For Picrust2

##Create biom format object for PICRUSt2
asvmat<- t(as.matrix(PS.Asc@otu_table))
biom.Asc<- make_biom(asvmat, matrix_element_type = "int")
write_biom(biom.Asc,"/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/biom_Asc.biom") #-> For Picrust2

##Select sequences from the ASV in PS.PA
keep <- data.frame(name = rownames(asvmat))
names(dna)
dna.Asc<- dna[keep$name]
writeXStringSet(dna.Asc, "/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/ASV_Asc.fasta") #-> For Picrust2
