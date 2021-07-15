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

reRun <- FALSE

##Load data 
if(!exists("PS.2")){
  if(isTRUE(reRun)){
    source("R/0_Sequence_preprocessing.R") ## Run the script at base directory of repository!   
  } else {
    PS<- readRDS(file = "/fast/AG_Forslund/Victor/data/Ascaris/tmp/PhyloSeq_Ascaris_Main.Rds") ##Work just with full run for now
    sample <- read.csv("Data/Pig_Ascaris_16S_Samples_P2.csv", dec=",", stringsAsFactors=FALSE)
  }
}

##Check phyloseq 
summarize_phyloseq(PS)

## HOW many ASVs for off-target eukaryotes and archaea
table(tax_table(PS)[, "Kingdom"], exclude = NULL) ## ---> README results summary

## HOW many reads for off-target eukaryotes and archaea
by((otu_table(PS)), tax_table(PS)[, "Kingdom"], sum) ## --->  README results summary

###Before further filtering first check rarefaction curves
##species richness (q = 0), Shannon diversity (q = 1, the exponential of Shannon entropy)
##iNEXT uses the observed sample of abundance or incidence data 
##to compute diversity estimates and the associated 95% confidence intervals for the following two types of rarefaction and extrapolation (R/E)
rare<- iNEXT(asvmat2,  q=c(0, 1), datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)


##Filtering 
## 1) Eliminate "empty" samples 
PS <- prune_samples(sample_sums(PS)>0, PS)

##2) Sample filtering: Filtering samples with low counts  
PS <- prune_samples(sample_sums(PS)>=2000, PS)
summarize_phyloseq(PS)

##2) Taxa filtering: Remove "uncharachterized" ASVs
PS<- subset_taxa(PS, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) ##--> for alpha diversity

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

phyla2Filter<- c("Aquificota", "Dependentiae")
PS<- subset_taxa(PS3, !Phylum %in% phyla2Filter)

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
wh0 <- genefilter_sample(PS, filterfun_sample(function(x) x > 5), A=0.01*nsamples(PS))
PS<- prune_taxa(wh0, PS)
hist(readcount(PS))

##4) Transform to even sampling depth
## Rarefy without replacement to the min sequencing depth 
#PS4<- rarefy_even_depth(PS3, rngseed=2020, sample.size=min(sample_sums(PS3)), replace=F)
##Normalization transformation to an even sample size
PS4<- transform_sample_counts(PS3, function(x) 1E6 * x/sum(x))


##4) Transform to even sampling depth
## Rarefy without replacement
vegan::rarecurve(t(otu_table(PS)), step=50, cex=0.5)

PS<- rarefy_even_depth(PS, rngseed=2020, sample.size=min(sample_sums(PS)), replace=F)
readcount(PS)

saveRDS(PS, file="/SAN/Victors_playground/Ascaris_Microbiome/output/PhyloSeqRare.Rds")
## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case Phylum and Family)
PS.Fam<-  tax_glom(PS, "Family", NArm = F)
summarize_phyloseq(PS.Fam)

PS.Gen<-  tax_glom(PS, "Genus", NArm = T)
summarize_phyloseq(PS.Gen)

PS.Phy<-  tax_glom(PS, "Phylum", NArm = F)
summarize_phyloseq(PS.Phy)

plot_bar(PS.Phy, fill="Phylum") + facet_wrap(~Compartment, scales= "free_x", nrow=1)

##Alpha diversity (rarefied)
plot_richness(PS, x= "Compartment", color = "Compartment" , measures = c("Observed","Chao1", "Shannon")) +
  #geom_boxplot()+
  geom_jitter(alpha= 0.005)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

alphadiv<- estimate_richness(PS)

alphadiv%>%
  rownames_to_column()->tmp1

as.tibble(sample)%>%
  mutate(rowname= paste0("Sample", 1:nrow(sample)))->tmp2

tmp1<-inner_join(tmp1, tmp2, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
alphadiv<- tmp1
rm(tmp1,tmp2)

table(alphadiv$System, alphadiv$Compartment) ## ---> README sample overview (post filtering)

##Prune samples for different questions
##Pig samples (Just GI compartment)
PS.Pig<- subset_samples(PS, !(Compartment%in%c("Negative", "Ascaris")))
sdt.pig <- data.table(as(sample_data(PS.Pig), "data.frame"), keep.rownames = T)
alphadiv.pig<- estimate_richness(PS.Pig) ###Estimate alpha diversity values 
alphadiv.pig%>%
  rownames_to_column()->tmp1

row.names(sdt.pig)<- sdt.pig$rn
names(sdt.pig)[names(sdt.pig) == "rn"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.pig, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.pig<- tmp1
rm(tmp1,alphadiv.pig)

##Merge compartment data
##Pig samples (no faeces) 
sample_data(PS.Pig)$Replicates<- paste(sdt.pig$System, sdt.pig$Compartment, sep = ".")
sdt.pig$Replicate<- paste(sdt.pig$System, sdt.pig$Compartment, sep = ".")
sdt.pig%>%
  select(InfectionStatus,AnimalSpecies,WormSex,Live,Compartment,System, Replicate)%>%
  distinct()->sdt.pig2
sdt.pig2<- sample_data(sdt.pig2)
sample_names(sdt.pig2) <- sdt.pig2$Replicate

PS.pig2<-merge_samples(PS3.Pig, "Replicates")
sample_data(PS.pig2)<- sdt.pig2

alphadiv.pig2<- estimate_richness(PS.pig2) ###Estimate alpha diversity values
alphadiv.pig2%>%
  rownames_to_column()->tmp1

row.names(sdt.pig2)<- sdt.pig2$Replicate
names(sdt.pig2)[names(sdt.pig2) == "Replicate"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.pig2, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.pig2<- tmp1
rm(tmp1,alphadiv.pig2)


##Pig samples compartment and Ascaris (not SH)
PS.PA<- subset_samples(PS, !(Compartment%in%c("Negative", "Faeces")))
PS.PA<- subset_samples(PS.PA, !(System%in%c("SH")))
sdt.PA <- data.table(as(sample_data(PS.PA), "data.frame"), keep.rownames = T)
alphadiv.PA<- estimate_richness(PS.PA) ###Estimate alpha diversity values
alphadiv.PA%>%
  rownames_to_column()->tmp1

row.names(sdt.PA)<- sdt.PA$rn
names(sdt.PA)[names(sdt.PA) == "rn"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.PA, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.PA<- tmp1
rm(tmp1,alphadiv.PA)

##Merge compartment and Ascaris data
sample_data(PS.PA)$Replicates<- paste(sdt.PA$System, sdt.PA$Compartment, sep = ".")
sdt.PA$Replicate<- paste(sdt.PA$System, sdt.PA$Compartment, sep = ".")
sdt.PA%>%
  select(InfectionStatus,AnimalSpecies,Live,Compartment,System, Replicate)%>%
  distinct()->sdt.PA2
sdt.PA2$InfectionStatus[is.na(sdt.PA2$InfectionStatus)] <- "Worm" ##Remove NAs
sdt.PA2<- sample_data(sdt.PA2)
sample_names(sdt.PA2) <- sdt.PA2$Replicate

PS.PA2<-merge_samples(PS.PA, "Replicates")
sample_data(PS.PA2)<- sdt.PA2

alphadiv.PA2<- estimate_richness(PS.PA2) ###Estimate alpha diversity values
alphadiv.PA2%>%
  rownames_to_column()->tmp1

row.names(sdt.PA2)<- sdt.PA2$Replicate
names(sdt.PA2)[names(sdt.PA2) == "Replicate"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.PA2, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.PA2<- tmp1
rm(tmp1,alphadiv.PA2)

##Ascaris samples
PS.Asc<- subset_samples(PS, Compartment%in%c("Ascaris"))
sdt.Asc <- data.table(as(sample_data(PS.Asc), "data.frame"), keep.rownames = T)
alphadiv.Asc<- estimate_richness(PS.Asc) ###Estimate alpha diversity values
alphadiv.Asc%>%
  rownames_to_column()->tmp1

row.names(sdt.Asc)<- sdt.Asc$rn
names(sdt.Asc)[names(sdt.Asc) == "rn"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.Asc, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.Asc<- tmp1
rm(tmp1,alphadiv.Asc)