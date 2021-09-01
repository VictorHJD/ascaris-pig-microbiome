##Project: Ascaris - Pig Microbiome
##Aim: Differential abundant taxa
##Author: Víctor Hugo Jarquín-Díaz
##Root repo setwd("../Ascaris/ascaris/")

##Load packages
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(microbiome)
library(vegan)
library(picante)
library(ALDEx2)
library(metagenomeSeq)
library(HMP)
library(dendextend)
library(selbal)
library(rms)
library(breakaway)

library(data.table)
require(ggpubr)
require(RColorBrewer)
require(rstatix)
library(cowplot)
library(ggsci)

##For taxa 
tax.palette<- c("Taxa less represented" = "black",  "Unassigned"="lightgray", "Streptococcus"= "#925E9FFF", 
                "Lactobacillus"=  "#631879FF", "Clostridium sensu stricto 1"= "#00468BFF","Bifidobacterium" = "#3C5488FF",
                "Roseburia" = "#0072B5FF", "Escherichia-Shigella"= "#E762D7FF", "Pseudomonas" = "#ED0000FF",
                "Agathobacter" = "#0099B4FF" , "Prevotella" = "#E64B35FF", "Veillonella"= "#3B4992FF", 
                "Turicibacter" = "#AD002AFF", "Terrisporobacter"  = "#00A087FF",  "Romboutsia" = "#F39B7FFF", 
                "Prevotellaceae NK3B31 group"= "#8491B4FF", "Megasphaera"= "#CD534CFF", "Anaerovibrio" = "#FAFD7CFF",
                "Intestinibacter"="#6F99ADFF", "Parasutterella" = "#BC3C29FF", "Helicobacter"  = "#E18727FF",
                "Succinivibrio"= "#7876B1FF", "Prevotellaceae UCG-003" = "#FFDC91FF", 
                "Klebsiella"  = "#EE4C97FF",  "Aliterella"= "#42B540FF", "Succinivibrionaceae UCG-001" = "#925E9FFF", 
                "Erysipelotrichaceae UCG-002"= "#008B45FF", "Alloprevotella"= "#4DBBD5FF",  "Olsenella"= "#B09C85FF", 
                "Pseudoscardovia"= "#BB0021FF", "Aeromonas"= "#FFCD00FF")

##Color palette for compartment and system ##

pal.compartment <- c("Ascaris"="#1B9E77","Cecum"= "#D95F02","Colon"= "#7570B3",
                     "Duodenum"= "#E7298A","Ileum"= "#66A61E","Jejunum"="#E6AB02")

pal.system <- c("Pig1"= "#A6761D","Pig2"= "#666666","Pig3"= "#A6CEE3","Pig4"= "#1F78B4",
                "Pig5"= "#B2DF8A","Pig6"= "#33A02C","Pig7"= "#FB9A99","Pig8"="#E31A1C","Pig9"= "#FDBF6F",
                "Pig10"= "#FF7F00","Pig11"= "#CAB2D6","Pig12"= "#6A3D9A","Pig13"= "#FFFF99",  "Pig14"= "#3B3B3BFF", "SH" = "#BB0021FF")

##Load data 
##General
PS.alpha<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.alpha.Rds") ##Row alpha diversity 
PS.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Norm.Rds") ##All not merged samples but normalized for beta diversity

##Question specific
PS.pig<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.pig.Rds") ## Data just merged pigs not normalized for alpha diversity plots  
PS.pig.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.pig.Norm.Rds") ## Data just merged pigs normalized for beta diversity plots 
PS.Asc<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Asc.Rds") ## Data all Ascaris not normalized for alpha diversity plots 
PS.Asc.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Asc.Norm.Rds") ## Data all Ascaris normalized for beta diversity plots 
PS.PA<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Rds") ## Data merged pigs and Ascaris (not SH) not normalized for alpha diversity plots 
PS.PA.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Norm.Rds") ## Data merged pigs and Ascaris (not SH) normalized for beta diversity plots 

##Alpha diverisity tables with sample information
alphadiv<- readRDS("Data/alphadiv.rds")
alphadiv.pig<- readRDS("Data/alphadiv.pig.rds")
alphadiv.Asc<- readRDS("Data/alphadiv.Asc.rds")
alphadiv.PA<- readRDS("Data/alphadiv.PA.rds")

#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ InfectionStatus, data = df)
}

wilcox_pval <- function(df){
  wilcox.test(abund ~ InfectionStatus, data = df)$p.value
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#Generate data.frame with ASVs and metadata
PS.Rel<- transform_sample_counts(PS.pig, function(x) 100 * x/sum(x)) ##--> Adjust Rel abund

stat.pig <- data.frame(data.frame(phyloseq::otu_table(PS.Rel)))
stat.pig$InfectionStatus <- phyloseq::sample_data(PS.Rel)$InfectionStatus

#Create nested data frames by ASV and loop over each using map 
stat.pig %>%
  gather(key = ASV, value = abund, -InfectionStatus) %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))-> wilcox.pig                      

#Unnesting
wilcox.pig %>%
  dplyr::select(ASV, p_value) %>%
  unnest()-> wilcox.pig 

#Adding taxonomic information
taxa.info <- data.frame(tax_table(PS.Rel))
taxa.info %>% 
  rownames_to_column(var = "ASV")-> taxa.info

#Computing bonferroni and FDR corrected p-values
wilcox.pig%>%
  full_join(taxa.info) %>%
  arrange(p_value) %>%
  ungroup()%>%
  dplyr::mutate(BH_FDR= rstatix::adjust_pvalue(p_value, method = "BH")) %>%
  dplyr::mutate(BON_Adj= rstatix::adjust_pvalue(p_value, method ="bonferroni")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(ASV, p_value, BH_FDR, BON_Adj, everything())-> wilcox.pig


##Analysis for Pigs
##Subset just site of infection
tmp<- row.names(PS.pig.Norm@sam_data)
tmp<- alphadiv.pig[rownames(alphadiv.pig)%in%tmp, ]

tmp%>%
  dplyr::filter(Compartment== "Jejunum")%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.pig.Jej<- subset_samples(PS.pig.Norm, Replicate%in%Inf.Keep)

DS.pig.Jej <- phyloseq_to_deseq2(PS.pig.Jej, ~InfectionStatus)

geoMeans <- apply(counts(DS.pig.Jej), 1, gm_mean)

DS.pig.Jej <- estimateSizeFactors(DS.pig.Jej , geoMeans = geoMeans)
DS.pig.Jej <- DESeq(DS.pig.Jej)

res <- results(DS.pig.Jej, cooksCutoff = FALSE)

sigtab <- cbind(as(res,"data.frame"), as(tax_table(PS.pig.Jej)[rownames(res), ], "matrix"))
##Remove rows with columns that contain NA
sigtab <- sigtab[complete.cases(sigtab), ]
head(sigtab,20)

##Volcano plot to detect differential genes in Jejunum between Non infected vs Infected
#The significantly deferentially abundant genes are the ones found upper-left and upper-right corners
##Add a column to the data specifying if they are highly (positive) or lowly abundant (negative)
## Considering the comparison Non infected vs Infected

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

##Merge Genus and species 
sigtab%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)-> sigtab

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
require("ggrepel")
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(size=3, alpha= 0.9, aes(fill= AbundLev), shape= 21, color= "black")+
  scale_fill_manual(values=c("#00A087B2", "#800000FF", "#767676FF"), 
                    labels = c("High (Non Infected)", "Low (Non Infected)", "Not Significant"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "A)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Gene \n abundance")+
  theme_bw()+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  geom_text_repel(data = subset(sigtab, AbundLev=="Low"),
                  aes(label = Bacteria_name),
                  size = 3,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = 15)+
  geom_text_repel(data = subset(sigtab, AbundLev=="High"),
                  aes(label = Bacteria_name),
                  size = 3,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"),
                  max.overlaps = 12)+
  theme(text = element_text(size=16))-> A

##Jejunum Infected and Ascaris 
tmp<- row.names(PS.PA.Norm@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  dplyr::filter(InfectionStatus!= ("Non_infected"))%>%
  dplyr::filter(Compartment%in% c("Jejunum", "Ascaris"))%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.PA.Jej<- subset_samples(PS.PA.Norm, Replicate%in%Inf.Keep)

DS.PA.Jej <- phyloseq_to_deseq2(PS.PA.Jej, ~Compartment)

geoMeans <- apply(counts(DS.PA.Jej), 1, gm_mean)

DS.PA.Jej <- estimateSizeFactors(DS.PA.Jej , geoMeans = geoMeans)
DS.PA.Jej <- DESeq(DS.PA.Jej)

res <- results(DS.PA.Jej, cooksCutoff = FALSE)

sigtab <- cbind(as(res,"data.frame"), as(tax_table(PS.PA.Jej)[rownames(res), ], "matrix"))
##Remove rows with columns that contain NA
sigtab <- sigtab[complete.cases(sigtab), ]
head(sigtab,20)

##Volcano plot to detect differential genes in Jejunum between Non infected vs Infected
#The significantly deferentially abundant genes are the ones found upper-left and upper-right corners
##Add a column to the data specifying if they are highly (positive) or lowly abundant (negative)
## Considering the comparison Non infected vs Infected

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

##Merge Genus and species 
sigtab%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)-> sigtab

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
require("ggrepel")
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(size=3, alpha= 0.9, aes(fill= AbundLev), shape= 21, color= "black")+
  scale_fill_manual(values=c("#00A087B2", "#800000FF", "#767676FF"), 
                    labels = c("High (Non Infected)", "Low (Non Infected)", "Not Significant"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "B)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Gene \n abundance")+
  theme_bw()+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  geom_text_repel(data = subset(sigtab, AbundLev=="Low"),
                  aes(label = Bacteria_name),
                  size = 3,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = 15)+
  geom_text_repel(data = subset(sigtab, AbundLev=="High"),
                  aes(label = Bacteria_name),
                  size = 3,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"),
                  max.overlaps = 12)+
  theme(text = element_text(size=16))-> B

##ggsave(file = "results/figures/Q1_Gene_Abundance_black.pdf", plot = D, width = 10, height = 8)
#ggsave(file = "results/figures/Q1_Gene_Abundance_black.png", plot = D, width = 10, height = 8)