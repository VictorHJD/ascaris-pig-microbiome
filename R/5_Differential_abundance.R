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

###DeSeq2 pipeline 
##Analysis for Pigs
##Subset just site of infection
tmp<- row.names(PS.PA.Norm@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  dplyr::filter(Compartment== "Jejunum")%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.pig.Jej<- subset_samples(PS.PA.Norm, Replicate%in%Inf.Keep)

DS.pig.Jej <- phyloseq_to_deseq2(PS.pig.Jej, ~InfectionStatus)

geoMeans <- apply(counts(DS.pig.Jej), 1, gm_mean)

DS.pig.Jej <- estimateSizeFactors(DS.pig.Jej , geoMeans = geoMeans)
DS.pig.Jej <- estimateDispersions(DS.pig.Jej)
DS.pig.Jej <- DESeq(DS.pig.Jej, test = "Wald", fitType= "mean")

res <- results(DS.pig.Jej, cooksCutoff = FALSE)
##Remove rows with columns that contain NA
res <- res[complete.cases(res), ]

sigtab <- cbind(as(res,"data.frame"), as(tax_table(PS.pig.Jej)[rownames(res), ], "matrix"))

head(sigtab,20)

sigtab[complete.cases(sigtab), ]
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
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("UCG-005", "Oscillospiraceae UCG-005", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> sigtab

##Save this data
write.csv(sigtab, "Tables/Q2_DiffAbund_Jej_InfvNonInf.csv")

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
require("ggrepel")
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(size=3, alpha= 0.5, position=position_jitter(0.2), aes(fill= AbundLev), shape= 21, color= "black")+
  scale_fill_manual(values=c("#1a9850", "#800000FF", "#767676FF"), 
                    labels = c("High (Jejunum Non Infected)", "Low (Jejunum Non Infected)", "Not Significant"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "A)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Abundance level")+
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
                  max.overlaps = 15)+
  theme(text = element_text(size=16))-> A

##Alternative 
sigtab%>%
  dplyr::filter(AbundLev!= "NS")%>%
  rownames_to_column()%>%
  ggplot(aes(x = reorder(rowname, -log2FoldChange), y = log2FoldChange)) +
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Family), color= "black") + 
  labs(y = "\nLog2 Fold-Change for Non infected vs Infected", x = "") +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "none") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dotted")

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
DS.PA.Jej <- estimateDispersions(DS.PA.Jej, fitType= "mean")
DS.PA.Jej <- DESeq(DS.PA.Jej, test = "Wald", fitType= "mean")

res <- results(DS.PA.Jej, cooksCutoff = FALSE)
##Remove rows with columns that contain NA
res <- res[complete.cases(res), ]

sigtab <- cbind(as(res,"data.frame"), as(tax_table(PS.PA.Jej)[rownames(res), ], "matrix"))
head(sigtab,20)

##Volcano plot to detect differential genes in Jejunum between Non infected vs Infected
#The significantly deferentially abundant genes are the ones found upper-left and upper-right corners
##Add a column to the data specifying if they are highly (positive) or lowly abundant (negative)
## Considering the comparison Jejunum vs Ascaris

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

##Merge Genus and species 
sigtab%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("UCG-005", "Oscillospiraceae UCG-005", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> sigtab

##Save this data
write.csv(sigtab, "Tables/Q2_DiffAbund_PA_JejvAsc.csv")

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
require("ggrepel")
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(size=3, alpha= 0.5, position=position_jitter(0.2), aes(fill= AbundLev), shape= 21, color= "black")+
  scale_fill_manual(values=c("#00A087B2", "#800000FF", "#767676FF"), 
                    labels = c("High (Jejunum Infected)", "Low (Jejunum Infected)", "Not Significant"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "B)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Abundance level")+
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

##Ascaris FU vs SH 
alphadiv.Asc%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(Location = case_when(Origin %in% c("Experiment_1", "Experiment_2")  ~ "FU",
                                     Origin == "Slaughterhouse" ~ "SH"))%>%
  dplyr::mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))-> alphadiv.Asc

PS.Asc.Norm@sam_data<- sample_data(alphadiv.Asc)

DS.Asc <- phyloseq_to_deseq2(PS.Asc.Norm, ~Location)

geoMeans <- apply(counts(DS.Asc), 1, gm_mean)

DS.Asc <- estimateSizeFactors(DS.Asc, geoMeans = geoMeans)
DS.Asc <- estimateDispersions(DS.Asc, fitType= "mean")
DS.Asc <- DESeq(DS.Asc, test = "Wald", fitType= "mean")

res <- results(DS.Asc, cooksCutoff = FALSE)
##Remove rows with columns that contain NA
res <- res[complete.cases(res), ]

sigtab <- cbind(as(res,"data.frame"), as(tax_table(PS.Asc.Norm)[rownames(res), ], "matrix"))
head(sigtab,20)

##Volcano plot to detect differential genes in Jejunum between Non infected vs Infected
#The significantly deferentially abundant genes are the ones found upper-left and upper-right corners
##Add a column to the data specifying if they are highly (positive) or lowly abundant (negative)
## Considering the comparison Jejunum vs Ascaris

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

##Merge Genus and species 
sigtab%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("UCG-005", "Oscillospiraceae UCG-005", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> sigtab

##Save this data
write.csv(sigtab, "Tables/Q2_DiffAbund_Asc_FUvSH.csv")

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
require("ggrepel")
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(size=3, alpha= 0.5, position=position_jitter(0.2), aes(fill= AbundLev), shape= 21, color= "black")+
  scale_fill_manual(values=c("#00A087B2", "#800000FF", "#767676FF"), 
                    labels = c("High (Ascaris SH)", "Low (Ascaris SH)", "Not Significant"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "C)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Abundance level")+
  theme_bw()+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  geom_text_repel(data = subset(sigtab, AbundLev=="Low"),
                  aes(label = Bacteria_name),
                  size = 3,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = 25)+
  geom_text_repel(data = subset(sigtab, AbundLev=="High"),
                  aes(label = Bacteria_name),
                  size = 3,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"),
                  max.overlaps = 12)+
  theme(text = element_text(size=16))-> C

##Ascaris SH Female vs Male
tmp<- row.names(PS.Asc.Norm@sam_data)
tmp<- alphadiv.Asc[rownames(alphadiv.Asc)%in%tmp, ]

tmp%>%
  dplyr::filter(Location!= ("FU"))%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.Asc.SH<- subset_samples(PS.Asc.Norm, Replicate%in%Inf.Keep)

DS.Asc.SH <- phyloseq_to_deseq2(PS.Asc.SH, ~WormSex)

geoMeans <- apply(counts(DS.Asc.SH), 1, gm_mean)

DS.Asc.SH <- estimateSizeFactors(DS.Asc.SH, geoMeans = geoMeans)
DS.Asc.SH <- estimateDispersions(DS.Asc.SH, fitType= "mean")
DS.Asc.SH <- DESeq(DS.Asc.SH, test = "Wald", fitType= "mean")

res <- results(DS.Asc.SH, cooksCutoff = FALSE)
##Remove rows with columns that contain NA
res <- res[complete.cases(res), ]

sigtab <- cbind(as(res,"data.frame"), as(tax_table(PS.Asc.SH)[rownames(res), ], "matrix"))
head(sigtab,20)

##Volcano plot to detect differential genes in SH Ascaris between Female vs Male
#The significantly deferentially abundant genes are the ones found upper-left and upper-right corners
##Add a column to the data specifying if they are highly (positive) or lowly abundant (negative)
## Considering the comparison Female vs Male

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 0.6 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -0.6 & sigtab$padj < 0.001] <- "Low"

##Merge Genus and species 
sigtab%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("UCG-005", "Oscillospiraceae UCG-005", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> sigtab

##Save this data
write.csv(sigtab, "Tables/Q2_DiffAbund_Asc_FemvMal.csv")

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
require("ggrepel")
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(size=3, alpha= 0.5, position=position_jitter(0.2), aes(fill= AbundLev), shape= 21, color= "black")+
  scale_fill_manual(values=c("#00A087B2", "#800000FF", "#767676FF"), 
                    labels = c("High (Ascaris Male)", "Low (Ascaris Male)", "Not Significant"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "D)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Abundance level")+
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
                  max.overlaps = 15)+
  theme(text = element_text(size=16))-> D

##Save plots 
##ggsave(file = "Figures/Q1_Diff_Abundance_JejInfNonInf.pdf", plot = A, width = 12, height = 8)
#ggsave(file = "Figures/Q1_Diff_Abundance_JejInfNonInf.png", plot = A, width = 12, height = 8)

##ggsave(file = "Figures/Q1_Diff_Abundance_JejAsc.pdf", plot = B, width = 12, height = 8)
#ggsave(file = "Figures/Q1_Diff_Abundance_JejAsc.png", plot = B, width = 12, height = 8)

##ggsave(file = "Figures/Q1_Diff_Abundance_AscFUSH.pdf", plot = C, width = 12, height = 8)
#ggsave(file = "Figures/Q1_Diff_Abundance_AscFUSH.png", plot = C, width = 12, height = 8)

##ggsave(file = "Figures/Q1_Diff_Abundance_AscFM.pdf", plot = D, width = 12, height = 8)
#ggsave(file = "Figures/Q1_Diff_Abundance_AscFM.png", plot = D, width = 12, height = 8)

##Predictions
bray_dist<- phyloseq::distance(PS.pig.Jej, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.pig.Jej,
                      method="PCoA", distance="bray")


tmp<- row.names(PS.pig.Jej@sam_data)
tmp<- alphadiv.pig[rownames(alphadiv.pig)%in%tmp, ]

tmp%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

Jej.adonis<- vegan::adonis(bray_dist~ InfectionStatus + System + Origin,
                                   permutations = 999, data = tmp, na.action = T)

####
## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, tmp$InfectionStatus, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2 
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("InfectionStatus","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

##Add sample data
tmp%>%
  dplyr::select(!c(InfectionStatus))%>%
  cbind(seg.data)-> seg.data

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, shape= InfectionStatus, fill= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F)+
  labs(tag= "A)", fill  = "Infection status", color= "Origin of samples")+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> A

###Compare distances at PCo1 and PCo2 among groups
# Horizontal marginal boxplot - to appear at the top of the chart
##PCo1
seg.data%>%
  wilcox_test(v.PCoA1 ~ InfectionStatus)%>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Beta_Compartment_PCo1_JejInf.csv")

seg.data%>%
  ggplot(aes(x=InfectionStatus, y= v.PCoA1))+
  geom_boxplot(aes(fill= InfectionStatus), color= "black", alpha= 0.5)+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p.signif}", coord.flip = T)+
  coord_flip()+ 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(1, 0.2, -0.5, 0.5), "lines"))-> xplot

# Vertical marginal boxplot - to appear at the right of the chart
##PCo2
seg.data%>%
  wilcox_test(v.PCoA2 ~ InfectionStatus)%>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Beta_Compartment_PCo1_JejInf.csv")

seg.data%>%
  ggplot(aes(x=InfectionStatus, y= v.PCoA2))+
  geom_boxplot(aes(fill= InfectionStatus), color= "black", alpha= 0.5)+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p.signif}")+
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0.2, 1, 0.5, -0.5), "lines"))->  yplot 

require(cowplot)
p1<- insert_xaxis_grob(A, xplot, grid::unit(.25, "null"), position = "top")
p2<- insert_yaxis_grob(p1, yplot, grid::unit(.25, "null"), position = "right")
A2<- ggdraw(p2)

##Save plots 
##ggsave(file = "Figures/Q1_Composition_JejInfNonInf.pdf", plot = A2, width = 10, height = 8)
#ggsave(file = "Figures/Q1_Composition_JejInfNonInf.png", plot = A2, width = 10, height = 8)

##Subset just site of infection and worms
tmp<- row.names(PS.PA.Norm@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  dplyr::filter(Compartment%in% c("Jejunum", "Ascaris"))%>%
  #dplyr::filter(System!= "Pig14")%>% #No Jejunum 
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.PA.Jej<- subset_samples(PS.PA.Norm, Replicate%in%Inf.Keep)

bray_dist<- phyloseq::distance(PS.PA.Jej, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.PA.Jej,
                      method="PCoA", distance="bray")


tmp<- row.names(PS.PA.Jej@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

JejAsc.adonis<- vegan::adonis(bray_dist~ InfectionStatus + System + Origin,
                           permutations = 999, data = tmp, na.action = T)

####
## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, tmp$InfectionStatus, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2 
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("InfectionStatus","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

##Add sample data
tmp%>%
  dplyr::select(!c(InfectionStatus))%>%
  cbind(seg.data)-> seg.data

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, shape= InfectionStatus, fill= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24, 25, 21), labels = c("Infected", "Non infected", "Ascaris"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF", "#fdae61"), labels = c("Infected", "Non infected", "Ascaris"))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25, 21))), shape= F)+
  labs(tag= "B)", fill  = "Infection status", color= "Origin of samples")+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> B

###Compare distances at PCo1 and PCo2 among groups
# Horizontal marginal boxplot - to appear at the top of the chart
##PCo1
seg.data%>%
  wilcox_test(v.PCoA1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Beta_Compartment_PCo1_JejAscInf.csv")

seg.data%>%
  ggplot(aes(x=InfectionStatus, y= v.PCoA1))+
  geom_boxplot(aes(fill= InfectionStatus), color= "black", alpha= 0.5)+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF", "#fdae61"), labels = c("Infected", "Non infected", "Ascaris"))+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p.adj.signif}", coord.flip = T)+
  coord_flip()+ 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(1, 0.2, -0.5, 0.5), "lines"))-> xplot

# Vertical marginal boxplot - to appear at the right of the chart
##PCo2
seg.data%>%
  wilcox_test(v.PCoA2 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Beta_Compartment_PCo2_JejAscInf.csv")

stats.test%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(y.position= c(0.35))-> stats.test

seg.data%>%
  ggplot(aes(x=InfectionStatus, y= v.PCoA2))+
  geom_boxplot(aes(fill= InfectionStatus), color= "black", alpha= 0.5)+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF", "#fdae61"), labels = c("Infected", "Non infected", "Ascaris"))+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p.adj.signif}")+
  theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        legend.position = "none",
        plot.margin = unit(c(0.2, 1, 0.5, -0.5), "lines"))->  yplot 

require(cowplot)
p1<- insert_xaxis_grob(B, xplot, grid::unit(.25, "null"), position = "top")
p2<- insert_yaxis_grob(p1, yplot, grid::unit(.25, "null"), position = "right")
B2<- ggdraw(p2)

##Save plots 
##ggsave(file = "Figures/Q1_Composition_JejInfNonInfAsc.pdf", plot = B2, width = 10, height = 8)
#ggsave(file = "Figures/Q1_Composition_JejInfNonInfAsc.png", plot = B2, width = 10, height = 8)

#Generate data.frame
seg.data$Status <- ifelse(seg.data$InfectionStatus == "Non_infected", 0, 1)
head(seg.data)

#Specify a datadist object (for rms)
seg.data$WormSex<- NULL
dd <- rms::datadist(seg.data)
options(datadist = "dd")

#Plot the unconditional associations
ggplot(seg.data, aes(x = v.PCoA1, y = Status)) +
  Hmisc::histSpikeg(Status ~ v.PCoA1, lowess = TRUE, data = seg.data) +
  labs(x = "\nPCoA1", y = "Pr(Infection status)\n")

ggplot(seg.data, aes(x = v.PCoA2, y = Status)) +
  Hmisc::histSpikeg(Status ~ v.PCoA2, lowess = TRUE, data = seg.data) +
  labs(x = "\nPCoA2", y = "Pr(Infection status)\n")

#Fit full model with splines (3 knots each)
m1 <- rms::lrm(Status ~ rms::rcs(v.PCoA1, 3) + rms::rcs(v.PCoA2, 3), data = seg.data, x = TRUE, y = TRUE)

#Grid search for penalties
rms::pentrace(m1, list(simple = c(0, 1, 2), nonlinear = c(0, 100, 200)))

pen_m1 <- update(m1, penalty = list(simple = 0, nonlinear = 200))

#Obtain optimism corrected estimates
val <- rms::validate(pen_m1)

#Plot calibration
cal <- rms::calibrate(pen_m1, B = 200)
plot(cal) ##Crap model 

##Use selbal
##forward-selection method for the identification of taxa whose relative abundance, or balance, is associated with the response variable 
#we only use 1 repeat of 5-fold cross-validation to tune the selections. 
#To run the function we need to specify two input objects:
#x is a matrix with the microbiome information. It represents the number of counts or reads for each sample (row) and each taxon (column).
#y is a vector with the response variable. It should be specified as a factor if the response variable is dichotomous and numeric if it is continuous.

##Run for Jejunum infected vs non infected failed due to low number of non infected samples :S
##Let's do it for Jejunum infected vs Ascaris 
cv_sebal <- selbal::selbal.cv(x = data.frame(phyloseq::otu_table(PS.PA.Jej)), 
                             y = phyloseq::sample_data(PS.PA.Jej)$InfectionStatus, 
                             n.fold = 5, n.iter = 1, logit.acc = "AUC", zero.rep = "one") 
