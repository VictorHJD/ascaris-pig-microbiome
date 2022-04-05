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
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

###DeSeq2 pipeline 
##Jejunum Infected and Ascaris 
tmp<- row.names(PS.PA@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  dplyr::filter(InfectionStatus!= ("Non_infected"))%>%
  dplyr::filter(Compartment%in% c("Jejunum", "Ascaris"))%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.PA.Jej<- subset_samples(PS.PA, Replicate%in%Inf.Keep)

##Prevalence and total abundance
Prev.JA<- apply(X = t(otu_table(PS.PA.Jej)),
                MARGIN = 1,
                FUN = function(x){sum(x > 0)})
Prev.JA<- data.frame(Prevalence = Prev.JA,
                     TotalAbundance = taxa_sums(PS.PA.Jej),
                     tax_table(PS.PA.Jej))
Prev.JA%>%
  dplyr::select(c(Prevalence,TotalAbundance, Genus))-> Prev.JA

#Prevalence Jejunum
PS.tmp<- subset_samples(PS.PA.Jej, PS.PA.Jej@sam_data$Compartment=="Jejunum")
tmp1<- apply(X = t(otu_table(PS.tmp)),
            MARGIN = 1,
            FUN = function(x){sum(x > 0)})
tmp1<- data.frame(Prevalence_Jej = tmp1,
                 JejAbundance = taxa_sums(PS.tmp),
                 tax_table(PS.tmp))

tmp1%>%
  dplyr::select(c(Prevalence_Jej, JejAbundance))%>%
  cbind(Prev.JA)-> Prev.JA

#Prevalence Ascaris
PS.tmp<- subset_samples(PS.PA.Jej, PS.PA.Jej@sam_data$Compartment=="Ascaris")
tmp1<- apply(X = t(otu_table(PS.tmp)),
            MARGIN = 1,
            FUN = function(x){sum(x > 0)})
tmp1<- data.frame(Prevalence_Ascaris = tmp1,
                  AscAbundance = taxa_sums(PS.tmp),
                 tax_table(PS.tmp))

tmp1%>%
  dplyr::select(c(Prevalence_Ascaris, AscAbundance))%>%
  cbind(Prev.JA)-> Prev.JA

rm(tmp, PS.tmp)

Prev.JA$TotalReads<- 1178219

Prev.JA%>%
  dplyr::mutate(Rel_abund= (TotalAbundance/TotalReads)*100)%>%
  dplyr::mutate(Rel_abund_Ascaris= (AscAbundance/TotalReads)*100)%>%
  dplyr::mutate(Rel_abund_Jej= (JejAbundance/TotalReads)*100)%>%
  dplyr::mutate(Rel_prev= (Prevalence/56)*100)%>%
  dplyr::mutate(Rel_prev_Ascaris= (Prevalence_Ascaris/47)*100)%>%
  dplyr::mutate(Rel_prev_Jej= (Prevalence_Jej/9)*100)-> Prev.JA

##Now differential abundance run
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

##Subset prevalence for those significant
Prev.JA<- Prev.JA[rownames(Prev.JA)%in%rownames(sigtab), ]

##Save this data
write.csv(sigtab, "Tables/Q2_DiffAbund_PA_JejvAsc.csv")
write.csv(Prev.JA, "Tables/Q2_Prev_PA_JejvAsc.csv")

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
require("ggrepel")
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(size=7, alpha= 0.5, position=position_jitter(0.2), aes(fill= AbundLev), shape= 21, color= "black")+
  scale_fill_manual(values=c("#D55E00","#E69F00", "#767676FF"), 
                    labels = c("High (Jejunum Infected)", "High (Ascaris)", "Not Significant"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(x= expression(Ascaris~parasite~~~~~~~~~~~~~~~~~~log[2]~Fold~change~~~~~~~~~~~~~~~~~~Pig~host),
       y= expression(-log[10]~(p~Adj)), fill= "Abundance level", tag = "A)")+
  theme_bw()+
  guides(fill = F)+
  geom_text_repel(data = subset(sigtab, AbundLev=="Low"),
                  aes(label = Bacteria_name),
                  size = 5,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),
                  max.overlaps = 15)+
  geom_text_repel(data = subset(sigtab, AbundLev=="High"),
                  aes(label = Bacteria_name),
                  size = 5,
                  box.padding = unit(0.3, "lines"),
                  point.padding = unit(0.3, "lines"),
                  max.overlaps = 12)+
  theme(text = element_text(size=16))-> B

##Get those significantly high or low abundant 
sigtab%>%
  dplyr::filter(AbundLev%in%c("High", "Low"))->tmp1

##Relative abundances 
otu <- as(otu_table(PS.PA.Jej), "matrix")
# transpose if necessary
if(taxa_are_rows(PS.PA.Jej)){otu <- t(otu)}
otu <- otu_table(otu, taxa_are_rows = F)
tax <- tax_table(PS.PA.Jej)
# Coerce to data.frame
n <- as.data.frame(tax)
n%>%
  rownames_to_column()%>%
  dplyr::rename(ASV = rowname)-> n
  
j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
j2 <- lapply(j1,'[[',"x") # select for Names
  
m <- data.frame(unlist(j2))
  
m%>%
  rownames_to_column()%>%
  dplyr::filter(unlist.j2.!=0)%>%
  dplyr::mutate(rowname = gsub(".ASV", "_ASV", rowname))%>%
  separate(rowname, sep = "_", c("Replicate","ASV"))%>%
  dplyr::group_by(Replicate)%>%
  dplyr::rename(Abundance = unlist.j2.)%>%
  dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
  left_join(n, by="ASV")%>%
  arrange(Replicate, desc(Genus))->m

m$Genus[is.na(m$Genus)]<- "Unassigned" ##Change NA's into Unassigned 
rm(otu, tax, j1, j2, n)

tmp<- PS.PA@sam_data

m%>%
  left_join(tmp, by= "Replicate")%>%
  dplyr::filter(ASV%in%rownames(tmp1))%>%
  dplyr::mutate(TaxaID= paste0(ASV, "-", Genus))%>%
  dplyr::mutate(Genus = fct_relevel(Genus, "Agathobacter", "Alloprevotella", "Anaerosporobacter", "Asaccharospora", 
                                    "Bifidobacterium", "Clostridium sensu stricto 1",
                                    "Coriobacteriaceae UCG-002", "Dialister",  "Escherichia-Shigella", "Lachnospira", 
                                    "Lactobacillus", "Megasphaera", "Peptococcus",
                                    "Prevotella", "Prevotellaceae UCG-001", "Prevotellaceae NK3B31 group", 
                                    "Pseudomonas", "Pseudoscardovia", "Roseburia", "Ruminococcus",
                                    "Staphylococcus", "Streptococcus"))%>%
  dplyr::mutate(OrderID = case_when(Genus =="Agathobacter"~1, Genus =="Alloprevotella"~2, Genus =="Anaerosporobacter"~3,
                                    Genus =="Asaccharospora"~4, Genus =="Bifidobacterium"~5, Genus =="Clostridium sensu stricto 1"~6,
                                    Genus =="Coriobacteriaceae UCG-002"~7, Genus =="Dialister"~8,  Genus =="Escherichia-Shigella"~9, 
                                    Genus =="Lachnospira"~10, Genus =="Lactobacillus"~11, Genus =="Megasphaera"~12, 
                                    Genus =="Peptococcus"~13, Genus =="Prevotella"~14, Genus =="Prevotellaceae UCG-001"~15, 
                                    Genus =="Prevotellaceae NK3B31 group"~15, Genus =="Pseudomonas"~16, Genus =="Pseudoscardovia"~17, 
                                    Genus =="Roseburia"~18, Genus =="Ruminococcus"~19,
                                    Genus =="Staphylococcus"~20, Genus =="Streptococcus"~21))%>%
  dplyr::mutate(TaxaID = fct_reorder(TaxaID, -OrderID))%>%
  ggplot(aes(x=Replicate, y=TaxaID)) + 
  geom_raster(aes(fill=Abundance)) +
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA, breaks=c(min("Abundance"),0.1, 0.4, max("Abundance")),
                      limits=c(0.00001,1))+
  theme_bw()+
  facet_grid(~Compartment, scales = "free_x")+
  labs(fill="Relative\nabundance (%)", tag = "B)")+
  theme(text = element_text(size=16), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank())
  
  
###Make a bubble plot for those ASV deferentially abundant
Ndom<- c("Agathobacter", "Alloprevotella", "Anaerosporobacter", "Asaccharospora", "Bifidobacterium",
  "Coriobacteriaceae UCG-002", "Dialister", "Lachnospira", "Megasphaera", "Peptococcus",
   "Prevotellaceae UCG-001", "Prevotellaceae NK3B31 group", "Pseudomonas", "Pseudoscardovia", "Roseburia",
  "Ruminococcus","Staphylococcus")

Prev.JA%>%
  rownames_to_column("ASV")%>%
  dplyr::select(c("ASV",  "Rel_prev_Ascaris", "Rel_abund_Ascaris", "Rel_prev_Jej", "Rel_abund_Jej", "Genus"))%>%
  dplyr::filter(ASV%in%rownames(tmp1))%>%
  dplyr::mutate(TaxaID= paste0(ASV, "-", Genus))%>%
  gather("Rel_abund_Ascaris", "Rel_abund_Jej", key = Host, value = Rel_abund, -ASV, -Genus, -Rel_prev_Ascaris, -Rel_prev_Jej)%>%
  gather("Rel_prev_Ascaris", "Rel_prev_Jej", key = Host2, value = Prevalence)%>%
  dplyr::mutate(Host= case_when(Host=="Rel_abund_Ascaris"~ "Ascaris",
                                Host=="Rel_abund_Jej"~"Pig-Jejunum"))%>%
  dplyr::mutate(Host2= case_when(Host2=="Rel_prev_Ascaris"~ "Ascaris",
                                Host2=="Rel_prev_Jej"~"Pig-Jejunum"))%>%
  dplyr::filter(Host == Host2)%>%
  dplyr::mutate(Genus = fct_relevel(Genus, "Agathobacter", "Alloprevotella", "Anaerosporobacter", "Asaccharospora", 
                                    "Bifidobacterium", "Clostridium sensu stricto 1",
                                    "Coriobacteriaceae UCG-002", "Dialister",  "Escherichia-Shigella", "Lachnospira", 
                                    "Lactobacillus", "Megasphaera", "Peptococcus",
                                    "Prevotella", "Prevotellaceae UCG-001", "Prevotellaceae NK3B31 group", 
                                    "Pseudomonas", "Pseudoscardovia", "Roseburia", "Ruminococcus",
                                    "Staphylococcus", "Streptococcus"))%>%
  dplyr::mutate(OrderID = case_when(Genus =="Agathobacter"~1, Genus =="Alloprevotella"~2, Genus =="Anaerosporobacter"~3,
                                    Genus =="Asaccharospora"~4, Genus =="Bifidobacterium"~5, Genus =="Clostridium sensu stricto 1"~6,
                                    Genus =="Coriobacteriaceae UCG-002"~7, Genus =="Dialister"~8,  Genus =="Escherichia-Shigella"~9, 
                                    Genus =="Lachnospira"~10, Genus =="Lactobacillus"~11, Genus =="Megasphaera"~12, 
                                    Genus =="Peptococcus"~13, Genus =="Prevotella"~14, Genus =="Prevotellaceae UCG-001"~15, 
                                    Genus =="Prevotellaceae NK3B31 group"~15, Genus =="Pseudomonas"~16, Genus =="Pseudoscardovia"~17, 
                                    Genus =="Roseburia"~18, Genus =="Ruminococcus"~19,
                                    Genus =="Staphylococcus"~20, Genus =="Streptococcus"~21))%>%
  dplyr::mutate(Genus = case_when(Genus%in%Ndom~ "Not dominant taxa",
                T~ as.character(Genus)))%>%
  dplyr::mutate(TaxaID = fct_reorder(TaxaID, -OrderID))%>%
  ggplot(aes(x = TaxaID, y = Prevalence))+
  geom_point(shape= 21, aes(size =as.factor(round(Rel_abund, digits = 2)), fill= Genus), color= "black", alpha= 0.75)+
  scale_fill_manual(values=c("Clostridium sensu stricto 1"= "#00468BFF","Escherichia-Shigella"= "#E762D7FF", 
                               "Lactobacillus"=  "#631879FF", "Prevotella" = "#E64B35FF", "Romboutsia" = "#F39B7FFF",
                               "Streptococcus"= "#925E9FFF", "Not dominant taxa"= "lightgray"))+
  coord_flip()+
  theme_bw()+
  facet_grid(~Host)+
  labs(size ="Relative abundance (%)", fill= "Genus",  y= "Prevalence (%)", tag = "B)")+
  guides(fill = guide_legend(override.aes=list(shape=c(21), size= 3)), size= guide_legend(nrow = 10), color= "none")+
  theme(text = element_text(size=16), axis.title.y = element_blank())-> C

##Ascaris SH Female vs Male
tmp<- row.names(PS.PA.Norm@sam_data)
tmp<- alphadiv.PA.rare[rownames(alphadiv.PA.rare)%in%tmp, ]

tmp%>%
  dplyr::filter(InfectionStatus!= "Non_infected")%>%
  dplyr::filter(Compartment%in% c("Ascaris"))%>%
  dplyr::group_by(System)%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.Asc.Norm<- subset_samples(PS.PA.Norm, Replicate%in%Inf.Keep)

DS.Asc.SH <- phyloseq_to_deseq2(PS.Asc.Norm, ~WormSex)

geoMeans <- apply(counts(DS.Asc.SH), 1, gm_mean)

DS.Asc.SH <- estimateSizeFactors(DS.Asc.SH, geoMeans = geoMeans)
DS.Asc.SH <- estimateDispersions(DS.Asc.SH, fitType= "mean")
DS.Asc.SH <- DESeq(DS.Asc.SH, test = "Wald", fitType= "mean")

res <- results(DS.Asc.SH, cooksCutoff = FALSE)
##Remove rows with columns that contain NA
res <- res[complete.cases(res), ]

sigtab <- cbind(as(res,"data.frame"), as(tax_table(PS.Asc.Norm)[rownames(res), ], "matrix"))
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
  scale_fill_manual(values=c("#00BDC4", "#F8766D", "#767676FF"), 
                    labels = c("High (Ascaris Male)", "High (Ascaris Female)", "Not Significant"))+
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

##Save plots for manuscript
ggsave(file = "Figures/Q1_Diff_Abundance_JejAsc.pdf", plot = B, width = 12, height = 12, dpi = 600)
ggsave(file = "Figures/Q1_Diff_Abundance_JejAsc.png", plot = B, width = 12, height = 12, dpi = 600)
ggsave(file = "Figures/Q1_Diff_Abundance_JejAsc.svg", plot = B, width = 12, height = 12, dpi = 600)
saveRDS(B, "Figures/Q1_Diff_Abundance_JejAsc.RDS") ##to compile with other

ggsave(file = "Figures/Q1_Diff_Abundance_AscFM.pdf", plot = D, width = 12, height = 8, dpi = 600)
ggsave(file = "Figures/Q1_Diff_Abundance_AscFM.png", plot = D, width = 12, height = 8, dpi = 600)
ggsave(file = "Figures/Q1_Diff_Abundance_AscFM.svg", plot = D, width = 12, height = 8, dpi = 600)
saveRDS(D, "Figures/Q1_Diff_Abundance_AscFM.RDS") ##to compile with other

###
Plot1<- cowplot::plot_grid(B, C, align = "vh", nrow = 2)
ggsave(file = "Figures/Figure_4_DiffAbund_adj.pdf", plot = Plot1, width = 12, height = 20)
ggsave(file = "Figures/Figure_4_DiffAbund_adj.png", plot = Plot1, width = 12, height = 20)
ggsave(file = "Figures/Figure_4_DiffAbund_adj.svg", plot = Plot1, width = 12, height = 20)

##Additional
RUN_also<- FALSE

if(RUN_also){
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
#write.csv(sigtab, "Tables/Q2_DiffAbund_Jej_InfvNonInf.csv")

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
#plot adding up all layers we have seen so far
require("ggrepel")
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(size=3, alpha= 0.5, position=position_jitter(0.2), aes(fill= AbundLev), shape= 21, color= "black")+
  scale_fill_manual(values=c("#009E73", "#D55E00", "#767676FF"), 
                    labels = c("High (Jejunum Non Infected)", "High (Jejunum Infected)", "Not Significant"))+
  geom_vline(xintercept=c(-0.6, 0.6), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "A)", x= "log2 Fold change", y= "-log10 (p Adjusted)", fill= "Abundance level")+
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
  theme(text = element_text(size=16))#-> A

#ggsave(file = "Figures/Q1_Diff_Abundance_JejInfNonInf.pdf", plot = A, width = 12, height = 8, dpi = 600)
#ggsave(file = "Figures/Q1_Diff_Abundance_JejInfNonInf.png", plot = A, width = 12, height = 8, dpi = 600)
#ggsave(file = "Figures/Q1_Diff_Abundance_JejInfNonInf.svg", plot = A, width = 12, height = 8, dpi = 600)
#saveRDS(A, "Figures/Q1_Diff_Abundance_JejInfNonInf.RDS") ##to compile with other

}