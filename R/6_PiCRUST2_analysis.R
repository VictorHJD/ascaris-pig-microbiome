##Project: Ascaris - Pig Microbiome
##Aim: Differential abundant taxa
##Author: Víctor Hugo Jarquín-Díaz
##Root repo setwd("../Ascaris/ascaris/")

###Analysis of predicted metagenomes generated in PICRUSt2
library(lme4)
library(lmtest)
library(ggplot2)
library(reshape2)
library(vegan)
library(gtools)
library(ggpubr)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(ALDEx2)
library(tidyverse)
library(viridis)
library("car")
library("merTools")
library(sjPlot)
library(sjlabelled)
library(sjmisc)

##Scaling abundances function
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

##Functions geometric mean 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


##1) Sample data
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
##1)General
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

##2)Predicted metagenomes
##Enzyme classification (EC) abundances per sample (Count table)
PredMet.PA<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_PA/EC_metagenome_out/pred_metagenome_unstrat.tsv", header = T, sep = "\t")
PredMet.Asc<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_Asc/EC_metagenome_out/pred_metagenome_unstrat.tsv", header = T, sep = "\t")

##3) Predicted KEGG onthology (KO)
##How the ASVs contribute to KEGG onthology (KO) abundances in each sample
PredKO.PA<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_PA/KO_metagenome_out/pred_metagenome_unstrat.tsv", header = T, sep = "\t")
PredKO.Asc<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_Asc/KO_metagenome_out/pred_metagenome_unstrat.tsv", header = T, sep = "\t")
##Descriptions
PredKO.PA.des<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_PA/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv", header = T, sep = "\t")
PredKO.Asc.des<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_Asc/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv", header = T, sep = "\t")

##4)Predicted pathway functions
##Metabolic pathway abundances per sample
PredPath.PA<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_PA/pathways_out/path_abun_unstrat.tsv", header = T, sep = "\t")
PredPath.Asc<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_Asc/pathways_out/path_abun_unstrat.tsv", header = T, sep = "\t")
##Descriptions
PredPath.PA.des<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_PA/pathways_out/path_abun_unstrat_descrip.tsv", header = T, sep = "\t")
PredPath.Asc.des<- read.table("/fast/AG_Forslund/Victor/data/Ascaris/PiCRUST/output_Asc/pathways_out/path_abun_unstrat_descrip.tsv", header = T, sep = "\t")

##Adjust tables
#Make functions/pathways the row names 
##PA
###Pathways
rownames(PredPath.PA)<-PredPath.PA$pathway
PredPath.PA$pathway<- NULL
PredPath.PA<- as.matrix(PredPath.PA)
PredPath.PA<- otu_table(PredPath.PA, taxa_are_rows = T)
sample_names(PredPath.PA)

###Descriptions
PredPath.PA.des%>% 
  dplyr::select(c(pathway,description))%>%
  dplyr::rename(Pathway= pathway)->PredPath.PA.des

###KO functions
rownames(PredKO.PA)<-PredKO.PA$function.
PredKO.PA$function.<- NULL
PredKO.PA<- as.matrix(PredKO.PA)
PredKO.PA<- otu_table(PredKO.PA, taxa_are_rows = T)
sample_names(PredKO.PA)
###Descriptions
PredKO.PA.des%>% 
  dplyr::select(c(function.,description))%>%
  dplyr::rename(KO= function.)->PredKO.PA.des

##Ascaris
###Pathway
rownames(PredPath.Asc)<-PredPath.Asc$pathway
PredPath.Asc$pathway<- NULL
PredPath.Asc<- as.matrix(PredPath.Asc)
PredPath.Asc<- otu_table(PredPath.Asc, taxa_are_rows = T)
sample_names(PredPath.Asc) 
###Descriptions
PredPath.Asc.des%>% 
  dplyr::select(c(pathway,description))%>%
  dplyr::rename(Pathway= pathway)->PredPath.Asc.des

###KO functions
rownames(PredKO.Asc)<-PredKO.Asc$function.
PredKO.Asc$function.<- NULL
PredKO.Asc<- as.matrix(PredKO.Asc)
PredKO.Asc<- otu_table(PredKO.Asc, taxa_are_rows = T)
sample_names(PredKO.Asc) 
###Descriptions
PredKO.Asc.des%>% 
  dplyr::select(c(function.,description))%>%
  dplyr::rename(KO= function.)->PredKO.Asc.des

#Prepare sample data 
##PA
###Sample information
alphadiv.PA%>%
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                                     "Experiment_1", "Experiment_2"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> alphadiv.PA
##Ascaris
###Sample information
alphadiv.Asc%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(Location = case_when(Origin %in% c("Experiment_1", "Experiment_2")  ~ "FU",
                                     Origin == "Slaughterhouse" ~ "SH"))%>%
  dplyr::mutate(Location = fct_relevel(Location, "FU", "SH"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3","Pig4",
                                     "Pig5","Pig6","Pig7","Pig8","Pig9",
                                     "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> alphadiv.Asc

#Make Phyloseq objects for pathways and KO
##PA
###Pathways
PS.Path.PA <- merge_phyloseq(PredPath.PA, sample_data(alphadiv.PA))
###KO
PS.KO.PA <- merge_phyloseq(PredKO.PA, sample_data(alphadiv.PA))

##Ascaris
###Pathways
PS.Path.Asc <- merge_phyloseq(PredPath.Asc, sample_data(alphadiv.PA))
###KO
PS.KO.Asc <- merge_phyloseq(PredKO.Asc, sample_data(alphadiv.Asc))

###Let's start the analysis 
##################Generate PCAs#####################
##KO
##Transform dataset to determine contributors
PS.KO.PA.clr <- microbiome::transform(PS.KO.PA, "clr") #Centered log ratio transformation
Ord.KO.PA.clr <- phyloseq::ordinate(PS.KO.PA.clr, "RDA") #principal components analysis

#Examine eigenvalues and % prop. variance explained
head(Ord.KO.PA.clr$CA$eig)
sapply(Ord.KO.PA.clr$CA$eig[1:6], function(x) x / sum(Ord.KO.PA.clr$CA$eig))

##KOs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.KO.PA.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))

ind_cont_PCA1 %>% 
  rownames_to_column("KO") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 KO contribute for the 4.09% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("KO") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 KO contribute for the 2.13% of the variation in PC2

ind_cont_PCA_top.PA <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge description
ind_cont_PCA_top.PA%>%
  left_join(PredKO.PA.des, by= "KO")-> ind_cont_PCA_top.PA

##Estimate alpha diversity of KO Richness
alpha.picrusts.PA<- microbiome::alpha(PS.KO.PA, c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))

alphadiv.PA%>%
  select(10:16,)%>%
  cbind(alpha.picrusts.PA)-> alpha.picrusts.PA

### Infected vs Non Infected (Richness)
alpha.picrusts.PA%>% 
  dplyr::filter(InfectionStatus!= "Worm")%>%
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                          "Duodenum", "Jejunum", "Ileum", 
                                          "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(chao1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
  write.csv(x, "Tables/Q1_Alpha_PICRUST_Infection.csv")

##Plot 
alpha.picrusts.PA%>% 
    dplyr::filter(InfectionStatus!= "Worm")%>%
    dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                            "Duodenum", "Jejunum", "Ileum", 
                                            "Cecum", "Colon"))%>%
    dplyr::group_by(Compartment)%>%
    mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= chao1))+
  geom_boxplot(aes(color= InfectionStatus, fill= InfectionStatus), outlier.shape=NA)+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  xlab("GI compartment")+
  ylab("Predicted KO Richness (Chao1 Index)")+
  labs(tag= "A)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -2, step.increase = 0.05, hide.ns = T,
                     tip.length = 0)-> A

alpha.picrusts.PA%>% 
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                          "Duodenum", "Jejunum", "Ileum", 
                                          "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(InfectionStatus== "Infected")%>%
  dplyr::group_by(System)%>%
  dplyr::filter(n()==5)%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

##Plot 
alpha.picrusts.PA%>%
  dplyr::filter(InfectionStatus== "Infected")%>%
  dplyr::filter(Replicate%in%Inf.Keep)%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_line(aes(group = System), colour= "gray")+
  geom_point(shape= 21, size=3, aes(fill= System), color= "black")+
  scale_fill_manual(values = pal.system)+
  xlab("GI compartment")+
  ylab("Predicted KO Richness (Chao1 Index)")+
  labs(tag= "B)", 
       shape = "Infection status", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+ #, caption = get_pwc_label(stats.test)
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())-> B

##Comparison Alpha diversity between worms and site of infection
alpha.picrusts.PA%>%
  dplyr::filter(Compartment%in% c("Jejunum","Ascaris"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  wilcox_test(chao1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_PICRUST_Jej_Ascaris.csv")

##Plot 
alpha.picrusts.PA%>%
  dplyr::filter(Compartment%in% c("Jejunum", "Ascaris"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= InfectionStatus, y= chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, aes(fill= System, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 25, 21), labels = c("Infected", " Non infected", "Worm"))+
  scale_fill_manual(values = pal.system)+
  xlab("Infection status")+
  ylab("Predicted KO Richness (Chao1 Index)")+
  labs(tag= "C)", caption = get_pwc_label(stats.test), 
       shape = "Infection status", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -0.2, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> C

plot1<-plot_grid(A, B, ncol=1, align="v")

D<- grid.arrange(plot1, C, widths = c(4, 2.5),
                 layout_matrix = rbind(c(1, 2)))

ggsave(file = "Figures/Q2_Alpha_PICRUST_Compartment.pdf", plot = D, width = 12, height = 8, dpi = 600)
ggsave(file = "Figures/Q2_Alpha_PICRUST_Compartment.png", plot = D, width = 12, height = 8, dpi = 600)

rm(A,B,C,D, plot1)

## Bray Curtis
##Just infected and not infected 
PS.KO.pig<- subset_samples(PS.KO.PA, InfectionStatus!= "Worm")
##Normalize to RA
PS.KO.pig.ra <- microbiome::transform(PS.KO.pig, "compositional")

BC_dist<- phyloseq::distance(PS.KO.pig.ra, method="bray", weighted=F)
Ord.KO.PA <- ordinate(PS.KO.pig.ra, method="PCoA", distance="bray") 

plot_ordination(PS.KO.pig.ra, ordination = Ord.KO.PA)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= InfectionStatus, shape= InfectionStatus), color= "black")+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  labs(tag= "A)", fill  = "Infection status", color= "Origin of samples")+
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  labs(title = "Bray-Curtis dissimilariy (KO prediction)",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F, color= F)+
  xlab(paste0("PCo 1 [", round(Ord.KO.PA$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(Ord.KO.PA$values[2,2]*100, digits = 2), "%]")) -> A

##Normalize to RA
##Just infected and Ascaris
PS.KO.PA.ra<- subset_samples(PS.KO.PA, InfectionStatus%in%c("Infected", "Worm"))

PS.KO.PA.ra <- microbiome::transform(PS.KO.PA.ra, "compositional")

BC_dist<- phyloseq::distance(PS.KO.PA.ra, method="bray", weighted=F)
Ord.KO.PA <- ordinate(PS.KO.PA.ra, method="PCoA", distance="bray") 

plot_ordination(PS.KO.PA.ra, ordination = Ord.KO.PA)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Compartment, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 21), labels = c("Infected", "Ascaris"))+
  scale_fill_manual(values = pal.compartment)+
  stat_ellipse(aes(color= Compartment), linetype = 2)+
  scale_color_manual(values = pal.compartment)+
  labs(title = "Bray-Curtis dissimilariy (KO prediction)",tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  labs(fill = "Compartment")+
  labs(shape = "Infection Status")+
  xlab(paste0("PCo 1 [", round(Ord.KO.PA$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(Ord.KO.PA$values[2,2]*100, digits = 2), "%]")) -> B

tmp<- row.names(PS.KO.PA.ra@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

BC.KO.PA<- vegan::adonis(BC_dist~ InfectionStatus + Origin,
                         permutations = 999, data = tmp, na.action = T, strata = tmp$System)
##Store data
foo<- as.data.frame(BC.KO.PA$aov.tab)
#write.csv(foo, file = "Tables/Q2_Adonis_PICRUST_Compartments_Ascaris.csv")

####Subset pigs with before, after, site of infection and worm data
tmp%>%
  dplyr::filter(Compartment%in% c("Jejunum", "Duodenum", "Ileum", "Ascaris"))%>%
  dplyr::filter(System%in% c("Pig1", "Pig3", "Pig5", "Pig10", "Pig12", "Pig13"))%>%
  dplyr::select(Replicate)-> Inf.Keep
  
Inf.Keep<- Inf.Keep$Replicate

PS.KO.Jej<- subset_samples(PS.KO.PA, Replicate%in%Inf.Keep)

tmp<- row.names(PS.KO.Jej@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

####
PS.KO.Jej.ra <- microbiome::transform(PS.KO.Jej, "compositional")
bray_dist<- phyloseq::distance(PS.KO.Jej.ra, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.KO.Jej.ra,
                      method="PCoA", distance="bray")

JejAsc.adonis<- vegan::adonis(bray_dist~ InfectionStatus + Origin,
                              permutations = 999, data = tmp, na.action = T, strata = tmp$System)
##Store data
foo<- as.data.frame(JejAsc.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q2_Adonis_PICRUST_Jejunum_Ascaris.csv")

## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, tmp$Compartment, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2 
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("Compartment","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

##Add sample data
tmp%>%
  dplyr::select(!c(Compartment))%>%
  cbind(seg.data)-> seg.data

##PCoA
ggplot() + 
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2, color= grps, group=grps), size=4, shape= 4) +
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Compartment, shape= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24,21), labels= c("Infected Pig",  "Ascaris"))+
  scale_fill_manual(values = pal.compartment)+
  scale_color_manual(values = pal.compartment)+
  labs(tag= "C)", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Compartment), linetype = 2)+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> C

##Just Site of infection and ascaris
tmp<- row.names(PS.KO.PA@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  dplyr::filter(Compartment%in% c("Jejunum", "Ascaris"))%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.KO.Jej<- subset_samples(PS.KO.PA, Replicate%in%Inf.Keep)

tmp<- row.names(PS.KO.Jej@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

####
PS.KO.Jej.ra <- microbiome::transform(PS.KO.Jej, "compositional")
bray_dist<- phyloseq::distance(PS.KO.Jej.ra, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.KO.Jej.ra,
                      method="PCoA", distance="bray")

JejAsc.adonis<- vegan::adonis(bray_dist~ InfectionStatus + Origin,
                              permutations = 999, data = tmp, na.action = T, strata = tmp$System)
##Store data
#foo<- as.data.frame(JejAsc.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q2_Adonis_PICRUST_Jejunum_Ascaris.csv")

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
  labs(tag= "D)", fill  = "Infection status", color= "Origin of samples")+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> D

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
write.csv(x, "Tables/Q1_Beta_PICRUST_PCo1_JejAscInf.csv")

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
write.csv(x, "Tables/Q1_Beta_PICRUST_PCo2_JejAscInf.csv")

stats.test%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(y.position= c(0.25))-> stats.test

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
p1<- insert_xaxis_grob(D, xplot, grid::unit(.25, "null"), position = "top")
p2<- insert_yaxis_grob(p1, yplot, grid::unit(.25, "null"), position = "right")
D2<- ggdraw(p2)

E<- grid.arrange(A,B, C, D2)

ggsave(file = "Figures/Q2_Beta_PICRUST_Compartment.pdf", plot = E, width = 12, height = 8, dpi = 600)
ggsave(file = "Figures/Q2_Beta_PICRUST_Compartment.png", plot = E, width = 12, height = 8, dpi = 600)

rm(A,B,C,D2, E)

##Pathways
PS.Path.PA.clr <- microbiome::transform(PS.Path.PA, "clr")  
Ord.Path.PA.clr <- phyloseq::ordinate(PS.Path.PA.clr, "RDA")
#Examine eigenvalues and % prop. variance explained
head(Ord.Path.PA.clr$CA$eig)
sapply(Ord.Path.PA.clr$CA$eig[1:6], function(x) x / sum(Ord.Path.PA.clr$CA$eig))

##KOs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.Path.PA.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))

ind_cont_PCA1 %>% 
  rownames_to_column("Pathway") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 KO contribute for the 42.6% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("Pathway") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 KO contribute for the 28.3% of the variation in PC2

ind_cont_PCA_top.ptw.PA <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge description
ind_cont_PCA_top.ptw.PA%>%
  left_join(PredPath.PA.des, by= "Pathway")-> ind_cont_PCA_top.ptw.PA

##Plot
plot_ordination(PS.Path.PA.clr, ordination = Ord.Path.PA.clr)+
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Compartment, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 25, 21))+
  scale_fill_manual(values = pal.compartment)+
  labs(title = "PCA (centered-log ratio pathway prediction)",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Compartment")+
  labs(shape = "Infection Status")+
  xlab(paste0("PC 1 [", round(Ord.Path.PA.clr$CA$eig[1] / sum(Ord.Path.PA.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.Path.PA.clr$CA$eig[2] / sum(Ord.Path.PA.clr$CA$eig)*100, digits = 2), "%]"))

## Bray Curtis
PS.Path.PA.ra <- microbiome::transform(PS.Path.PA, "compositional")  
BC_dist<- phyloseq::distance(PS.Path.PA.ra, method="bray", weighted=F)
Ord.path.PA <- ordinate(PS.Path.PA.ra, method="PCoA", distance="bray") 

plot_ordination(PS.Path.PA.ra, ordination = Ord.path.PA)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(fill= Compartment, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 25, 21), labels = c("Infected", "Non infected", "Ascaris"))+
  scale_fill_manual(values = pal.compartment)+
  labs(title = "Bray-Curtis dissimilarity (pathway prediction)",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  labs(fill = "Compartment")+
  labs(shape = "Infection Status")+
  xlab(paste0("PCo 1 [", round(Ord.path.PA$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(Ord.path.PA$values[2,2]*100, digits = 2), "%]")) -> A

tmp<- row.names(PS.Path.PA.ra@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

BC.path.PA<- vegan::adonis(BC_dist~ InfectionStatus + Compartment + Origin,
                               permutations = 999, data = tmp, na.action = F, strata = tmp$System)

##Store data
foo<-BC.path.PA$aov.tab#-> Report 
#write.csv(foo, file = "Tables/Q2_Adonis_PICRUST_Path_Compartment.csv")

ggsave(file = "Figures/Q2_path_pred_ord_PA.pdf", plot = A, width = 10, height = 10)
ggsave(file = "Figures/Q2_path_pred_ord_PA.png", plot = A, width = 10, height = 10)

####Heat map pathways#####
##Matrix Pathways sample
PathM<- PS.Path.PA.ra@otu_table

top.ptw<- ind_cont_PCA_top.ptw.PA[complete.cases(ind_cont_PCA_top.ptw.PA), ]

top.ptw<- top.ptw$Pathway

PathM<- PathM[rownames(PathM) %in% top.ptw]

P.clust <- hclust(dist(t(PathM)), method = "complete") ##Dendogram

as.dendrogram(P.clust) %>%
  plot(horiz = T)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))

P.col<- cbind(P.col, alphadiv.PA)

col_groups <- P.col %>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2"))%>%
  dplyr::select(c(System, Compartment, Origin, InfectionStatus, Replicate)) ##Here It is possible to add the other characteristics

col_groups$Replicate<- NULL

colour_groups <- list(System= pal.system, Compartment= pal.compartment, 
                      InfectionStatus=  c("Infected"="#ED0000FF", "Non_infected"= "#008B45FF", "Worm" ="#fdae61"))

require(pheatmap)
pheatmap(PathM, cluster_rows = F, cluster_cols = T,
                           color = colorRampPalette(c("#67001f","#f7f7f7","#053061"))(100), 
                           border_color = NA,
                           annotation_col = col_groups, 
                           annotation_colors = colour_groups,
                           show_rownames = T,
                           show_colnames = F)

###Differentially abundant functions
##Infected and non-infected pigs

##Merge Genus and species 
sigtab%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("UCG-005", "Oscillospiraceae UCG-005", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> sigtab

##Save this data
write.csv(sigtab, "Tables/Q2_DiffAbund_Jej_InfvNonInf.csv")

##Analysis for Infected vs Non-infected pigs in site of infection
DS.Func.Inf <- phyloseq_to_deseq2(PS.KO.Jej, ~InfectionStatus)

geoMeans <- apply(counts(DS.Func.Inf), 1, gm_mean)

DS.Func.Inf<- estimateSizeFactors(DS.Func.Inf, geoMeans = geoMeans)
DS.Func.Inf<- DESeq(DS.Func.Inf)

res <- results(DS.Func.Inf, cooksCutoff = FALSE)

PredKO.PA.des%>%
  column_to_rownames("KO")-> tmp.ko

sigtab <- as.data.frame(res)

##Remove rows with columns that contain NA
sigtab <- sigtab[complete.cases(sigtab$padj), ]
head(sigtab,25)

##Volcano plot to detect differential genes in Non infected vs Infected

ggplot(sigtab, aes(x=log2FoldChange, y= -log10(padj))) +
  geom_point(size=1) +
  theme_bw()+
  geom_vline(xintercept = c(-1, 1), col= "red", linetype= "dashed")+
  geom_hline(yintercept = -log10(0.001), col="red", linetype= "dashed")

#The significantly deferentially abundant genes are the ones found upper-left and upper-right corners
##Add a column to the data specifying if they are highly (positive) or lowly abundant (negative)
## Considering the comparison Wildling vs Black 6

# add a column of Non-significant
sigtab$AbundLev <- "NS"
# if log2Foldchange > 1 and pvalue < 0.01, set as "High" 
sigtab$AbundLev[sigtab$log2FoldChange > 1 & sigtab$padj < 0.001] <- "High"
# if log2Foldchange < -1 and pvalue < 0.01, set as "DOWN"
sigtab$AbundLev[sigtab$log2FoldChange < -1 & sigtab$padj < 0.001] <- "Low"

#Organize the labels nicely using the "ggrepel" package and the geom_text_repel() function

#plot adding up all layers we have seen so far
sigtab%>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(size=3, alpha= 0.5, aes(fill= AbundLev), shape= 21, color= "black")+
  scale_fill_manual(values = c("#008B45FF", "#ED0000FF", "gray"), labels = c("High (Non Infected)", "High (Infected)", "Not Significant"))+
  geom_vline(xintercept=c(-1, 1), col="black", linetype= "dashed") +
  geom_hline(yintercept=-log10(0.001), col="black", linetype= "dashed") +
  labs(tag= "A)", x= "log2 Fold change", y= "-Log10 (p Adjusted)", fill= "Gene \n abundance")+
  theme_bw()+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  geom_text_repel(data = subset(sigtab, AbundLev=="Low"),
                  aes(label = Module_C_name),
                  size = 3,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"))+
  geom_text_repel(data = subset(sigtab, AbundLev=="High"),
                  aes(label = Module_C_name),
                  size = 3,
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"))+
  theme(text = element_text(size=12))-> A

#ggsave(file = "results/figures/Q4_Gene_Abundance_WildB6_bin.pdf", plot = A, width = 10, height = 8)
#ggsave(file = "results/figures/Q4_Gene_Abundance_WildB6_bin.png", plot = A, width = 10, height = 8)

##Extract Highly and lowly abundant genes between Wildlings from d21 and Black 6 from both collections
sigtab%>%
  dplyr::select(c(baseMean, log2FoldChange, padj, Module_A_name, 
                  Module_B_name, Module_C_name, AbundLev))%>%
  dplyr::filter(AbundLev!= "NS")%>%
  rownames_to_column(var = "GenID")%>%
  left_join(gene.names, by= "GenID")%>%
  dplyr::select(c(GenID, log2FoldChange, padj, AbundLev, Module_A_name, Module_B_name, Module_C_name, Kegg_Accession))%>%
  left_join(KO.ref, by= "Kegg_Accession")%>%
  distinct(GenID, .keep_all = TRUE)-> Genes.wild

write.csv(Genes.wild, file = "results/Abundant_Genes_WildB6_bin.csv", row.names=FALSE)

##Just Site of infection and ascaris
tmp<- row.names(PS.Path.PA@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  dplyr::filter(Compartment%in% c("Jejunum", "Ascaris"))%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.Path.Jej<- subset_samples(PS.Path.PA, Replicate%in%Inf.Keep)

tmp<- row.names(PS.Path.Jej@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

####
PS.path.Jej.ra <- microbiome::transform(PS.Path.Jej, "compositional")
bray_dist<- phyloseq::distance(PS.path.Jej.ra, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.path.Jej.ra,
                      method="PCoA", distance="bray")

JejAsc.path.adonis<- vegan::adonis(bray_dist~ InfectionStatus + Origin,
                              permutations = 999, data = tmp, na.action = T, strata = tmp$System)

##Pathways
PS.Path.Jej.clr <- microbiome::transform(PS.Path.Jej, "clr")  
Ord.Path.Jej.clr <- phyloseq::ordinate(PS.Path.Jej.clr, "RDA")
#Examine eigenvalues and % prop. variance explained
head(Ord.Path.Jej.clr$CA$eig)
sapply(Ord.Path.Jej.clr$CA$eig[1:6], function(x) x / sum(Ord.Path.Jej.clr$CA$eig))

##KOs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.Path.Jej.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))

ind_cont_PCA1 %>% 
  rownames_to_column("Pathway") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 KO contribute for the 42.6% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("Pathway") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 KO contribute for the 28.3% of the variation in PC2

ind_cont_PCA_top.ptw.Jej <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge description
ind_cont_PCA_top.ptw.Jej%>%
  left_join(PredPath.PA.des, by= "Pathway")-> ind_cont_PCA_top.ptw.Jej

##Matrix Pathways sample
PathM<- PS.path.Jej.ra@otu_table

top.ptw<- ind_cont_PCA_top.ptw.Jej[complete.cases(ind_cont_PCA_top.ptw.Jej), ]

top.ptw<- top.ptw$Pathway

PathM<- PathM[rownames(PathM) %in% top.ptw]

P.clust <- hclust(dist(t(PathM)), method = "complete") ##Dendogram

as.dendrogram(P.clust) %>%
  plot(horiz = T)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))

P.col<- cbind(P.col, alphadiv.PA)

col_groups <- P.col %>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2"))%>%
  dplyr::select(c(System, Compartment, Origin, InfectionStatus, Replicate)) ##Here It is possible to add the other characteristics

col_groups$Replicate<- NULL

colour_groups <- list(System= pal.system, Compartment= pal.compartment, 
                      InfectionStatus=  c("Infected"="#ED0000FF", "Non_infected"= "#008B45FF", "Worm" ="#fdae61"))

require(pheatmap)
pheatmap(PathM, cluster_rows = F, cluster_cols = T,
         color = colorRampPalette(c("#67001f","#f7f7f7","#053061"))(100), 
         border_color = NA,
         annotation_col = col_groups, 
         annotation_colors = colour_groups,
         show_rownames = T,
         show_colnames = F)

##Just Site of infection and ascaris
tmp<- row.names(PS.Path.PA@sam_data)
tmp<- alphadiv.pig[rownames(alphadiv.pig)%in%tmp, ]

tmp%>%
  dplyr::filter(Compartment%in% c("Jejunum"))%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.Path.Jej<- subset_samples(PS.Path.PA, Replicate%in%Inf.Keep)

tmp<- row.names(PS.Path.Jej@sam_data)
tmp<-  alphadiv.pig[rownames( alphadiv.pig)%in%tmp, ]

tmp%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

####
PS.path.Jej.ra <- microbiome::transform(PS.Path.Jej, "compositional")
bray_dist<- phyloseq::distance(PS.path.Jej.ra, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.path.Jej.ra,
                      method="PCoA", distance="bray")

JejAsc.path.adonis<- vegan::adonis(bray_dist~ InfectionStatus + Origin,
                                   permutations = 999, data = tmp, na.action = T, strata = tmp$System)

##Pathways
PS.Path.Jej.clr <- microbiome::transform(PS.Path.Jej, "clr")  
Ord.Path.Jej.clr <- phyloseq::ordinate(PS.Path.Jej.clr, "RDA")
#Examine eigenvalues and % prop. variance explained
head(Ord.Path.Jej.clr$CA$eig)
sapply(Ord.Path.Jej.clr$CA$eig[1:6], function(x) x / sum(Ord.Path.Jej.clr$CA$eig))

##KOs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.Path.Jej.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))

ind_cont_PCA1 %>% 
  rownames_to_column("Pathway") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 KO contribute for the 42.6% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("Pathway") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 KO contribute for the 28.3% of the variation in PC2

ind_cont_PCA_top.ptw.Jej <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge description
ind_cont_PCA_top.ptw.Jej%>%
  left_join(PredPath.PA.des, by= "Pathway")-> ind_cont_PCA_top.ptw.Jej

##Matrix Pathways sample
PathM<- PS.path.Jej.ra@otu_table

top.ptw<- ind_cont_PCA_top.ptw.Jej[complete.cases(ind_cont_PCA_top.ptw.Jej), ]

top.ptw<- top.ptw$Pathway

PathM<- PathM[rownames(PathM) %in% top.ptw]

P.clust <- hclust(dist(t(PathM)), method = "complete") ##Dendogram

as.dendrogram(P.clust) %>%
  plot(horiz = T)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))

P.col<- cbind(P.col, alphadiv.pig)

col_groups <- P.col %>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2"))%>%
  dplyr::select(c(System, Compartment, Origin, InfectionStatus, Replicate)) ##Here It is possible to add the other characteristics

col_groups$Replicate<- NULL

colour_groups <- list(System= pal.system, Compartment= pal.compartment, 
                      InfectionStatus=  c("Infected"="#ED0000FF", "Non_infected"= "#008B45FF", "Worm" ="#fdae61"))

require(pheatmap)
pheatmap(PathM, cluster_rows = F, cluster_cols = T,
         color = colorRampPalette(c("#67001f","#f7f7f7","#053061"))(100), 
         border_color = NA,
         #annotation_col = col_groups, 
         #annotation_colors = colour_groups,
         show_rownames = T,
         show_colnames = F)
