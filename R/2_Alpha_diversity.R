##Project: Ascaris - Pig Microbiome
##Aim: Alpha diversity
##Author: Víctor Hugo Jarquín-Díaz
##Root repo setwd("../Ascaris/ascaris/")

##Load libraries 
library(phyloseq)
library(microbiome)
library(tidyverse)
require(ggpubr)
require(RColorBrewer)
require(rstatix)
library(cowplot)
library(gridExtra)
library(grid)
library(ggsci)
library(microbiome)

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

##Color palette for compartment and sytem ##

pal.compartment <- c("Ascaris"="#1B9E77","Cecum"= "#D95F02","Colon"= "#7570B3",
                     "Duodenum"= "#E7298A","Ileum"= "#66A61E","Jejunum"="#E6AB02")

pal.system <- c("Pig1"= "#A6761D","Pig2"= "#666666","Pig3"= "#A6CEE3","Pig4"= "#1F78B4",
           "Pig5"= "#B2DF8A","Pig6"= "#33A02C","Pig7"= "#FB9A99","Pig8"="#E31A1C","Pig9"= "#FDBF6F",
           "Pig10"= "#FF7F00","Pig11"= "#CAB2D6","Pig12"= "#6A3D9A","Pig13"= "#FFFF99",  "Pig14"= "#3B3B3BFF", "SH" = "#767676FF")

##Functions 
##Find dominant taxa per samples
find.top.asv <- function(x, taxa, num){
  require(phyloseq)
  require(magrittr)
  
  top.taxa <- tax_glom(x,taxa)
  otu <- as(otu_table(top.taxa), "matrix")
  # transpose if necessary
  if(taxa_are_rows(top.taxa)){otu <- t(otu)}
  otu <- otu_table(otu, taxa_are_rows = F)
  tax <- tax_table(top.taxa)
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
    separate(rowname, c("SampleID", "ASV"))%>%
    dplyr::group_by(SampleID)%>%
    slice_max(order_by = unlist.j2., n = num)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")->m
  
  rm(top.taxa, otu, tax, j1, j2, n)
  return(m)
}

##Get data frame for bar plot at genus level 
count.high.genus <- function(x, num){
  require(phyloseq)
  require(magrittr)
  #x is a phyloseq object glomed to Genus
  #num is the threshold of Relative abundance desired 
  otu <- as(otu_table(x), "matrix")
  # transpose if necessary
  if(taxa_are_rows(x)){otu <- t(otu)}
  otu <- otu_table(otu, taxa_are_rows = F)
  tax <- tax_table(x)
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
    separate(rowname, c("SampleID", "ASV"))%>%
    dplyr::group_by(SampleID)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")%>%
    mutate(Main_taxa= Abundance>= num)%>%
    dplyr::mutate(Genus= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(Genus)))%>%
    arrange(SampleID, desc(Genus))->m
  
  m$Genus[is.na(m$Genus)]<- "Unassigned" ##Change NA's into Unassigned 
  m$Species<- NULL
  
  rm(otu, tax, j1, j2, n)
  return(m)
}


##Transform abundance into relative abundance
Rel.abund_fun <- function(df){
  df2 <- sapply(df, function(x) (x/1E6)*100)  
  colnames(df2) <- colnames(df)
  rownames(df2) <- rownames(df)
  df2<- as.data.frame(df2)
  return(df2)
}

#########Question 1:
###How does Ascaris impact the porcine microbiome and does this differ in different gut regions?
###General comparison between infected and non infected pigs (all compartments, merged replicates) ###################
###Group comparisons

### Infected vs Non Infected
alphadiv.pig%>% 
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(Chao1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Infected_NonInfected_Compartment.csv")

alphadiv.pig%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_effsize(Chao1 ~ InfectionStatus)

##Plot 
alphadiv.pig%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                                   "Pig1","Pig2","Pig3","Pig4",
                                   "Pig5","Pig6","Pig7","Pig8","Pig9",
                                   "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= Chao1))+
  geom_boxplot(aes(color= InfectionStatus), alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, shape= 21,aes(fill= InfectionStatus), color= "black")+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  xlab("GI compartment")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test), 
       fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -2, step.increase = 0.05, hide.ns = T,
                     tip.length = 0)-> A

##Different among compartment
#Just Infected with points for all compartments

alphadiv.pig%>%
  dplyr::filter(InfectionStatus== "Infected")%>%
  dplyr::group_by(System)%>%
  dplyr::filter(n()==5)%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

alphadiv.pig%>%
  dplyr::filter(InfectionStatus== "Infected")%>%
  dplyr::filter(Replicate%in%Inf.Keep)%>%
  wilcox_test(Chao1 ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Infected_Compartment.csv")

alphadiv.pig%>% 
  dplyr::filter(InfectionStatus== "Infected")%>%
  dplyr::filter(Replicate%in%Inf.Keep)%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  wilcox_effsize(Chao1 ~ Compartment)

##Plot 
alphadiv.pig%>%
  dplyr::filter(InfectionStatus== "Infected")%>%
  dplyr::filter(Replicate%in%Inf.Keep)%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(shape= 21, size=3, aes(fill= System), color= "black")+
  scale_fill_manual(values = pal.system)+
  geom_line(aes(group = System), colour= "gray")+
  xlab("GI compartment")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test), 
       shape = "Infection status", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -0.2, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> B

C<-grid.arrange(A,B)

ggsave(file = "Figures/Q1_Alpha_Compartment.pdf", plot = C, width = 10, height = 8)
ggsave(file = "Figures/Q1_Alpha_Compartment.png", plot = C, width = 10, height = 8)


