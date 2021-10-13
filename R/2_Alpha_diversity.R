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

##Color palette for compartment and system ##

pal.compartment <- c("Ascaris"="#1B9E77","Cecum"= "#D95F02","Colon"= "#7570B3",
                     "Duodenum"= "#E7298A","Ileum"= "#66A61E","Jejunum"="#E6AB02")

pal.system <- c("Pig1"= "#A6761D","Pig2"= "#666666","Pig3"= "#A6CEE3","Pig4"= "#1F78B4",
           "Pig5"= "#B2DF8A","Pig6"= "#33A02C","Pig7"= "#FB9A99","Pig8"="#E31A1C","Pig9"= "#FDBF6F",
           "Pig10"= "#FF7F00","Pig11"= "#CAB2D6","Pig12"= "#6A3D9A","Pig13"= "#FFFF99",  "Pig14"= "#3B3B3BFF", "SH" = "#BB0021FF")

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
    dplyr::mutate(rowname = gsub(".ASV", "_ASV", rowname))%>%
    separate(rowname, sep = "_", c("Replicate","ASV"))%>%
    dplyr::group_by(Replicate)%>%
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
    dplyr::mutate(rowname = gsub(".ASV", "_ASV", rowname))%>%
    separate(rowname, sep = "_", c("Replicate","ASV"))%>%
    dplyr::group_by(Replicate)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")%>%
    mutate(Main_taxa= Abundance>= num)%>%
    dplyr::mutate(Genus= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(Genus)))%>%
    arrange(Replicate, desc(Genus))->m
  
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

### Infected vs Non Infected (Richness)
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
  geom_boxplot(aes(color= InfectionStatus, fill= InfectionStatus),outlier.shape=NA)+
  #geom_point(pch = 21, position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  xlab("GI compartment")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "A)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -2, step.increase = 0.05, hide.ns = T,
                     tip.length = 0)+
  scale_y_continuous(limits=c(0, 1000))-> A

##Phylogenetic richness
require("btools")
foo<- estimate_pd(PS.pig)

alphadiv.pig<- cbind(alphadiv.pig, foo)

### Infected vs Non Infected (Diversity)
alphadiv.pig%>% 
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                          "Duodenum", "Jejunum", "Ileum", 
                                          "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(Shannon ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Infected_NonInfected_Compartment_Shannon.csv")

##Plot 
alphadiv.pig%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= Shannon))+
  geom_boxplot(aes(color= InfectionStatus, fill= InfectionStatus), outlier.shape=NA)+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  xlab("GI compartment")+
  ylab("ASV Diversity (Shannon Index)")+
  labs(tag= "A)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, step.increase = 0.05, hide.ns = T,
                     tip.length = 0)+
  scale_y_continuous(limits=c(0, 5))-> Sup1A

##Phylogenetic diversity
alphadiv.pig%>% 
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                          "Duodenum", "Jejunum", "Ileum", 
                                          "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(PD ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Infected_NonInfected_Compartment_PD.csv")

##Plot 
alphadiv.pig%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= PD))+
  geom_boxplot(aes(color= InfectionStatus, fill= InfectionStatus), outlier.shape=NA)+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  xlab("GI compartment")+
  ylab("Phylogenetic diverstiy (Faith's Index)")+
  labs(tag= "B)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, step.increase = 0.05, hide.ns = T,
                     tip.length = 0)-> Sup1B

##Difference among compartment
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
  add_xy_position(x = "Compartment", dodge = 0.8)-> stats.test 

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
  geom_line(aes(group = System), colour= "gray")+
  geom_point(shape= 21, size=3, aes(fill= System), color= "black")+
  scale_fill_manual(values = pal.system)+
  xlab("GI compartment")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "B)", 
       shape = "Infection status", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+ #, caption = get_pwc_label(stats.test)
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -1000, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> B

plot1<-plot_grid(A, B, ncol=1, align="v")

#ggsave(file = "Figures/Q1_Alpha_Compartment.pdf", plot = plot1, width = 10, height = 8, dpi = 400)
#ggsave(file = "Figures/Q1_Alpha_Compartment.png", plot = plot1, width = 10, height = 8, dpi = 400)

Sup1<-ggarrange(Sup1A, Sup1B, ncol=2, common.legend = T)

#ggsave(file = "Figures/Q1_Alpha_Sup1.pdf", plot = Sup1, width = 12, height = 8, dpi = 400)
#ggsave(file = "Figures/Q1_Alpha_Sup1.png", plot = Sup1, width = 12, height = 8, dpi = 400)

##With Phylogenetic diversity
alphadiv.pig%>%
  dplyr::filter(InfectionStatus== "Infected")%>%
  dplyr::filter(Replicate%in%Inf.Keep)%>%
  wilcox_test(PD ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Infected_Compartment_PD.csv")

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
  ggplot(aes(x= Compartment, y= PD))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_line(aes(group = System), colour= "gray")+
  geom_point(shape= 21, size=3, aes(fill= System), color= "black")+
  scale_fill_manual(values = pal.system)+
  xlab("GI compartment")+
  ylab("Faiths phylogenetic diverstiy")+
  labs(tag= "B)", 
       shape = "Infection status", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+ #, caption = get_pwc_label(stats.test)
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -1000, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)

###Logistic regression 
alphadiv.pig%>%
  dplyr::mutate(InfectionStatus = case_when(InfectionStatus == "Infected"  ~ 1,
                                            InfectionStatus == "Non_infected" ~ 0))-> tmp

log.model.pig <- glm(InfectionStatus ~ Chao1, data = tmp, family = binomial)
summary(log.model.pig)$coef

alphadiv.pig%>%
  dplyr::mutate(Infection = case_when(InfectionStatus == "Infected"  ~ 1,
                                            InfectionStatus == "Non_infected" ~ 0))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment%in%c("Cecum"))%>%
  ggplot(aes(Chao1, Infection)) +
  geom_point(size=3, aes(color= InfectionStatus, fill= InfectionStatus), shape=21, color= "black")+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  xlab("ASV Richness (Chao1 Index)")+
  ylab("Infection status")+
  labs(tag= "B)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16))+
  scale_shape_manual(values = c(21, 24, 22))+
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = F)

##Is there any difference at Phylum level between infected and not infected?
#Agglomerate to phylum-level and rename
PS.pig.Phy <- phyloseq::tax_glom(PS.pig.Norm, "Phylum")
phyloseq::taxa_names(PS.pig.Phy) <- phyloseq::tax_table(PS.pig.Phy)[, "Phylum"]
phyloseq::otu_table(PS.pig.Phy)[, 1:15]
##Fusobacteriota, Patascibacteria, Planctomycetota, Synergistota have low counts.
##Subset high just for better plotting

PS.subset <- subset_taxa(PS.pig.Phy, rownames(tax_table(PS.pig.Phy)) %in% c("Bacteroidota",
                                                                            "Firmicutes",  "Actinobacteriota",
                                                                            "Proteobacteria"))
##Changes by compartment
phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(y.position= c(55, 60, 65, 70, 75, 80,
  110, 115, 120, 125, 130, 100))-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Compartments.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>%
  ggplot(data = ., aes(x = Compartment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = Compartment, shape= InfectionStatus), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(24,25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = pal.compartment)+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Relative Abundance (%)", shape = "Infection status") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)+
  facet_wrap(~ OTU, scales = "free")-> A

##Just in site of infection 
##Changes by compartment
phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment=="Jejunum")%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Jejunum_Infection.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment=="Jejunum")%>%
  mutate(Abundance = (Abundance/1E6)*100)%>%
  ggplot(data = ., aes(x = InfectionStatus, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = System), shape= 21, height = 0, width = .2, size= 3, color= "black") +
  scale_fill_manual(values = pal.system)+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Abundance (Sequencing reads)", fill= "Individual") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  facet_wrap(~ OTU, scales = "free")+
  scale_x_discrete(labels=c("Infected" = "Infected", 
                            "Non_infected" = "Non infected"))-> B

plot2<-plot_grid(A, B, ncol=1, align="v")

#ggsave(file = "Figures/Q1_Phylum_Compartment_Infection.pdf", plot = plot2, width = 12, height = 10, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Compartment_Infection.png", plot = plot2, width = 12, height = 10, dpi = 450)

###Save them individually
#ggsave(file = "Figures/Q1_Phylum_Compartment.png", plot = A, width = 12, height = 10, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Infection.png", plot = B, width = 12, height = 10, dpi = 450)

### Compared to Non-infected 
alphadiv.pig%>%
  dplyr::filter(InfectionStatus== "Non_infected")%>%
  #dplyr::filter(Replicate%in%Inf.Keep)%>%
  wilcox_test(Chao1 ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Non_Infected_Compartment.csv")

alphadiv.pig%>% 
  dplyr::filter(InfectionStatus== "Non_infected")%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  wilcox_effsize(Chao1 ~ Compartment)

##Plot 
alphadiv.pig%>%
  dplyr::filter(InfectionStatus== "Non_infected")%>%
  #dplyr::filter(Replicate%in%Inf.Keep)%>%
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
  labs(tag= "C)", caption = get_pwc_label(stats.test), 
       shape = "Infection status", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -1000, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)

##Comparison Alpha diversity between worms and site of infection

alphadiv.PA%>%
  dplyr::filter(Compartment%in% c("Jejunum","Ascaris"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  wilcox_test(Chao1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Infected_Jej_Ascaris_all.csv")

alphadiv.PA%>%
  dplyr::filter(Compartment%in% c("Jejunum","Ascaris"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  wilcox_effsize(Chao1 ~ InfectionStatus)

##Plot 
alphadiv.PA%>%
  dplyr::filter(Compartment%in% c("Jejunum", "Ascaris"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= InfectionStatus, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, aes(fill= System, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 25, 21), labels = c("Infected", " Non infected", "Worm"))+
  scale_fill_manual(values = pal.system)+
  xlab("Infection status")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "C)", caption = get_pwc_label(stats.test), 
       shape = "Infection status", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -0.2, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> C


D<- grid.arrange(plot1,C, widths = c(4, 2.5),
                 layout_matrix = rbind(c(1, 2)))

ggsave(file = "Figures/Q1_Alpha_Compartment.pdf", plot = D, width = 12, height = 8, dpi = 600)
ggsave(file = "Figures/Q1_Alpha_Compartment.png", plot = D, width = 12, height = 8, dpi = 600)

##Is there any difference at Phylum level between infected and not infected?
#Agglomerate to phylum-level and rename
PS.PA.Phy <- phyloseq::tax_glom(PS.PA.Norm, "Phylum")
phyloseq::taxa_names(PS.PA.Phy) <- phyloseq::tax_table(PS.PA.Phy)[, "Phylum"]
phyloseq::otu_table(PS.PA.Phy)[, 1:15]
##Fusobacteriota, Patascibacteria, Planctomycetota, Synergistota have low counts.
##Subset high just for better plotting

PS.subset <- subset_taxa(PS.PA.Phy, rownames(tax_table(PS.PA.Phy)) %in% c("Bacteroidota",
                                                                            "Firmicutes",  "Actinobacteriota",
                                                                            "Proteobacteria"))

PS.subset <- subset_samples(PS.subset, InfectionStatus%in%c("Infected", "Worm"))

##Changes by compartment
phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(y.position= c(55, 60, 65, 70, 75, 80, 85, 90,
                              110, 115, 120))-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Compartments_Ascaris.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>%
  ggplot(data = ., aes(x = Compartment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = Compartment, shape= InfectionStatus), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(24, 21), labels = c("Infected", "Worm"))+
  scale_fill_manual(values = pal.compartment)+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Relative Abundance (%)", shape = "Host-Parasite") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                    tip.length = 0)+
  facet_wrap(~ OTU, scales = "free")-> A

##Just in site of infection vs worm
PS.subset <- subset_taxa(PS.PA.Phy, rownames(tax_table(PS.PA.Phy)) %in% c("Bacteroidota",
                                                                          "Firmicutes",  "Actinobacteriota",
                                                                          "Proteobacteria"))
phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment%in%c("Jejunum","Ascaris"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Jejunum_Ascaris.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment%in%c("Jejunum","Ascaris"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>%
  ggplot(data = ., aes(x = InfectionStatus, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = System), shape= 21, height = 0, width = .2, size= 3, color= "black") +
  scale_fill_manual(values = pal.system)+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Abundance (Sequencing reads)", fill= "Individual") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  facet_wrap(~ OTU, scales = "free")+
  scale_x_discrete(labels=c("Infected" = "Jejunum (Infected)", 
                            "Non_infected"= "Jejunum (Non infected)",
                            "Worm" = "Ascaris"))-> B
##Check with duodenum
phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment%in%c("Duodenum","Ascaris"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Duodenum_Ascaris.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment%in%c("Duodenum","Ascaris"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>%
  ggplot(data = ., aes(x = InfectionStatus, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = System), shape= 21, height = 0, width = .2, size= 3, color= "black") +
  scale_fill_manual(values = pal.system)+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Abundance (Sequencing reads)", fill= "Individual") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  facet_wrap(~ OTU, scales = "free")+
  scale_x_discrete(labels=c("Infected" = "Duodenum (Infected)", 
                            "Non_infected"= "Duodenum (Non infected)",
                            "Worm" = "Ascaris"))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> C

###Save them individually
#ggsave(file = "Figures/Q1_Phylum_Compartment_Ascaris.png", plot = A, width = 12, height = 10, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Jejunum_Ascaris.png", plot = B, width = 12, height = 10, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Duodenum_Ascaris.png", plot = C, width = 12, height = 10, dpi = 450)

plot3<-plot_grid(A, B, C, ncol=1, align="v")

#ggsave(file = "Figures/Q1_Phylum_Compartment_Infection_Ascaris.pdf", plot = plot3, width = 15, height = 14, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Compartment_Infection_Ascaris.png", plot = plot3, width = 15, height = 14, dpi = 450)

##Comparison Alpha diversity between worms experiments and Slaughterhouse
foo<- estimate_pd(PS.Asc)

alphadiv.Asc<- cbind(alphadiv.Asc, foo)

alphadiv.Asc%>%
  mutate(Origin = fct_relevel(Origin, 
                                   "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  wilcox_test(Chao1 ~ Origin)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Origin", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Infected_Ascaris_all.csv")

alphadiv.Asc%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  wilcox_effsize(Chao1 ~ Origin)

##Plot 
alphadiv.Asc%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  ggplot(aes(x= Origin, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, aes(fill= System, shape=WormSex), color= "black")+
  scale_shape_manual(values = c(23,22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  xlab("Origin")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test), 
       shape = "Worm Sex", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  scale_x_discrete(labels=c("Experiment_1" = "Exp. 1", 
                            "Experiment_2" = "Exp. 2",
                            "Slaughterhouse" = "Slaughterhouse"))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -0.2, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> A

###Sex difference 
alphadiv.Asc%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::group_by(Origin)%>%
  wilcox_test(Chao1 ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Sex_Ascaris_Origin.csv")

# New facet label names for supp variable
supp.labs <- c("Exp. 1", "Exp. 2", "Slaughterhouse")
names(supp.labs) <- c("Experiment_1", "Experiment_2", "Slaughterhouse")

##Plot 
alphadiv.Asc%>%
  dplyr::group_by(Origin)%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3",
                              "Pig5","Pig10","Pig11", 
                              "Pig12", "Pig13", "Pig14", "SH"))%>%
  ggplot(aes(x= WormSex, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, aes(fill= System, shape=WormSex), color= "black")+
  scale_shape_manual(values = c(23,22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  facet_grid(~Origin, scales = "free", space = "free", labeller = labeller(Origin = supp.labs))+
  xlab("Worm sex")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test), 
       shape = "Worm Sex", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -0.2, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> B

###Pool all experiments toghether to assess sex difference 
alphadiv.Asc%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(Location = case_when(Origin %in% c("Experiment_1", "Experiment_2")  ~ "FU",
                                     Origin == "Slaughterhouse" ~ "SH"))-> alphadiv.Asc
alphadiv.Asc%>%
  dplyr::group_by(Location)%>%
  wilcox_test(Chao1 ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Sex_Ascaris_Location.csv")

# New facet label names for supp variable
supp.labs <- c("Local housing", "Slaughterhouse")
names(supp.labs) <- c("FU", "SH")

##Plot 
alphadiv.Asc%>%
  dplyr::group_by(Location)%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14", "SH"))%>%
  ggplot(aes(x= WormSex, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, aes(fill= System, shape=WormSex), color= "black")+
  scale_shape_manual(values = c(23,22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  facet_grid(~Location, scales = "free", space = "free", labeller = labeller(Location = supp.labs))+
  xlab("Worm sex")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "C)", caption = get_pwc_label(stats.test), 
       shape = "Worm Sex", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -0.2, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> C

###Pooling all worms independently from origin and location
alphadiv.Asc%>%
  wilcox_test(Chao1 ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Sex_Ascaris_all.csv")

alphadiv.Asc%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14", "SH"))%>%
  ggplot(aes(x= WormSex, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, aes(fill= System, shape=WormSex), color= "black")+
  scale_shape_manual(values = c(23,22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  xlab("Worm sex")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "D)", caption = get_pwc_label(stats.test), 
       shape = "Worm Sex", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -0.2, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> D

##Check for richness difference by location
alphadiv.Asc%>%
  wilcox_test(Chao1 ~ Location)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Location", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Ascaris_Location.csv") ##NS do not plot 

##Is there any difference at Phylum level between Sexes of worms and from different locations?
#Agglomerate to phylum-level and rename
PS.Asc.Norm.Phy <- phyloseq::tax_glom(PS.Asc.Norm, "Phylum")
phyloseq::taxa_names(PS.Asc.Norm.Phy) <- phyloseq::tax_table(PS.Asc.Norm.Phy)[, "Phylum"]
phyloseq::otu_table(PS.Asc.Norm.Phy)[, 1:15]
##Fusobacteriota, Patascibacteria, Planctomycetota, Synergistota have low counts.
##Subset high just for better plotting

PS.subset <- subset_taxa(PS.Asc.Norm.Phy, rownames(tax_table(PS.Asc.Norm.Phy)) %in% c("Bacteroidota",
                                                                          "Firmicutes",  "Actinobacteriota",
                                                                          "Proteobacteria"))
tmp<- row.names(PS.subset@sam_data)
PS.subset@sam_data<- sample_data(alphadiv.Asc[rownames(alphadiv.Asc)%in%tmp, ])

##Changes by housing, experiment, sex
#1) Experiments and SH
phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                                     "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                       "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ Origin)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Origin")%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
    dplyr::mutate(y.position= c(11, 12, 40, 45, 111, 116, 51, 101))-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Ascaris_Origin.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  ggplot(data = ., aes(x = Origin, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = System, shape= WormSex), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Relative Abundance (%)", shape = "Worm Sex", fill= "Host") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)+
  scale_x_discrete(labels=c("Experiment_1" = "Exp. 1", 
                            "Experiment_2" = "Exp. 2",
                            "Slaughterhouse" = "Slaughterhouse"))+
  facet_wrap(~ OTU, scales = "free")-> E

#2) Location FU vs SH
phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ Location)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Location")%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(y.position= c(12, 40, 111, 101))-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Ascaris_Location.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  ggplot(data = ., aes(x = Location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = Location, shape= WormSex), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  scale_fill_manual(values = c("#84BD00FF", "#BB0021FF"), labels = c("Local housing", "Slaughterhouse"))+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+
  labs(x = "", y = "Relative Abundance (%)", shape = "Worm Sex", fill= "Location of Host") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)+
  facet_wrap(~ OTU, scales = "free")-> G

#3) Sex all worms pooled
phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Ascaris_Sex_All.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  ggplot(data = ., aes(x = WormSex, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = WormSex, shape= WormSex), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  #scale_fill_manual(values = pal.system)+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16),  axis.title.x=element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+
  labs(x = "", y = "Relative Abundance (%)",fill= "Worm Sex") +
  guides(fill = guide_legend(override.aes=list(shape=c(23, 22))), color= FALSE, shape= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)+
  facet_wrap(~ OTU, scales = "free")-> H

#3.1) Sex FU worms
phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  filter(Location== "FU")%>%
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Ascaris_Sex_FU.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  filter(Location== "FU")%>%
  dplyr::group_by(OTU)%>%
  ggplot(data = ., aes(x = WormSex, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = WormSex, shape= WormSex), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16),  axis.title.x=element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+
  labs(x = "", y = "Relative Abundance (%)",fill= "Worm Sex") +
  guides(fill = guide_legend(override.aes=list(shape=c(23, 22))), color= FALSE, shape= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)+
  facet_wrap(~ OTU, scales = "free")-> Asc.FU.sex

#3.2) Sex SH worms
phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  filter(Location== "SH")%>%
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_Phylum_Ascaris_Sex_SH.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  mutate(Location = fct_relevel(Location, 
                                "FU", "SH"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  filter(Location== "SH")%>%
  dplyr::group_by(OTU)%>%
  ggplot(data = ., aes(x = WormSex, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = WormSex, shape= WormSex), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16),  axis.title.x=element_blank(), axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+
  labs(x = "", y = "Relative Abundance (%)",fill= "Worm Sex") +
  guides(fill = guide_legend(override.aes=list(shape=c(23, 22))), color= FALSE, shape= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)+
  facet_wrap(~ OTU, scales = "free")-> I

##Save the individual plots
#ggsave(file = "Figures/Q1_Alpha_Worms_Origin.png", plot = A, width = 10, height = 8, dpi = 450)
#ggsave(file = "Figures/Q1_Alpha_Worms_Sex_Origin.png", plot = B, width = 10, height = 8, dpi = 450)
#ggsave(file = "Figures/Q1_Alpha_Worms_Sex_Location.png", plot = C, width = 10, height = 8, dpi = 450)
#ggsave(file = "Figures/Q1_Alpha_Worms_Sex.png", plot = D, width = 10, height = 8, dpi = 450)

#ggsave(file = "Figures/Q1_Phylum_Ascaris_Origin.png", plot = E, width = 12, height = 10, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Ascaris_Location.png", plot = G, width = 12, height = 10, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Ascaris_Sex_all.png", plot = H, width = 12, height = 10, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Ascaris_Sex_SH.png", plot = I, width = 12, height = 10, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Ascaris_Sex_FU.png", plot = Asc.FU.sex, width = 12, height = 10, dpi = 450)

###Al together
Plot1<- ggarrange(A,B,C,D, ncol=2, nrow=2, common.legend = TRUE, legend="right")

#ggsave(file = "Figures/Q1_Alpha_Worm.pdf", plot = Plot1, width = 12, height = 10, dpi = 450)
#ggsave(file = "Figures/Q1_Alpha_Worm.png", plot = Plot1, width = 12, height = 10, dpi = 450)

E+
  guides(shape= F)-> E

I+
  guides(fill= F)-> I

plot4<-plot_grid(E, G, H, I, ncol=1, align="v", axis = "lr")

#ggsave(file = "Figures/Q1_Phylum_Ascaris.pdf", plot = plot4, width = 12, height = 22, dpi = 450)
#ggsave(file = "Figures/Q1_Phylum_Ascaris.png", plot = plot4, width = 12, height = 22, dpi = 450)