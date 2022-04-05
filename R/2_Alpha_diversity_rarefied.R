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
library(lme4)

##Load data 
PS.PA<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Rds") ## Data merged pigs and Ascaris (not SH) not normalized for alpha diversity plots 
PS.PA.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Norm.Rds") ## Data merged pigs and Ascaris (not SH) normalized for beta diversity plots 
PS.PA.rare<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.rare.Rds") ## Data merged pigs and Ascaris (not SH) rarefied for alpha diversity plots 

##Alpha diversity tables with sample information
alphadiv.PA.rare<- readRDS("Data/alphadiv.PA.rare.rds")

##Color palette for compartment and system ##

pal.compartment <- c("Ascaris"="#1B9E77","Cecum"= "#D95F02","Colon"= "#7570B3",
                     "Duodenum"= "#E7298A","Ileum"= "#66A61E","Jejunum"="#E6AB02")

pal.system <- c("Pig1"= "#A6761D","Pig2"= "#666666","Pig3"= "#A6CEE3","Pig4"= "#1F78B4",
                "Pig5"= "#B2DF8A","Pig6"= "#33A02C","Pig7"= "#FB9A99","Pig8"="#E31A1C","Pig9"= "#FDBF6F",
                "Pig10"= "#FF7F00","Pig11"= "#CAB2D6","Pig12"= "#6A3D9A","Pig13"= "#FFFF99",  "Pig14"= "#3B3B3BFF", "SH" = "#BB0021FF")

pal.infection<- c("Infected"= "#D55E00","Non_infected" = "#009E73")

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

###Add some additional stuff
## 1)Phylogenetic richness
require("btools")
foo<- estimate_pd(PS.PA.rare)
alphadiv.PA.rare<- cbind(alphadiv.PA.rare, foo)

## 2) Original library size
alphadiv.PA.rare$Library_Size<- sample_sums(PS.PA)

### Infected vs Non Infected (Richness)
alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus!= "Worm")%>%
  #dplyr::mutate(InfectionStatus = case_when(System == "Pig4"  ~ "Non_infected",
  #                                          TRUE ~ as.character(InfectionStatus)))%>%
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
write.csv(x, "Tables/Q1_Alpha_Rare_Infection_Richness.csv")

##Plot 
alphadiv.PA.rare%>%
  dplyr::filter(InfectionStatus!= "Worm")%>%
  #dplyr::mutate(InfectionStatus = case_when(System == "Pig4"  ~ "Non_infected",
  #                                          TRUE ~ as.character(InfectionStatus)))%>%
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= Chao1, color= InfectionStatus, fill= InfectionStatus))+
  geom_boxplot(aes(),outlier.shape=NA)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = pal.infection)+
  xlab("GI compartment")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "B)", fill= "Infection status")+
  guides(fill = "none", color= "none")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), panel.border = element_blank())-> A

###GLM
##Model selection do it with glm
alphadiv.PA.rare%>%
  dplyr::filter(InfectionStatus!= "Worm")%>%
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                          "Duodenum", "Jejunum", "Ileum", 
                                          "Cecum", "Colon"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3","Pig4",
                                     "Pig5","Pig6","Pig7","Pig8","Pig9",
                                     "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

full.model<- glm(Chao1 ~ InfectionStatus * Compartment + Library_Size, data = tmp) ##Full model
summary(full.model)

# Stepwise regression model
step.model <- MASS::stepAIC(full.model, direction = "both", 
                            trace = FALSE)
summary(step.model)

##Infection status is not a significant predictor for alpha diversity 
tr0<- lmer(Chao1 ~ (1 | System), data = tmp) ##Null model
tr1<-lmer(Chao1 ~ InfectionStatus * Compartment + (1 | System), data = tmp)

summary(tr1)

require("lmtest")
lrtest(tr0, tr1)

#Likelihood ratio test

#Model 1: Chao1 ~ (1 | System)
#Model 2: Chao1 ~ InfectionStatus * Compartment + (1 | System)
#Df  LogLik Df  Chisq Pr(>Chisq)    
#1   3 -375.04                         
#2  12 -277.71  9 194.66  < 2.2e-16 ***

##Plot model 
##For analysis with linear models
require("merTools")
est.plot<- plotREsim(REsim(tr1))  ## plot the interval estimates
est.plot$data%>%
  dplyr::mutate(groupID = fct_relevel(groupID, 
                                      "Pig1","Pig2","Pig3","Pig4",
                                      "Pig5","Pig6","Pig7","Pig8","Pig9",
                                      "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> est.plot$data
est.plot+
  geom_point(shape= 21, size=2.5, aes(fill= groupID), color= "black")+
  scale_fill_manual(values = pal.system)+
  xlab(label = NULL)+
  ylab(label = "Estimates")+
  labs(title = NULL, tag= "A)", fill= "Individual")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=16))-> est.plot

##PLot model  
require("sjPlot")
plot_model(tr1, p.adjust = "BH", vline.color = "gray", show.p = T, sort.est = TRUE)+
  geom_point(shape= 21, size=2.5, aes(fill= group), color= "black")+
  labs(title = NULL, tag= "B)")+
  theme_classic()+
  theme(text = element_text(size=16))-> mod.plot

## Infected vs Non Infected (Shannon Diversity)
alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus!= "Worm")%>%
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
write.csv(x, "Tables/Q1_Alpha_Rare_Infection_Shannon.csv")

##Plot 
alphadiv.PA.rare%>%
  dplyr::filter(InfectionStatus!= "Worm")%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= Shannon, color= InfectionStatus, fill= InfectionStatus))+
  geom_boxplot(aes(),outlier.shape=NA)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = pal.infection, labels = c("Infected", "Non infected"))+
  xlab("GI compartment")+
  ylab("ASV Diversity (Shannon Index)")+
  labs(tag= "A)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  scale_y_continuous(limits=c(0, 4.5))+
  annotate("text", x = 5, y = 4.4, label = '"*"', parse = TRUE)+
  annotate("segment", x = 4.8, xend = 5.2, y = 4.3, yend = 4.3, colour = "black")-> Sup1A

##Infected vs Non Infected (Phylogenetic Diversity) 
alphadiv.PA.rare%>%
  dplyr::filter(InfectionStatus!= "Worm")%>%
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
write.csv(x, "Tables/Q1_Alpha_Rare_Infection_PD.csv")

##Plot 
alphadiv.PA.rare%>%
  dplyr::filter(InfectionStatus!= "Worm")%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= PD, color= InfectionStatus, fill= InfectionStatus))+
  geom_boxplot(aes(),outlier.shape=NA)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = pal.infection, labels = c("Infected", "Non infected"))+
  xlab("GI compartment")+
  ylab("Phylogenetic diverstiy (Faith's Index)")+
  labs(tag= "B)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  scale_y_continuous(limits=c(100, 850))-> Sup1B

Sup1<-ggarrange(Sup1A, Sup1B, ncol=2, common.legend = T)

#ggsave(file = "Figures/Q1_Alpha_Rare_Sup1.pdf", plot = Sup1, width = 12, height = 8, dpi = 400)
#ggsave(file = "Figures/Q1_Alpha_Rare_Sup1.png", plot = Sup1, width = 12, height = 8, dpi = 400)

alphadiv.PA.rare%>%
  dplyr::filter(InfectionStatus!= "Worm")-> alphadiv.pig.rare

wilcox.test(Chao1 ~ InfectionStatus, data = alphadiv.pig.rare, exact = FALSE, conf.int = TRUE)

##Compare alpha diversity for sections and ascaris 
alphadiv.PA.rare%>%
  dplyr::filter(Compartment%in% c("Jejunum","Ascaris"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  wilcox_test(Chao1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test 

##Manual adjustment of position 
stats.test%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(xmin= c(1.8, 2.2))%>%
  dplyr::mutate(xmax= c(3, 3))-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Jej_Ascaris_rare.csv")

##Plot 
alphadiv.PA.rare%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ascaris","Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment%in%c("Jejunum", "Ascaris"))%>%
  ggplot(aes(x= Compartment, y= Chao1, color= InfectionStatus, fill= InfectionStatus))+
  geom_boxplot(aes(),outlier.shape=NA)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black", "black"))+
  scale_fill_manual(values = c("#D55E00","#009E73","#E69F00"), labels = c("Infected", "Non infected", "Ascaris"))+
  xlab("GI compartment")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "C)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  scale_y_continuous(limits=c(0, 250))+
  annotate("text", x = 1.6, y = 210, label = '"*"', parse = TRUE)+
  annotate("text", x = 1.4, y = 228, label = '"***"', parse = TRUE)+
  annotate("segment", x = 0.8, xend = 2, y = 226, yend = 226, colour = "black")+
  annotate("segment", x = 1.2, xend = 2, y = 208, yend = 208, colour = "black")-> B

fig.1<- grid.arrange(A,B, widths = c(4, 4, 4, 4),
                     layout_matrix = rbind(c(NA, NA, 2, 2),
                                           c( 1, 1, 2, 2)))

ggsave(file = "Figures/Q1_Alpha_Compartment_rare.pdf", plot = fig.1, width = 12, height = 8, dpi = 600)
ggsave(file = "Figures/Q1_Alpha_Compartment_rare.png", plot = fig.1, width = 12, height = 8, dpi = 600)
ggsave(file = "Figures/Q1_Alpha_Compartment_rare.svg", plot = fig.1, width = 12, height = 8, dpi = 600)

###Further adjustments in inkscape 

##With Shannon diversity 
Sup1A+
  guides(fill= F)-> Sup1A

alphadiv.PA.rare%>%
  dplyr::filter(Compartment%in% c("Jejunum","Ascaris"))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  wilcox_test(Shannon ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus") ##No difference

##Plot 
alphadiv.PA.rare%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ascaris","Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment%in%c("Jejunum", "Ascaris"))%>%
  ggplot(aes(x= Compartment, y= Shannon, color= InfectionStatus, fill= InfectionStatus))+
  geom_boxplot(aes(),outlier.shape=NA)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black", "black"))+
  scale_fill_manual(values = c("#D55E00","#009E73","#E69F00"), labels = c("Infected", "Non infected", "Ascaris"))+
  xlab("GI compartment")+
  ylab("ASV Diversity (Shannon Index)")+
  labs(tag= "B)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  scale_y_continuous(limits=c(0, 4.5))-> Sup1.1A

fig.1.1<- grid.arrange(Sup1A,Sup1.1A, widths = c(4, 4, 4, 4),
                     layout_matrix = rbind(c(NA, NA, 2, 2),
                                           c( 1, 1, 2, 2)))

##Comparison Alpha diversity between worms experiments 
alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus== "Worm")%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2"))%>%
  wilcox_test(Chao1 ~ Origin)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Origin", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Infected_Ascaris_all.csv")

alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus== "Worm")%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2"))%>%
  wilcox_effsize(Chao1 ~ Origin)

##Plot 
set.seed(2021)
alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus== "Worm")%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Origin, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, aes(fill= System, shape=WormSex), color= "black")+
  scale_shape_manual(values = c(23,22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  xlab("Origin")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "A)", 
       shape = "Worm Sex", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  scale_x_discrete(labels=c("Experiment_1" = "Exp. 1", 
                            "Experiment_2" = "Exp. 2"))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -0.2, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> Sup2A

###Sex difference within experiments
alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus== "Worm")%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2"))%>%
  dplyr::group_by(Origin)%>%
  wilcox_test(Chao1 ~ WormSex)%>%
  #adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Sex_Ascaris_Origin.csv")

##Plot 
alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus== "Worm")%>%
  dplyr::group_by(Origin)%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                                     "Experiment_1", "Experiment_2"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= Chao1, color= WormSex, fill= WormSex))+
  geom_boxplot(aes(),outlier.shape=NA)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black"))+
  facet_grid(~Origin, scales = "free", space = "free")+
  xlab("Worm sex")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test), 
       fill= "Worm Sex")+
  guides(color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())-> Sup2B

###Supplementary figure alpha diversity worms
Sup2<- grid.arrange(Sup2A,Sup2B)

ggsave(file = "Figures/Sup_Alpha_worms_rare.pdf", plot = Sup2, width = 8, height = 10, dpi = 600)
ggsave(file = "Figures/Sup_Alpha_worms_rare.png", plot = Sup2, width = 8, height = 10, dpi = 600)
ggsave(file = "Figures/Sup_Alpha_worms_rare.svg", plot = Sup2, width = 8, height = 10, dpi = 600)

###Pool all experiments together to assess sex difference 

alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus== "Worm")%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2"))%>%
  wilcox_test(Chao1 ~ WormSex)%>%
  add_significance()%>%
  add_xy_position(x = "WormSex", dodge = 0.8)-> stats.test 

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Alpha_Sex_Ascaris_Location.csv")

##Plot 
set.seed(2021)
alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus== "Worm")%>%
  mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= WormSex, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(position=position_jitter(0.3), size=3, aes(fill= System, shape=WormSex), color= "black")+
  scale_shape_manual(values = c(23,22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  xlab("Worm sex")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test), 
       shape = "Worm Sex", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -0.2, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> WormsA.V1

set.seed(2021)
alphadiv.PA.rare%>% 
  dplyr::filter(InfectionStatus== "Worm")%>%
  dplyr::group_by(Origin)%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                                     "Experiment_1", "Experiment_2"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= Chao1, color= WormSex, fill= WormSex))+
  geom_boxplot(aes(),outlier.shape=NA)+
  geom_point(position = position_jitterdodge())+
  scale_color_manual(values = c("black", "black"))+
  xlab("Worm sex")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test), 
       fill= "Worm Sex")+
  guides(color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank())-> WormsA.V2

# save figures as rds for further composition with differential abundance
saveRDS(WormsA.V1, "Figures/Alpha_Sex_Worms_A_V1.RDS") ##with individual coloring 
saveRDS(WormsA.V2, "Figures/Alpha_Sex_Worms_A_V2.RDS") ##with sex coloring 

###Add worm burden data
x<- c(5, 187, 42, 0, 1, 0, 0, 0, 65, 61, 108, 14)
z<- c("Pig1.Jejunum", "Pig2.Jejunum", "Pig3.Jejunum", "Pig4.Jejunum", "Pig5.Jejunum", "Pig7.Jejunum",
      "Pig8.Jejunum", "Pig9.Jejunum", "Pig10.Jejunum", "Pig11.Jejunum", "Pig12.Jejunum", "Pig13.Jejunum")

foo<- data.frame(z,x, row.names = z)
colnames(foo)<- c("Replicate", "Worm_load")

foo%>%
  dplyr::mutate(Replicate = fct_relevel(Replicate, 
                                        "Pig1.Jejunum", "Pig2.Jejunum", "Pig3.Jejunum", 
                                        "Pig4.Jejunum", "Pig5.Jejunum", "Pig7.Jejunum",
                                        "Pig8.Jejunum", "Pig9.Jejunum", "Pig10.Jejunum", 
                                        "Pig11.Jejunum", "Pig12.Jejunum", "Pig13.Jejunum"))-> foo

alphadiv.pig.rare%>%
  dplyr::filter(Compartment=="Jejunum")-> x

plyr::join(x, foo, by= "Replicate")-> foo

###Is the alpha diversity in the site of infection linked to the worm load
foo%>%
  ggplot(aes(x= Worm_load, y= Chao1))+
  geom_point(size=3, aes(shape= InfectionStatus, fill= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = pal.infection, labels = c("Infected", "Non infected"))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F)+
  xlab("Worm load")+
  ylab("ASV Richness (Chao1 Index)")+
  labs(tag= "A)", fill= "Infection status")+
  theme_bw()+
  geom_smooth(method=lm, se = T, color= "black")+
  stat_cor(label.x = 145, label.y = 150,
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  theme(text = element_text(size=16), axis.title.x = element_blank())-> ch1.load

foo%>%
  ggplot(aes(x= Worm_load, y= Shannon))+
  geom_point(size=3, aes(shape= InfectionStatus, fill= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = pal.infection, labels = c("Infected", "Non infected"))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F)+
  xlab("Worm load")+
  ylab("Diversity (Shannon Index)")+
  labs(tag= "B)", fill= "Infection status")+
  theme_bw()+
  geom_smooth(method=lm, se = T, color= "black")+
  stat_cor(label.x = 145, label.y = 2.6,
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+ # Add correlation coefficient
  theme(text = element_text(size=16), axis.title.x = element_blank())-> sha.load

foo%>%
  ggplot(aes(x= Worm_load, y= PD))+
  geom_point(size=3, aes(shape= InfectionStatus, fill= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = pal.infection, labels = c("Infected", "Non infected"))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F)+
  xlab("Worm load")+
  ylab("Phylogenetic Diversity (Faith's Index)")+
  labs(tag= "C)", fill= "Infection status")+
  theme_bw()+
  geom_smooth(method=lm, se = T, color= "black")+
  stat_cor(label.x = 140, label.y = 350,
           aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  theme(text = element_text(size=16))-> PD.load

print(summary (lm (data = foo, Chao1 ~ Worm_load)))
print(summary (lm (data = foo, Shannon ~ Worm_load)))
print(summary (lm (data = foo, PD ~ Worm_load)))


Plot1<- ggarrange(ch1.load, sha.load, PD.load, ncol=1, align = "v", common.legend = T)

ggsave(file = "Figures/Sup_Alpha_load.pdf", plot = Plot1, width = 9, height = 12, dpi = 450)
ggsave(file = "Figures/Sup_Alpha_load.png", plot = Plot1, width = 9, height = 12, dpi = 450)
ggsave(file = "Figures/Sup_Alpha_load.svg", plot = Plot1, width = 9, height = 12, dpi = 450)

rm(ch1.load, sha.load, PD.load, Plot1)

##Core microbiome analysis
##For Pigs and worms
###Infected Non Infected Jejunum Ascairs
infst<- c("Infected", "Non_infected", "Worm")

list_core <- c() # an empty object to store information

for (n in infst){ # for each variable n in InfectionStatus
  
  tmp<- subset_samples(PS.PA.rare, InfectionStatus==n)
  tmp <- microbiome::transform(tmp, "compositional")
  core_m <- core_members(tmp, # ps.sub is phyloseq selected with only samples from g 
            detection = 0.0001, # 0.001 in atleast 90% samples 
            prevalence = 0.5)

  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each InfectionStatus.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

plot(eulerr::venn(list_core),
     fills = c("#D55E00","#009E73","#E69F00"),
      alpha= 0.5,
     labels = c("Infected", "Non infected", "Ascaris"))-> D

cowplot::plot_grid(
 D, labels = "D)", 
  label_fontfamily = "sans",
  label_fontface = "plain",
 label_size = 16)-> D

ggsave(file = "Figures/Q1_Alpha_Compartment_rare_D.pdf", plot = D, width = 8, height = 8, dpi = 600)
ggsave(file = "Figures/Q1_Alpha_Compartment_rare_D.png", plot = D, width = 8, height = 8, dpi = 600)
ggsave(file = "Figures/Q1_Alpha_Compartment_rare_D.svg", plot = D, width = 8, height = 8, dpi = 600)

# get the taxonomy data
tax.mat <- tax_table(PS.PA.rare)
tax.df <- as.data.frame(tax.mat)

# add the OTUs to last column
tax.df$ASV <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.piginf <- dplyr::filter(tax.df, rownames(tax.df) %in% list_core[["Infected"]])
core.taxa.pigNinf <- dplyr::filter(tax.df, rownames(tax.df) %in% list_core[["Non_infected"]])
core.taxa.worm <- dplyr::filter(tax.df, rownames(tax.df) %in% list_core[["Worm"]])

core.asv<- union(list_core[["Infected"]], list_core[["Non_infected"]])
core.asv<- union(core.asv, list_core[["Worm"]])

##Store "core" ASVs
write.csv(core.taxa.piginf, "Tables/Core_ASVs_Inf_Jejunum.csv")
write.csv(core.taxa.pigNinf, "Tables/Core_ASVs_NonInf_Jejunum.csv")
write.csv(core.taxa.worm, "Tables/Core_ASVs_Ascaris.csv")