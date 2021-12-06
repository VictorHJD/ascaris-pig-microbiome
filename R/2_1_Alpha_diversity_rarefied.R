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

##Alpha diverisity tables with sample information
alphadiv.PA<- readRDS("Data/alphadiv.PA.rds")
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

alphadiv.pig%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_effsize(Chao1 ~ InfectionStatus)

##Plot 
alphadiv.PA.rare%>%
  dplyr::filter(InfectionStatus!= "Worm")%>%
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
  guides(fill = F, color= FALSE)+
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

full.model<- glm(Chao1 ~ InfectionStatus * Compartment, data = tmp) ##Full model
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
  labs(title = NULL, tag= "A)")+
  theme_classic()+
  theme(text = element_text(size=16))

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

###Logistic regression 
alphadiv.PA.rare%>%
  dplyr::filter(InfectionStatus!= "Worm")%>%
  dplyr::mutate(InfectionStatus = case_when(InfectionStatus == "Infected"  ~ 1,
                                            InfectionStatus == "Non_infected" ~ 0))-> tmp

log.model.pig <- glm(InfectionStatus ~ Chao1, data = tmp, family = binomial)
summary(log.model.pig)$coef

alphadiv.PA.rare%>%
  dplyr::filter(InfectionStatus!= "Worm")%>%
  dplyr::mutate(Infection = case_when(InfectionStatus == "Infected"  ~ 1,
                                      InfectionStatus == "Non_infected" ~ 0))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::filter(Compartment%in%c("Jejunum"))%>%
  ggplot(aes(Chao1, Infection)) +
  geom_point(size=3, aes(color= InfectionStatus, fill= InfectionStatus), shape=21, color= "black")+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values =  c("#D55E00", "#009E73"), labels = c("Infected", "Non infected"))+
  xlab("ASV Richness (Chao1 Index)")+
  ylab("Infection status")+
  labs(tag= "B)", fill= "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16))+
  scale_shape_manual(values = c(21, 24, 22))+
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = F)

###Ranked analysis 
#Just Infected with points for all compartments
alphadiv.PA.rare%>%
  dplyr::group_by(System)%>%
  dplyr::filter(n()==5)%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

alphadiv.PA.rare%>%
  dplyr::filter(Replicate%in%Inf.Keep)%>%
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                          "Duodenum", "Jejunum", "Ileum", 
                                          "Cecum", "Colon"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3","Pig4",
                                     "Pig5","Pig6","Pig7","Pig8","Pig9",
                                     "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::group_by(Compartment)%>%
  dplyr::mutate(Rank= rank(Chao1, ties.method = "first"))-> tmp

alphadiv.PA.rare%>%
  dplyr::filter(Replicate%in%Inf.Keep)%>%
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                          "Duodenum", "Jejunum", "Ileum", 
                                          "Cecum", "Colon"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3","Pig4",
                                     "Pig5","Pig6","Pig7","Pig8","Pig9",
                                     "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  dplyr::group_by(Compartment)%>%
  dplyr::mutate(Rank= rank(Shannon, ties.method = "first"))%>%
  ggplot(aes(x= Compartment, y= Rank))+
  #geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_line(aes(group = System), colour= "gray")+
  geom_point(shape= 21, size=3, aes(fill= System), color= "black")+
  scale_fill_manual(values = pal.system)+
  xlab("GI compartment")+
  ylab("Diversity rank (Chao1 Index)")+
  labs(tag= "B)", 
       shape = "Infection status", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+ #, caption = get_pwc_label(stats.test)
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())

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