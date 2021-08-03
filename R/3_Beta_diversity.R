##Project: Ascaris - Pig Microbiome
##Aim: Beta diversity
##Author: Víctor Hugo Jarquín-Díaz
##Root repo setwd("../Ascaris/ascaris/")

require(ggpubr)
require(RColorBrewer)
require(rstatix)
library(cowplot)
library(gridExtra)
library(grid)
library(doParallel)  
library(ggsci)
library(microbiome)
library(tidyverse)

##For taxa 
tax.palette<- c("Taxa less represented" = "#767676FF",  "Unassigned"="lightgray", "Prevotella"= "#3C5488FF", "Blautia" = "#AD002AFF",
                "Bacteroides" = "#00A087FF",  "Bifidobacterium" = "#E64B35FF", "Subdoligranulum"= "#F39B7FFF", "Faecalibacterium"= "#8491B4FF",   
                "Collinsella"= "#CD534CFF", "Alistipes" = "#FAFD7CFF", "Holdemanella"= "#7E6148FF","Lactobacillus"=  "#631879FF",
                "Ruminococcus"= "#BC3C29FF","Escherichia-Shigella" = "#0072B5FF", "Enterococcus" = "#E18727FF", "Roseburia"= "#E762D7FF",                  
                "Acidaminococcus"= "#7876B1FF", "Intestinibacter"="#6F99ADFF", "Butyricicoccus" = "#FFDC91FF", 
                "[Ruminococcus] gnavus group" = "#EE4C97FF","Clostridium sensu stricto 1"= "#00468BFF", "Agathobacter" = "#0099B4FF" , 
                "Lachnoclostridium"= "#42B540FF", "Pediococcus"= "#925E9FFF", "Streptococcus"= "#925E9FFF", "Staphylococcus"= "#008B45FF", 
                "Rothia" ="#FDAF91FF","Alloprevotella"= "#4DBBD5FF", "Veillonella"= "#3B4992FF", "Stenotrophomonas"= "#B09C85FF", 
                "Porphyromonas"= "#BB0021FF", "Granulicatella"= "#A20056FF","Gemella" = "#0073C2FF", 
                "Achromobacter"= "#EFC000FF" , "Pseudomonas" = "#ED0000FF", "Actinomyces" = "#91D1C2FF", 
                "Haemophilus" ="#7AA6DCFF"  ,  "Lachnoanaerobaculum"= "#003C67FF",   "Fusobacterium"= "#8F7700FF", 
                "Campylobacter"= "#A73030FF")

##Color palette for compartment and sytem ##

pal.compartment <- c("Ascaris"="#1B9E77","Cecum"= "#D95F02","Colon"= "#7570B3",
                     "Duodenum"= "#E7298A","Ileum"= "#66A61E","Jejunum"="#E6AB02")

pal.system <- c("Pig1"= "#A6761D","Pig2"= "#666666","Pig3"= "#A6CEE3","Pig4"= "#1F78B4",
                "Pig5"= "#B2DF8A","Pig6"= "#33A02C","Pig7"= "#FB9A99","Pig8"="#E31A1C","Pig9"= "#FDBF6F",
                "Pig10"= "#FF7F00","Pig11"= "#CAB2D6","Pig12"= "#6A3D9A","Pig13"= "#FFFF99",  "Pig14"= "#3B3B3BFF", "SH" = "#767676FF")

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

#########Question 1:
###How does Ascaris impact the porcine microbiome and does this differ in different gut regions?
###General comparison between infected and non infected pigs (all compartments, merged replicates) ###################
###Different compartments
bray_dist<- phyloseq::distance(PS.pig.Norm, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.pig.Norm,
                      method="PCoA", distance="bray")


tmp<- row.names(PS.pig.Norm@sam_data)
tmp<- alphadiv.pig[rownames(alphadiv.pig)%in%tmp, ]

tmp%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

compartment.adonis<- vegan::adonis(bray_dist~ Compartment + InfectionStatus + Origin,
                                permutations = 999, data = tmp, na.action = F, strata = tmp$System)

####
## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, tmp$InfectionStatus, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

plot(mvd)
anova(mvd)
boxplot(mvd)

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
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2,shape= grps), size=4, fill="black") + 
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  guides(shape= F)+
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, shape= InfectionStatus, fill= InfectionStatus), size=3) +
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  labs(tag= "A)", fill  = "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), color= F)+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> A

###Compare distances at PCo1 and PCo2 among groups
##PCo1
seg.data%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(v.PCoA1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Beta_Compartment_PCo1_Inf.csv")

##PCo2
seg.data%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(v.PCoA2 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Beta_Compartment_PCo2_Inf.csv")


## Non-metric multidimensional scaling
nmds.ordination<- ordinate(PS.pig.Norm,
                           method="NMDS", distance="bray")

nmds.scores<- as.data.frame(vegan::scores(nmds.ordination))
nmds.scores<- cbind(nmds.scores, tmp)

genus.scores<- as.data.frame(vegan::scores(nmds.ordination, "species"))
genus.data<- as.data.frame(PS.pig.Norm@tax_table)
genus.scores<- cbind(genus.scores, genus.data)

rm(genus.data)

find_hull<- function(df) df[chull(df$NMDS1, df$NMDS2), ]

nmds.scores%>%
  dplyr::filter(InfectionStatus== "Infected")-> nmds.scores.inf

nmds.scores%>%
  dplyr::filter(InfectionStatus!= "Infected")-> nmds.scores.Non.inf

hulls<- plyr::ddply(nmds.scores.Non.inf, "Compartment", find_hull)

nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  labs(tag= "A)", fill  = "Infection status")+
  scale_color_manual(values = c("#ED0000FF", "#008B45FF"))+
  geom_point(size=3, aes(fill= InfectionStatus, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  labs(tag= "B)", fill = "Infection status")+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F)+
  theme_bw()+
  theme(text = element_text(size=16))

##
mds.ordination <- ordinate(PS.pig.Norm, method="MDS", distance="euclidean")
eigen.vals <- mds.ordination$values$Eigenvalues 
# allows us to scale the axes according to their magnitude of separating apart the samples

##Compartments and Ascaris
bray_dist<- phyloseq::distance(PS.PA.Norm, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.PA.Norm,
                      method="PCoA", distance="bray")


tmp<- row.names(PS.PA.Norm@sam_data)

tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

pig.ascaris.adonis<- vegan::adonis(bray_dist~ Compartment + InfectionStatus + Origin, AnimalSpecies,
                                   permutations = 999, data = tmp, na.action = F, strata = tmp$System)

####
## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, tmp$Compartment, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

plot(mvd)
anova(mvd)
boxplot(mvd)
plot(TukeyHSD(mvd))

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

ggplot() + 
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2), shape= 21,size=4, fill="red") + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Compartment, shape= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24, 25, 21), labels= c("Infected", "Non infected", "Ascaris"))+
  scale_fill_manual(values = pal.compartment)+
  scale_color_manual(values = pal.compartment)+
  labs(tag= "B)", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> B

##Subset distances for infected pigs and their worms
tmp%>%
  dplyr::filter(InfectionStatus!= "Non_infected")%>%
  dplyr::filter(Compartment%in% c("Ascaris", "Jejunum"))%>%
  dplyr::filter(System!= "Pig14")%>%
  dplyr::group_by(System)%>%
  dplyr::filter(n()>=2)%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.JejAsc<- subset_samples(PS.PA.Norm, Replicate%in%Inf.Keep)

bray_dist<- phyloseq::distance(PS.JejAsc, 
                               method="bray", weighted=F)

ordination<- ordinate(PS.JejAsc,
                      method="PCoA", distance="bray")

tmp<- row.names(PS.JejAsc@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

jejunum.ascaris.adonis<- vegan::adonis(bray_dist~ Compartment + InfectionStatus + Origin, AnimalSpecies,
                                   permutations = 999, data = tmp, na.action = F, strata = tmp$System)

mvd<- vegan::betadisper(bray_dist, tmp$Compartment, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

plot(mvd)
anova(mvd)
boxplot(mvd)

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

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= System, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Pig", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  scale_color_manual(values = pal.system)+
  labs(tag= "C)", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> C

##Extract PCA axis
Vectors.PCA<- as.data.frame(ordination$vectors)
##Merge to the metadata 
Vectors.PCA<- cbind(Vectors.PCA, tmp)

###Create Boxplot to compare distances at PCo1 and PCo2 among groups
seg.data%>%
  wilcox_test(v.PCoA1 ~ Compartment)%>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_JejAsc_PCo1.csv")

seg.data%>%
  ggplot(aes(x=Compartment, y= v.PCoA1))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_point(size=3, aes(shape=InfectionStatus, fill= System)) +
  scale_shape_manual(values = c(24, 21), labels = c("Pig", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p} {p.signif}")+
  ylab("Bray-Curtis distance (PCo 1)")+
  labs(tag= "D)",  shape = "Host-Parasite", fill= "Individual")+ ##caption = get_pwc_label(stats.test), add in 
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())-> D

# Boxplot PCo2
seg.data%>%
  wilcox_test(v.PCoA2 ~ Compartment)%>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_JejAsc_PCo2.csv")

seg.data%>%
  ggplot(aes(x=Compartment, y= v.PCoA2))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_point(size=3, aes(shape=InfectionStatus, fill= System)) +
  scale_shape_manual(values = c(24, 21), labels = c("Pig", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p} {p.signif}")+
  ylab("Bray-Curtis distance (PCo 2)")+
  labs(tag= "E)",  shape = "Host-Parasite", fill= "Individual")+ ##caption = get_pwc_label(stats.test), add in 
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())-> E

D<- ggarrange(D, E, nrow = 1, common.legend = T)

E<- grid.arrange(A,B,C,D,
                 widths = c(4, 3, 3),
                 layout_matrix = rbind(c(1, 2, 2),
                                       c(3, 3, 3)))

ggsave(file = ".pdf", plot = E, width = 10, height = 8)
ggsave(file = ".png", plot = E, width = 10, height = 8)

rm(A,B,C, D, E)

##Extract PCA axis
Vectors.PCA<- as.data.frame(ordination$vectors)
##Merge to the metadata 
Vectors.PCA<- cbind(Vectors.PCA, tmp)

# Marginal density plot of x (top panel) and y (right panel)

xdens <- axis_canvas(B, axis = "x")+
  geom_density(data = Vectors.PCA, aes(x = Axis.1, fill = Compartment),
               alpha = 0.7, size = 0.2)+
  scale_fill_manual(values = pal.compartment)

ydens <- axis_canvas(B, axis = "y", coord_flip = TRUE)+
  geom_density(data = Vectors.PCA, aes(x = Axis.2, fill = Compartment),
               alpha = 0.7, size = 0.2)+
  coord_flip()+
  scale_fill_manual(values = pal.compartment)

p1 <- insert_xaxis_grob(B, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
B<- ggdraw(p2)

seg.data%>%
  dplyr::group_by(Compartment)%>%
  ggplot(aes(x=InfectionStatus, y= v.PCoA1))+
  geom_boxplot(color= "black", alpha= 0.5)+
  #geom_line(aes(group = System), colour= "gray")+
  geom_jitter(size=3, jitter= 0.5, aes(shape=InfectionStatus, fill= Compartment)) +
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = pal.compartment)+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p.adj.signif}")+
  xlab("Infection Status")+
  ylab("Bray-Curtis distance (PCo 1)")+
  labs(tag= "B)",  shape = "Infection status")+ ##caption = get_pwc_label(stats.test), add in 
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  facet_wrap(~InfectionStatus)+ 
  theme(text = element_text(size=16), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())-> C

# Boxplot PCo2
seg.data%>%
  dplyr::group_by(InfectionStatus)%>%
  wilcox_test(v.PCoA2 ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Beta_Compartment_PCo2_dist.csv")

seg.data%>%
  ggplot(aes(x=Compartment, y= v.PCoA2))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_line(aes(group = System), colour= "gray")+
  geom_point(size=3, aes(shape=InfectionStatus, fill= Compartment)) +
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = pal.compartment)+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p.adj.signif}")+
  xlab("Gastrointesinal compartment")+
  ylab("Bray-Curtis distance (PCo 2)")+
  labs(tag= "C)",  shape = "Infection status")+ ##caption = get_pwc_label(stats.test), add in 
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  facet_wrap(~InfectionStatus)+ 
  theme(text = element_text(size=16), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())-> D

D<- ggarrange(C, D, nrow = 1, common.legend = T)

E<- grid.arrange(B,D,
                 widths = c(3, 3),
                 layout_matrix = rbind(c(1, 1),
                                       c(2, 2)))

ggsave(file = "Figures/Q1_Composition_Compartments.pdf", plot = E, width = 10, height = 8, dpi = 650)
ggsave(file = "Figures/Q1_Composition_Compartments.png", plot = E, width = 10, height = 8,  dpi = 650)
