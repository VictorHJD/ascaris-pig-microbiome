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

##Load function
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
                "Pseudoscardovia"= "#BB0021FF")

##Color palette for compartment and system ##

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

##Store the result
foo<- as.data.frame(compartment.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Pig_Infection.csv")

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

find_hull<- function(df) df[chull(df$v.PCoA1, df$v.PCoA2), ]

hulls<- plyr::ddply(seg.data, "Origin", find_hull)

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, shape= InfectionStatus, fill= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F)+
  #geom_polygon(data= hulls, aes(x=v.PCoA1,y=v.PCoA2, color= Origin), alpha= 0.1)+
  #scale_color_manual(values = c("#4682B4", "#B47846"), labels= c("Experiment 1", "Experiment 2"))+
  labs(tag= "A)", fill  = "Infection status", color= "Origin of samples")+
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
##With phyloseq
nmds.ordination<- ordinate(PS.pig.Norm, method="NMDS", distance="bray", 
                           p.adjust.methods= "bonferroni", permutations = 999)

nmds.scores<- as.data.frame(vegan::scores(nmds.ordination))
nmds.scores<- cbind(nmds.scores, tmp)

genus.scores<- as.data.frame(vegan::scores(nmds.ordination, "species" ))
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
##Compartments infected and Ascaris
##Subset just the infected pigs 
tmp<- row.names(PS.PA.Norm@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  dplyr::filter(InfectionStatus!= "Non_infected")%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.InfAsc<- subset_samples(PS.PA.Norm, Replicate%in%Inf.Keep)

bray_dist<- phyloseq::distance(PS.InfAsc, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.InfAsc,
                      method="PCoA", distance="bray")

tmp<- row.names(PS.InfAsc@sam_data)

tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

pig.ascaris.adonis<- vegan::adonis(bray_dist~ Compartment + AnimalSpecies + Origin,
                                   permutations = 999, data = tmp, na.action = F, strata = tmp$System)

##Store data
foo<- as.data.frame(pig.ascaris.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Pig_Ascaris.csv")

####
## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, tmp$Compartment, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

anova(mvd)
TukeyHSD(mvd)

##Extract centroids, distances and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))
distances<-as.data.frame(data.frame(mvd$distances))

##Select Axis 1 and 2 
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("Compartment","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

##Add sample data
tmp%>%
  dplyr::select(!c(Compartment))%>%
  cbind(seg.data)-> seg.data

ggplot() + 
  geom_point(data=centroids[,1:4], aes(x=PCoA1,y=PCoA2, color= grps, group=grps), size=4, shape= 4) +
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Compartment, shape= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24,21), labels= c("Infected Pig",  "Ascaris"))+
  scale_fill_manual(values = pal.compartment)+
  scale_color_manual(values = pal.compartment)+
  labs(tag= "B)", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Compartment), linetype = 2)+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> B

##Extract distances to centroid for each group
distances%>%
  cbind(tmp)-> distances

distances%>%
  wilcox_test(mvd.distances ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Centroid_Distances_Inf_Ascaris.csv")

##Boxplot 
distances%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x= Compartment, y= mvd.distances))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(shape= 21, size=3, aes(fill= System), color= "black")+
  scale_fill_manual(values = pal.system)+
  xlab("GI compartment")+
  ylab("Distance to centroid")+
  labs(tag= "C)", caption = get_pwc_label(stats.test), fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> C

##Save them individually
ggsave(file = "Figures/Q1_Composition_Infected_Non_Infected.png", plot = A, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Composition_Infected_Compartment.png", plot = B, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Distances_Infected_Compartments.png", plot = C, width = 10, height = 8, dpi = 450)

##Al together
Plot1<- grid.arrange(A,B,C)

ggsave(file = "Figures/Q1_Beta_diversity_Infected_Non_Infected.pdf", plot = Plot1, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Q1_Beta_diversity_Infected_Non_Infected.png", plot = Plot1, width = 10, height = 12, dpi = 450)

rm(A,B,C, Plot1)

##Are the worms microbiomes closer to their host microbiome? 
##Check first at site of infection
##Subset distances for jejunum infected pigs and their worms
tmp%>%
  dplyr::filter(InfectionStatus!= "Non_infected")%>%
  dplyr::filter(Compartment%in% c("Ascaris", "Jejunum"))%>%
  dplyr::filter(System!= "Pig14")%>% #No Jejunum 
  dplyr::filter(System!= "Pig5")%>% #Just one ascaris
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

jejunum.ascaris.adonis<- vegan::adonis(bray_dist~ AnimalSpecies + Origin,
                                   permutations = 999, data = tmp, na.action = F, strata = tmp$System)
##Store data
foo<- as.data.frame(jejunum.ascaris.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Jejunum_Ascaris.csv")

###Higher amount of variance is explain by experiment of origin 
mvd<- vegan::betadisper(bray_dist, tmp$Compartment, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

anova(mvd)

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

hulls<- plyr::ddply(seg.data, "Origin", find_hull)

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= System, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", shape= "Host-Parasite", color= "Origin of samples", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> A
  
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
  scale_shape_manual(values = c(24, 21), labels = c("Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p} {p.signif}")+
  ylab("Bray-Curtis distance (PCo 1)")+
  labs(tag= "D)",  shape = "Host-Parasite", fill= "Individual")+ ##caption = get_pwc_label(stats.test), add in 
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())

# Boxplot PCo2
seg.data%>%
  wilcox_test(v.PCoA2 ~ Compartment)%>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_JejAsc_PCo2.csv") ## Do not plot, ns

###Extract distances within same pig 
### Are worms from each pig closer to the microbiome from its host?
x<- as.matrix(bray_dist)

tmp%>%
  dplyr::filter(Compartment%in% c("Jejunum"))%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.pig.Keep

Inf.pig.Keep<- Inf.pig.Keep$Replicate

tmp%>%
  dplyr::filter(Compartment%in% c("Ascaris"))%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.asc.Keep

Inf.asc.Keep<- Inf.asc.Keep$Replicate

BC.JejAsc<- as.data.frame(x[c(Inf.asc.Keep), c(Inf.pig.Keep)])

BC.JejAsc %>%
  rownames_to_column(var = "Replicate")%>%
  gather("Pig1.Jejunum", "Pig10.Jejunum", "Pig11.Jejunum", "Pig12.Jejunum",
         "Pig13.Jejunum", "Pig2.Jejunum", "Pig3.Jejunum", key = Host, value = BC_dist)%>%
  left_join(tmp, by= "Replicate")-> BC.JejAsc

BC.JejAsc%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Host = fct_relevel(Host, 
                            "Pig1.Jejunum","Pig2.Jejunum","Pig3.Jejunum","Pig4.Jejunum",
                            "Pig5.Jejunum","Pig6.Jejunum","Pig7.Jejunum","Pig8.Jejunum","Pig9.Jejunum",
                            "Pig10.Jejunum","Pig11.Jejunum", "Pig12.Jejunum", "Pig13.Jejunum", "Pig14.Jejunum"))%>%
  dplyr::group_by(Host)%>%
  wilcox_test(BC_dist ~ System)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "System")-> stats.test

BC.JejAsc%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Host = fct_relevel(Host, 
                            "Pig1.Jejunum","Pig2.Jejunum","Pig3.Jejunum","Pig4.Jejunum",
                            "Pig5.Jejunum","Pig6.Jejunum","Pig7.Jejunum","Pig8.Jejunum","Pig9.Jejunum",
                            "Pig10.Jejunum","Pig11.Jejunum", "Pig12.Jejunum", "Pig13.Jejunum", "Pig14.Jejunum"))%>%
  dplyr::group_by(Host)%>%
  ggplot(aes(x= System, y= BC_dist))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(size=3, aes(fill= System, shape= WormSex), color= "black")+
  scale_shape_manual(values = c(23, 22), labels= c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  facet_grid(~Host, scales = "free", space = "free")+
  ylab("Bray-Curtis dissimilarity \n between Ascaris and Jejunum")+
  labs(tag= "B)", fill= "Ascaris origin", shape= "Worm sex")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -1000, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> B

###Enterotype
###Summarize to Genus
PS.JejAsc.Gen<-  tax_glom(PS.JejAsc, "Genus", NArm = F)

##Dominant taxa per sample 
dom.gen.JejAsc<- find.top.asv(PS.JejAsc, "Genus", 1) #--> Genus

dom.gen.JejAsc%>%
  dplyr::select(c(1,3,9))%>%
  dplyr::mutate(Gen.Abund= Abundance)%>%
  dplyr::mutate(Gen.Dom= Genus)%>%
  dplyr::select(c(Replicate, Gen.Abund, Gen.Dom))-> dom.gen.JejAsc 

tmp%>%
  left_join(dom.gen.JejAsc, by="Replicate")-> tmp

dom.phy.JejAsc<- find.top.asv(PS.JejAsc, "Phylum", 1) #--> Phylum

dom.phy.JejAsc%>%
  dplyr::select(c(1,3,5))%>%
  dplyr::mutate(Phy.Abund= Abundance)%>%
  dplyr::mutate(Phy.Dom= Phylum)%>%
  dplyr::select(c(Replicate, Phy.Abund, Phy.Dom))-> dom.phy.JejAsc 

tmp%>%
  left_join(dom.phy.JejAsc, by="Replicate")-> tmp 

tmp%>%
  dplyr::select(c(Replicate, Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom)) -> tmp.Dom

rownames(tmp.Dom)<- tmp.Dom$Replicate 
  
tmp.Dom%>%
  dplyr::select(c(Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom))%>%
  cbind(seg.data)-> seg.data 

###Estimate centroids for enterotypes
###Higher amount of variance is explain by experiment of origin 
mvd<- vegan::betadisper(bray_dist, tmp$Gen.Dom, type = "centroid")
anova(mvd)

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Gen.Dom, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "C)", shape= "Host-Parasite", fill= "Enterotype")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  scale_color_manual(values = tax.palette)+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> C

jejunum.ascaris.dom.adonis<- vegan::adonis(bray_dist~ AnimalSpecies + Origin + Gen.Dom,
                                       permutations = 999, data = tmp, na.action = F, strata = tmp$System)

##Store data
foo<- as.data.frame(jejunum.ascaris.dom.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Jejunum_Ascaris_Enterotype.csv")

##Barplot by sample 
gen.JejAsc<- count.high.genus(x = PS.JejAsc.Gen, num = 1) ##Taxa less represented had less than 1% Relative abundance

tmp%>%
  left_join(gen.JejAsc, by="Replicate")-> gen.JejAsc

#set color palette to accommodate the number of genera
length(unique(gen.JejAsc$Genus))

#plot
gen.JejAsc%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  ggplot(aes(x=Replicate, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack", width=.75) + 
  facet_grid(~System, scales = "free", space = "free")+
  scale_fill_manual(values=tax.palette) + 
  theme_bw()+
  labs(tag= "A)")+
  ylab("Relative abundance (%)")+
  xlab("Sample ID")+
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=4))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))-> H

###Compare abundance of dominants 
gen.JejAsc<- count.high.genus(x = PS.JejAsc.Gen, num = 0) ##Taxa less represented had less than 1% Relative abundance

tmp%>%
  left_join(gen.JejAsc, by="Replicate")-> gen.JejAsc

###Based on PCoA we have 6 Enterotypes:
###Lactobacillus, Escherichia-Shigella, Prevotella, Streptococcus, Clostridium sensu stricto 1 & Romboutsia

##Create variable cluster 
gen.JejAsc<- count.high.genus(x = PS.JejAsc.Gen, num = 0) ##Taxa less represented had less than 1% Relative abundance

tmp%>%
  left_join(gen.JejAsc, by="Replicate")-> gen.JejAsc

gen.JejAsc%>%
  dplyr::filter(Genus%in%c("Clostridium sensu stricto 1", 
                           "Lactobacillus",  "Romboutsia", "Escherichia-Shigella", 
                           "Prevotella", "Streptococcus"))%>%
  dplyr::select(10:20, 22, 28)%>%
  mutate(Gen.Dom = fct_relevel(Gen.Dom,"Clostridium sensu stricto 1", "Lactobacillus",
                               "Escherichia-Shigella", "Prevotella", "Romboutsia", "Streptococcus"))-> Enterotype.abund

Enterotype.abund %>%
  group_by(Genus)%>%
  wilcox_test(Abundance ~ Gen.Dom)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "Gen.Dom")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_JejAsc_Enterotype_abundances.csv")

stats.test%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(y.position= c(100, 110, 90, 100, 60, 70))-> stats.test

# New facet label names for supp variable
#supp.labs <- c("Clostridium s. stricto 1", "Lactobacillus", "Romboutsia")
#names(supp.labs) <- c("Clostridium sensu stricto 1", "Lactobacillus", "Romboutsia")

Enterotype.abund%>%
  group_by(Genus)%>%
  ggplot(aes(x= Gen.Dom, y= Abundance))+
  facet_grid(~Genus, scales = "free", space = "free")+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_jitter(size=3, width = 0.25, aes(fill= System, shape= Compartment), color= "black")+
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  ylab("Relative abundance (%)")+
  labs(tag= "D)", fill= "Individual", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0, hide.ns = T,
                     tip.length = 0)-> D

##Try a real clustering unbiased 
seg.data%>%
  dplyr::select(c("v.PCoA1","v.PCoA2"))-> x

set.seed(2021)
# Determine number of clusters
#wss <- (nrow(bray_dist)-1)*sum(apply(bray_dist,2,var))
#for (i in 2:15) wss[i] <- sum(kmeans(mydata,
#                                     centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters",
#     ylab="Within groups sum of squares") 

###Run clusters
clusters<- kmeans(x, centers = 3, iter.max = 10, nstart = 1)

tmp2<- as.data.frame(clusters$cluster)
colnames(tmp2)<- "Cluster"

seg.data<- cbind(seg.data, tmp2)
seg.data$Cluster<- as.factor(seg.data$Cluster)

clu.mean<- as.data.frame(clusters$centers)
clu.mean%>%
  rownames_to_column("grps")-> clu.mean

ggplot() + 
  geom_point(data=clu.mean, aes(x=v.PCoA1,y=v.PCoA2, color= grps, group=grps), size=4, shape= 4) +
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Cluster, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  labs(tag= "A)", shape= "Host-Parasite", fill= "Cluster")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Cluster), linetype = 2)+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> E

##Check abundance by cluster 
tmp2%>%
  rownames_to_column("Replicate")-> tmp2

Enterotype.abund%>%
  left_join(tmp2, by= "Replicate")%>%
  dplyr::filter(Genus%in%c("Clostridium sensu stricto 1", 
                           "Lactobacillus",  "Romboutsia"))-> Enterotype.abund

Enterotype.abund$Cluster<- as.factor(Enterotype.abund$Cluster)

Enterotype.abund %>%
  group_by(Genus)%>%
  wilcox_test(Abundance ~ Cluster)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "Cluster")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_JejAsc_KClusters_abundances.csv")

Enterotype.abund%>%
  group_by(Genus)%>%
  ggplot(aes(x= Cluster, y= Abundance))+
  facet_grid(~Genus, scales = "free", space = "free")+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_jitter(size=3, width = 0.25, aes(fill= System, shape= Compartment), color= "black")+
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  ylab("Relative abundance (%)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test), fill= "Individual", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  scale_x_discrete(labels=c("Cluster_1" = "1", 
                            "Cluster_2" = "2",
                            "Cluster_3" = "3"))+
  theme(text = element_text(size=16))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0, hide.ns = T,
                     tip.length = 0)-> G
##Save them individually
ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris.png", plot = A, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_Distance.png", plot = B, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Jejunum_Ascaris.png", plot = C, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Jejunum_Ascaris_Abundance.png", plot = D, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Kmean_Jejunum_Ascaris.png", plot = E, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Kmean_Jejunum_Ascaris_Abundance.png", plot = G, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Composition_Jejunum_Ascaris_Abundance.png", plot = H, width = 10, height = 8, dpi = 450)

###Save figures 
##Al together
##Host
B+
  guides(fill = FALSE)-> B

Plot1<- grid.arrange(A,B)

ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_All.png", plot = Plot1, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_All.pdf", plot = Plot1, width = 10, height = 8, dpi = 450)

##Enterotype
D+
  guides(shape = FALSE)-> D

Plot2<- grid.arrange(C,D)

ggsave(file = "Figures/Q1_Enterotype_Jejunum_Ascaris_All.png", plot = Plot1, width = 12, height = 10, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Jejunum_Ascaris_All.pdf", plot = Plot1, width = 12, height = 10, dpi = 450)

##Kmean
E+
  guides(shape = FALSE)-> E

Plot3<- grid.arrange(E,G)

ggsave(file = "Figures/Q1_Kmean_Jejunum_Ascaris_All.png", plot = Plot1, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Q1_Kmean_Jejunum_Ascaris_All.pdf", plot = Plot1, width = 10, height = 12, dpi = 450)

rm(A,B,C,D,E, G, H, Plot1, Plot2, Plot3)

###Let's work with Duodenum and Ascaris 
##Are the worms dominated by clostridium closer to duodenum? 

##Subset just the infected pigs 
tmp<- row.names(PS.PA.Norm@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  dplyr::filter(InfectionStatus!= "Non_infected")%>%
  dplyr::filter(Compartment%in% c("Ascaris", "Duodenum"))%>%
  dplyr::filter(System!= "Pig13")%>% #No duodenum
  dplyr::filter(System!= "Pig4")%>% #No Ascaris
  dplyr::filter(System!= "Pig5")%>% #Just one ascaris
  dplyr::group_by(System)%>%
  dplyr::filter(n()>=2)%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.DuoAsc<- subset_samples(PS.PA.Norm, Replicate%in%Inf.Keep)

bray_dist<- phyloseq::distance(PS.DuoAsc, 
                               method="bray", weighted=F)

ordination<- ordinate(PS.DuoAsc,
                      method="PCoA", distance="bray")

tmp<- row.names(PS.DuoAsc@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

duodenum.ascaris.adonis<- vegan::adonis(bray_dist~ AnimalSpecies + Origin,
                                       permutations = 999, data = tmp, na.action = F, strata = tmp$System)

##Store data
foo<- as.data.frame(duodenum.ascaris.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Duodenum_Ascaris.csv")

###Higher amount of variance is explain by experiment of origin 
mvd<- vegan::betadisper(bray_dist, tmp$Compartment, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

anova(mvd)

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
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Duodenum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", shape= "Host-Parasite", color= "Origin of samples", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> A

###Create Boxplot to compare distances at PCo1 and PCo2 among groups
seg.data%>%
  wilcox_test(v.PCoA1 ~ Compartment)%>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_DuoAsc_PCo1.csv")

seg.data%>%
  ggplot(aes(x=Compartment, y= v.PCoA1))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_point(size=3, aes(shape=InfectionStatus, fill= System)) +
  scale_shape_manual(values = c(24, 21), labels = c("Pig (Duodenum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  theme_bw()+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0, label = "{p} {p.signif}")+
  ylab("Bray-Curtis distance (PCo 1)")+
  labs(tag= "D)",  shape = "Host-Parasite", fill= "Individual")+ ##caption = get_pwc_label(stats.test), add in 
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())

# Boxplot PCo2
seg.data%>%
  wilcox_test(v.PCoA2 ~ Compartment)%>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_DuoAsc_PCo2.csv") ##ns not plotting 

###Distance to host duodenum
###Extract distances within same pig 
### Are worms from each pig closer to the microbiome from its host?

x<- as.matrix(bray_dist)

tmp%>%
  dplyr::filter(Compartment%in% c("Duodenum"))%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.pig.Keep

Inf.pig.Keep<- Inf.pig.Keep$Replicate

tmp%>%
  dplyr::filter(Compartment%in% c("Ascaris"))%>%
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.asc.Keep

Inf.asc.Keep<- Inf.asc.Keep$Replicate

BC.DuoAsc<- as.data.frame(x[c(Inf.asc.Keep), c(Inf.pig.Keep)])

BC.DuoAsc %>%
  rownames_to_column(var = "Replicate")%>%
  gather("Pig1.Duodenum", "Pig10.Duodenum", "Pig11.Duodenum", "Pig12.Duodenum",
         "Pig14.Duodenum", "Pig2.Duodenum", "Pig3.Duodenum", key = Host, value = BC_dist)%>%
  left_join(tmp, by= "Replicate")-> BC.DuoAsc

#Are the ascaris closer to host duodenum

BC.DuoAsc%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Host = fct_relevel(Host, 
                            "Pig1.Duodenum",  "Pig2.Duodenum", "Pig3.Duodenum",
                            "Pig10.Duodenum", "Pig11.Duodenum", "Pig12.Duodenum","Pig14.Duodenum"))%>%
  dplyr::group_by(Host)%>%
  wilcox_test(BC_dist ~ System)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "System")-> stats.test

BC.DuoAsc%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Host = fct_relevel(Host, 
                            "Pig1.Duodenum",  "Pig2.Duodenum", "Pig3.Duodenum",
                            "Pig10.Duodenum", "Pig11.Duodenum", "Pig12.Duodenum","Pig14.Duodenum"))%>%
  dplyr::group_by(Host)%>%
  ggplot(aes(x= System, y= BC_dist))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(size=3, aes(fill= System, shape= WormSex), color= "black")+
  scale_shape_manual(values = c(23, 22), labels= c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  facet_grid(~Host, scales = "free", space = "free")+
  ylab("Bray-Curtis dissimilarity \n between Ascaris and Duodenum")+
  labs(tag= "B)", fill= "Ascaris origin", shape= "Worm sex")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -1000, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> B

###Effect of dominant taxa (Enterotype)
###Summarize to Genus
PS.DuoAsc.Gen<-  tax_glom(PS.DuoAsc, "Genus", NArm = F)

##Dominant taxa per sample 
dom.gen.DuoAsc<- find.top.asv(PS.DuoAsc, "Genus", 1) #--> Genus

dom.gen.DuoAsc%>%
  dplyr::select(c(1,3,9))%>%
  dplyr::mutate(Gen.Abund= Abundance)%>%
  dplyr::mutate(Gen.Dom= Genus)%>%
  dplyr::select(c(Replicate, Gen.Abund, Gen.Dom))-> dom.gen.DuoAsc 

tmp%>%
  left_join(dom.gen.DuoAsc, by="Replicate")-> tmp

dom.phy.DuoAsc<- find.top.asv(PS.DuoAsc, "Phylum", 1) #--> Phylum

dom.phy.DuoAsc%>%
  dplyr::select(c(1,3,5))%>%
  dplyr::mutate(Phy.Abund= Abundance)%>%
  dplyr::mutate(Phy.Dom= Phylum)%>%
  dplyr::select(c(Replicate, Phy.Abund, Phy.Dom))-> dom.phy.DuoAsc 

tmp%>%
  left_join(dom.phy.DuoAsc, by="Replicate")-> tmp 

tmp%>%
  dplyr::select(c(Replicate, Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom)) -> tmp.Dom

rownames(tmp.Dom)<- tmp.Dom$Replicate 

tmp.Dom%>%
  dplyr::select(c(Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom))%>%
  cbind(seg.data)-> seg.data 

###Estimate centroids for enterotypes
###Higher amount of variance is explain by experiment of origin 
mvd<- vegan::betadisper(bray_dist, tmp$Gen.Dom, type = "centroid")
anova(mvd)

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))

ggplot() + 
  #geom_point(data=centroids[c(1,4),1:4], aes(x=PCoA1,y=PCoA2, color= grps, group=grps), size=4, shape= 4) +
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Gen.Dom, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Duodenum)", "Ascaris"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "C)", shape= "Host-Parasite", fill= "Enterotype")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  #stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Gen.Dom), linetype = 2)+
  scale_color_manual(values = tax.palette)+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> C

duodenum.ascaris.dom.adonis<- vegan::adonis(bray_dist~ AnimalSpecies + Origin + Gen.Dom,
                                           permutations = 999, data = tmp, na.action = F, strata = tmp$System)

##Store data
foo<- as.data.frame(duodenum.ascaris.dom.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Duodenum_Ascaris_Dom.csv")

###Compare abundance of dominants 
gen.DuoAsc<- count.high.genus(x = PS.DuoAsc.Gen, num = 0) ##Taxa less represented had less than 1% Relative abundance

tmp%>%
  left_join(gen.DuoAsc, by="Replicate")-> gen.DuoAsc

###Based on PCoA we have 3 clusters:
###Cluster 1: Lactobacillus, Escherichia-Shigella, Prevotella, Streptococcus
###Cluster 2: Clostridium sensu stricto 1
###Cluster 3: Romboutsia, Clostridium sensu stricto 1

##Create variable cluster 
gen.DuoAsc%>%
  dplyr::filter(Genus%in%c("Clostridium sensu stricto 1", 
                           "Lactobacillus",  "Romboutsia", "Escherichia-Shigella", "Helicobacter", 
                           "Prevotella", "Streptococcus"))%>%
  dplyr::select(10:20, 22, 28)%>%
  mutate(Gen.Dom = fct_relevel(Gen.Dom,"Clostridium sensu stricto 1", "Lactobacillus",
                               "Escherichia-Shigella", "Prevotella", "Romboutsia", "Helicobacter"))-> Enterotype.abund

Enterotype.abund %>%
  group_by(Genus)%>%
  wilcox_test(Abundance ~ Gen.Dom)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "Gen.Dom")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_DuoAsc_Enterotype_abundances.csv")

stats.test%>%
  dplyr::filter(p.adj.signif!= "ns")-> stats.test 

#%>%
#dplyr::mutate(y.position= c(100, 110, 90, 100, 60, 70))

# New facet label names for supp variable
supp.labs <- c("Clostridium s. stricto 1","Lactobacillus",
               "Escherichia-Shigella", "Prevotella", "Romboutsia", "Helicobacter")
names(supp.labs) <- c("Clostridium sensu stricto 1", "Lactobacillus",
                      "Escherichia-Shigella", "Prevotella", "Romboutsia", "Helicobacter")

Enterotype.abund%>%
  group_by(Genus)%>%
  ggplot(aes(x= Gen.Dom, y= Abundance))+
  facet_grid(~Genus, scales = "free")+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_jitter(size=3, width = 0.25, aes(fill= System, shape= Compartment), color= "black")+
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Duodenum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  ylab("Relative abundance (%)")+
  labs(tag= "D)", fill= "Individual", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0, hide.ns = T,
                     tip.length = 0)-> D

###K-mean clustering
seg.data%>%
  dplyr::select(c("v.PCoA1","v.PCoA2"))-> x

set.seed(2021)
# Determine number of clusters
#wss <- (nrow(bray_dist)-1)*sum(apply(bray_dist,2,var))
#for (i in 2:15) wss[i] <- sum(kmeans(mydata,
#                                     centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters",
#     ylab="Within groups sum of squares") 

###Run clusters
clusters<- kmeans(x, centers = 3, iter.max = 10, nstart = 1)

tmp2<- as.data.frame(clusters$cluster)
colnames(tmp2)<- "Cluster"

seg.data<- cbind(seg.data, tmp2)
seg.data$Cluster<- as.factor(seg.data$Cluster)

clu.mean<- as.data.frame(clusters$centers)
clu.mean%>%
  rownames_to_column("grps")-> clu.mean

ggplot() + 
  geom_point(data=clu.mean, aes(x=v.PCoA1,y=v.PCoA2, color= grps, group=grps), size=4, shape= 4) +
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Cluster, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Duodenum)", "Ascaris"))+
  labs(tag= "A)", shape= "Host-Parasite", fill= "Cluster")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Cluster), linetype = 2)+
  theme_bw()+
  #geom_text(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, label=Replicate),hjust=0, vjust=0)+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))-> E

##Check abundance by cluster 
Enterotype.abund$Cluster<- NULL

tmp2%>%
  rownames_to_column("Replicate")-> tmp2

Enterotype.abund%>%
  left_join(tmp2, by= "Replicate")-> Enterotype.abund

Enterotype.abund$Cluster<- as.factor(Enterotype.abund$Cluster)

Enterotype.abund %>%
  group_by(Genus)%>%
  wilcox_test(Abundance ~ Cluster)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "Cluster")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_DuoAsc_KClusters_abundances.csv")

Enterotype.abund%>%
  group_by(Genus)%>%
  ggplot(aes(x= Cluster, y= Abundance))+
  facet_grid(~Genus, scales = "free", space = "free")+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_jitter(size=3, width = 0.25, aes(fill= System, shape= Compartment), color= "black")+
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Duodenum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  ylab("Relative abundance (%)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test), fill= "Individual", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  scale_x_discrete(labels=c("Cluster_1" = "1", 
                            "Cluster_2" = "2",
                            "Cluster_3" = "3"))+
  theme(text = element_text(size=16))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0, hide.ns = T,
                     tip.length = 0)-> G


##Save them individually
ggsave(file = "Figures/Q1_Host_Duodenum_Ascaris.png", plot = A, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Duodenum_Ascaris_Distance.png", plot = B, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Duodenum_Ascaris.png", plot = C, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Duodenum_Ascaris_Abundance.png", plot = D, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Kmean_Duodenum_Ascaris.png", plot = E, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Kmean_Duodenum_Ascaris_Abundance.png", plot = G, width = 10, height = 8, dpi = 450)

###Save figures 
##Al together
##Host
B+
  guides(fill = FALSE)-> B

Plot1<- grid.arrange(A,B)

ggsave(file = "Figures/Q1_Host_Duodenum_Ascaris_All.png", plot = Plot1, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Duodenum_Ascaris_All.pdf", plot = Plot1, width = 10, height = 8, dpi = 450)

##Enterotype
Plot2<- grid.arrange(C,D)

ggsave(file = "Figures/Q1_Enterotype_Duodenum_Ascaris_All.png", plot = Plot1, width = 12, height = 10, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Duodenum_Ascaris_All.pdf", plot = Plot1, width = 12, height = 10, dpi = 450)

##Kmean
E+
  guides(shape = FALSE)-> E

Plot3<- grid.arrange(E,G)

ggsave(file = "Figures/Q1_Kmean_Duodenum_Ascaris_All.png", plot = Plot1, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Q1_Kmean_Duodenum_Ascaris_All.pdf", plot = Plot1, width = 10, height = 12, dpi = 450)

rm(A,B,C,D,E, G, Plot1, Plot2, Plot3)
