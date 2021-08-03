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
                "Lactobacillus"=  "#631879FF", "Clostridium sensu stricto 1"= "#00468BFF","Prevotella"= "#3C5488FF",
                "Escherichia-Shigella" = "#0072B5FF", "Roseburia"= "#E762D7FF", "Pseudomonas" = "#ED0000FF",
                "Agathobacter" = "#0099B4FF" , "Bifidobacterium" = "#E64B35FF", "Veillonella"= "#3B4992FF", 
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

find_hull<- function(df) df[chull(df$v.PCoA1, df$v.PCoA2), ]

hulls<- plyr::ddply(seg.data, "Origin", find_hull)

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, shape= InfectionStatus, fill= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F)+
  geom_polygon(data= hulls, aes(x=v.PCoA1,y=v.PCoA2, color= Origin), alpha= 0.1)+
  scale_color_manual(values = c("#4682B4", "#B47846"), labels= c("Experiment 1", "Experiment 2"))+
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

###Higher amount of variance is explain by experiment of origin 

mvd<- vegan::betadisper(bray_dist, tmp$Compartment, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

plot(mvd, ellipse = T, hull = F)
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

hulls<- plyr::ddply(seg.data, "Origin", find_hull)

ggplot() + 
  geom_polygon(data=hulls, aes(x=v.PCoA1,y=v.PCoA2, color= Origin), alpha= 0.1)+
  scale_color_manual(values = c("#4682B4", "#B47846"), labels= c("Experiment 1", "Experiment 2"))+
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= System, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "C)", shape= "Host-Parasite", color= "Origin of samples", fill= "Individual")+
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
  scale_shape_manual(values = c(24, 21), labels = c("Pig (Jejunum)", "Ascaris"))+
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
  scale_shape_manual(values = c(24, 21), labels = c("Pig (Jejunum)", "Ascaris"))+
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

###Summarize to Genus
PS.JejAsc.Gen<-  tax_glom(PS.JejAsc, "Genus", NArm = F)

##Barplot by sample 
gen.JejAsc<- count.high.genus(x = PS.JejAsc.Gen, num = 1)

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
  geom_bar(aes(), stat="identity", position="stack", width=.75, color= "black") + 
  facet_grid(~System, scales = "free", space = "free")+
  scale_fill_manual(values=tax.palette) + 
  theme_bw()+
  labs(tag= "B)")+
  ylab("Relative abundance (%)")+
  xlab("Sample ID")+
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=4))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


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

x[c(Inf.asc.Keep), c(Inf.pig.Keep)]


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
