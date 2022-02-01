####################
##Project: Ascaris - Pig Microbiome
##Aim: DMM Enterotyping
##Author: Morgan Essex mod. by Víctor Hugo Jarquín-Díaz
##Root repo setwd("../Ascaris/ascaris/")
####################

library(tidyverse)
library(magrittr)
library(phyloseq)
library(DirichletMultinomial)
library(microbiome)
library(reshape2)

# read in phyloseq objects
PS.pig.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.pig.Norm.Rds") ## Data just merged pigs normalized for beta diversity plots 
PS.Asc.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Asc.Norm.Rds") ## Data all Ascaris normalized for beta diversity plots 
PS.PA.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Norm.Rds") ## Data merged pigs and Ascaris (not SH) normalized for beta diversity plots 

###For methods 
#Enterotype classifications were performed from the genus abundance matrix using the Dirichlet multinomial mixture (DMM) method as described in Holmes and al. [30] and implemented in the Dirichlet Multinomial R package. 
###For Pig and Ascaris toghether
# count abundances
abun.tbl <- PS.PA@otu_table %>%
  as.data.frame() %>%
  rownames_to_column('ASV') %>%
  as_tibble() %>%
  tidyr::gather('SampleID','Ab',-ASV) %>%
  spread('ASV','Ab') %>%
  column_to_rownames('SampleID')

# ASV:Genus name pairs
asv.names <- PS.PA@tax_table %>%
  as.data.frame() %>%
  rownames_to_column('ASV') %>% as_tibble() %>%
  tidyr::gather('Rank','Value',-ASV) %>%
  dplyr::filter(Rank != 'Species') %>%
  nest(data = -ASV) %>%
  mutate(TaxaID = map_chr(data, function(x) {
    x %>% 
      dplyr::filter(Value != '?') %>%
      use_series(Value) %>% 
      tail(1) })) %>% 
  mutate(TaxaUp = map_chr(data, function(x) {
    x %>% 
      dplyr::filter(Value != '?') %>%
      use_series(Value) %>% 
      tail(2) %>% head(1)
  })) %>%
  mutate(Duped = duplicated(TaxaID)) %>%
  mutate(TaxaID = ifelse(Duped, paste(TaxaUp, TaxaID), TaxaID)) %>%
  dplyr::select(-c(data, TaxaUp, Duped))

###Remove zero rowsand change direction, rows should be samples, columns should be ASVs
abun.tbl.nozero<- t(abun.tbl[apply(abun.tbl[,-1], 1, function(x) !all(x==0)),])
asv.names<- asv.names[asv.names$ASV%in%colnames(abun.tbl.nozero),]

# fit dirichlet multinomial models
future::plan(future::multisession, workers = 6)
 dmns <- c(1:6) %>%
     furrr::future_map(~ dmn(count = as.matrix(abun.tbl.nozero), k = .))
 
 saveRDS(dmns, 'Data/DMM-Enterotypes.rds')
 
 fit <- lapply(1:6, dmn, count = as.matrix(abun.tbl.nozero), verbose=TRUE)
 
 saveRDS(fit, 'Data/DMM-Enterotypes_2.rds')
 
# check optimal number of metacommunities
lplc <- sapply(dmns, laplace)
plot(lplc, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit")

# assign metacommunities
mixtures <- dmns %>%
  map(~ mixture(., assign = TRUE)) 

all.enterotypes <- mixtures %>%
  bind_cols() %>%
  set_names(c('k=1','k=2','k=3','k=4','k=5','k=6')) %>%
  add_column(SampleID = names(mixtures[[1]]), .before = 1)

# model comparison
lplc <- sapply(dmns, laplace)
BIC <- sapply(dmns, BIC)
AIC <- sapply(dmns, AIC)
# optimal number of metacommunities (Laplace information criterion)
dmns[[which.min(lplc)]] 
# optimal number of metacommunities (Bayesian information criterion)
dmns[[which.min(BIC)]] 
# optimal number of metacommunities (Akaike information criterion)
dmns[[which.min(AIC)]] 

##Pick the best model based on Laplace information criterion
best <- dmns[[which.min(unlist(lplc))]]

# show heatmap
heatmapdmn(as.matrix(abun.tbl.nozero), dmns[[1]], dmns[[2]], 10)

# which genera do these correspond to
named.types <- asv.names %>% 
  filter(ASV %in% c('ASV1','ASV11','ASV12','ASV4','ASV24','ASV14','ASV19', 'ASV25', 'ASV5', 'ASV18'))

##Taxa fitted to Dirichlet components 1 and 2
lattice::splom(log(fitted(best)))

#The posterior mean difference between the best and single-component Dirichlet multinomial model 
#measures how each component differs from the population average; the sum is a measure of total difference from the mean.
p0 <- fitted(fit[[1]], scale=TRUE) # scale by theta
p2 <- fitted(best, scale=TRUE)
colnames(p2) <- paste("m", 1:2, sep="")
(meandiff <- colSums(abs(p2 - as.vector(p0))))

# m1        m2 
#0.6568182 0.5140916 

sum(meandiff)
#1.17091

##Mixture parameters pi and theta
mixturewt(best)

###Result
#Enterotype pi          theta
#1          0.7155963   26.32667
#2          0.2844037   247.87276

##Top ten components 
diff <- rowSums(abs(p2 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- head(cbind(Mean=p0[o], p2[o,], diff=diff[o], cdiff), 10)

named.types%>%
  cbind(df)-> df

##Top drivers
driv2 <- dplyr::select(df, c(TaxaID,m2)) %>%
  # Arrange OTUs by assignment strength
  arrange(m2) %>%
  mutate(ASV= paste0(TaxaID, " ", row.names(df)))%>%
  mutate(ASV = factor(ASV, levels = unique(ASV))) %>%
  ggplot(aes(x =ASV, y = m2), color= ASV) +
  geom_bar(stat = "identity") +
  #scale_color_manual(values=c("#E64B35FF", "#E762D7FF", "#00468BFF", "#00A087FF", "#00468BFF",
  #                            "#E64B35FF","#E64B35FF","#7876B1FF", "#999933", "#CD534CFF"))%>%
  coord_flip() +
  labs(title = "Top drivers: community type 2", x= "Taxa identity", y= "Fitted value", tag = "B)")+
  theme_bw()+
  theme(text = element_text(size=16))

driv1 <- dplyr::select(df, c(TaxaID,m1)) %>%
  # Arrange OTUs by assignment strength
  arrange(m1) %>%
  mutate(ASV= paste0(TaxaID, " ", row.names(df)))%>%
  mutate(ASV = factor(ASV, levels = unique(ASV))) %>%
  ggplot(aes(x =ASV, y = m1), color= ASV) +
  geom_bar(stat = "identity") +
  #scale_color_manual(values=c("#E64B35FF", "#E762D7FF", "#00468BFF", "#00A087FF", "#00468BFF",
  #                            "#E64B35FF","#E64B35FF","#7876B1FF", "#999933", "#CD534CFF"))%>%
  coord_flip() +
  labs(title = "Top drivers: community type 1", x= "Taxa identity", y= "Fitted value", tag = "A)")+
  theme_bw()+
  theme(text = element_text(size=16))

#Sample-component assignments
ass <- as.data.frame(apply(mixture(best), 1, which.max))

colnames(ass)<- "Community_type"

##NMDS With enterotype 
tmp<- row.names(PS.PA.Norm@sam_data)
tmp<- alphadiv.PA.rare[rownames(alphadiv.PA.rare)%in%tmp, ]

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

tmp<- alphadiv.PA.rare[rownames(alphadiv.PA.rare)%in%tmp, ]

#ass[rownames(ass)%in%tmp, ] #select samples in tmp

tmp%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))-> tmp

pig.ascaris.adonis<- vegan::adonis(bray_dist~ InfectionStatus + Compartment +  System,
                                   permutations = 999, data = tmp, na.action = F, strata = tmp$Origin)

##Store data
foo<- as.data.frame(pig.ascaris.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Pig_Ascaris.csv")

####
## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, tmp$Compartment, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)
anova(mvd)

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
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))

##ANOSIM
##Compartments
#ANOSIM statistic R: 0.4553
#Significance: 0.001
#permutations = 999
compartment.anosim<- vegan::anosim(bray_dist, tmp$Compartment, permutations = 999, strata =tmp$Origin)
##Conclusion: there is difference between the microbial communities from the different compartments.

##Experiment 1 vs Experiment 2
summary(vegan::anosim(bray_dist, tmp$Origin, permutations = 999, strata =tmp$System))
#ANOSIM statistic R: 0.26
#Significance: 1
#permutations = 999
##Conclusion: there is no difference between the microbial communities of Experiment 1 or Experiment 2 

##Pig vs Ascaris
summary(vegan::anosim(bray_dist, tmp$InfectionStatus, permutations = 999, strata =tmp$Origin))
#ANOSIM statistic R: 0.3783
#Significance: 0.001
#permutations = 999
##Conclusion: there is no difference between the microbial communities of Infected pigs and Ascaris

##Individual
summary(vegan::anosim(bray_dist, tmp$System, permutations = 999, strata =tmp$Origin))
#ANOSIM statistic R: 0.2379 
#Significance: 0.006
#permutations = 999
##Conclusion: there is no difference between the microbial communities of Infected pigs and Ascaris

## Non-metric multidimensional scaling
##With phyloseq
nmds.ordination<- ordinate(PS.InfAsc, method="NMDS", distance="bray", trymax= 75,
                           p.adjust.methods= "bonferroni", permutations = 999)

nmds.scores<- as.data.frame(vegan::scores(nmds.ordination))
nmds.scores<- cbind(nmds.scores, tmp)

genus.scores<- as.data.frame(vegan::scores(nmds.ordination, "species"))
genus.data<- as.data.frame(PS.InfAsc@tax_table)
genus.scores<- cbind(genus.scores, genus.data)
rm(genus.data)

genus.scores%>%
  dplyr::filter(!is.na(NMDS1), !is.na(NMDS2)) %>%
  dplyr::filter(!is.na(Genus))-> genus.scores 
genus.scores %>%
  dplyr::filter(rownames(genus.scores)%in%c("ASV76", "ASV24", "ASV44", "ASV42", "ASV59", "ASV16", 
                                            "ASV29", "ASV28", "ASV35", "ASV50"))-> genus.scores ##Here use main contributors from below

nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= Compartment, shape= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24,21), labels= c("Infected Pig",  "Ascaris"))+
  scale_fill_manual(values = pal.compartment)+
  scale_color_manual(values = pal.compartment)+
  labs(tag= "A)", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  stat_ellipse(aes(color= Compartment), linetype = 2)+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= genus.scores, aes(x = 0, y = 0, xend = (NMDS1)*2.5, yend = (NMDS2)*2.5),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_text_repel(data = genus.scores, aes(x = (NMDS1)*3, y = (NMDS2)*3), label= genus.scores$Genus)+
  annotate("text", x = 1.5, y = 3, label= "ANOSIM (compartment) \n")+
  annotate("text", x = 1.5, y = 2.6, label= "stress= 0.1381")+
  annotate("text", x = 1.5, y = 2.9, label= paste0(label = "R = ", round(compartment.anosim$statistic, digits = 3),
                                                   ", p = ", compartment.anosim$signif), color = "black")

###Which samples belong to cluster 1 and which to cluster 2
ab<- as.data.frame(abun.tbl.nozero)

two.clusters <- ab %>% 
  rownames_to_column('SampleID') %>% 
  left_join(all.enterotypes) %>% 
  dplyr::rename(ET='k=2') %>% 
  dplyr::select(-starts_with('k')) 

two.clusters%>%
  dplyr::filter(ET==1)%>%
  dplyr::select(SampleID)-> Ent1

two.clusters%>%
  dplyr::filter(ET==2)%>%
  dplyr::select(SampleID)-> Ent2

#Most Ascaris samples but one, cluster together with Duodenum, Jejunum and Ileum samples.

#Contribution of each taxonomic group to each component

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("ASV", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange ASVs by assignment strength
    arrange(value) %>%
    mutate(ASV = factor(ASV, levels = unique(ASV))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = ASV, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}

two.clusters%>% 
  as_tibble() %>% 
  nest(data = -ET) %>%
  dplyr::mutate(ET = as.character(ET)) %>%
  left_join( separate(named.types, ASV, c('ASV','ET'), sep = '_') ) %>%
  dplyr::select(-ASV) %>%
  mutate(Prevo = map_dbl(data, ~ mean(.$OTU_1))) %>%
  mutate(Bacter = map_dbl(data, ~ mean(.$OTU_2))) %>%
  mutate(Faecali = map_dbl(data, ~ mean(.$OTU_3))) %>%
  mutate(Blautia = map_dbl(data, ~ mean(.$OTU_6))) %>%
  mutate(Rumino = map_dbl(data, ~ mean(.$OTU_9))) 

three.ranked <- three.clusters %>%
  gather('genus','abundance', -c(ET, data, TaxaID)) %>%
  group_by(TaxaID) %>%
  mutate(abundance = min_rank(-abundance)) %>% 
  ungroup() %>% 
  spread('genus','abundance')


