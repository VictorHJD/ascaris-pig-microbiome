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

##Mixture parameters pi and theta
mixturewt(best)

###Result
#Enterotype pi          theta
#1          0.7155963   26.32667
#2          0.2844037   247.87276

#Sample-component assignments

ass <- apply(mixture(best), 1, which.max)

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


