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

# read in phyloseq objects
PS.pig.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.pig.Norm.Rds") ## Data just merged pigs normalized for beta diversity plots 
PS.Asc.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.Asc.Norm.Rds") ## Data all Ascaris normalized for beta diversity plots 
PS.PA.Norm<- readRDS("/fast/AG_Forslund/Victor/data/Ascaris/PS/PS.PA.Norm.Rds") ## Data merged pigs and Ascaris (not SH) normalized for beta diversity plots 

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

###Remove zero rows
abun.tbl.nozero<- abun.tbl[apply(abun.tbl[,-1], 1, function(x) !all(x==0)),]
asv.names<- asv.names[asv.names$ASV%in%rownames(abun.tbl.nozero),]

# fit dirichlet multinomial models
# future::plan(future::multisession, workers = 6)
 dmns <- c(1:6) %>%
     furrr::future_map(~ dmn(count = as.matrix(abun.tbl.nozero), k = .))
# saveRDS(dmns, here('all-data','proc-data','dmm-enterotypes.rds'))
dmns <- here('all-data','proc-data','dmm-enterotypes-named.rds') %>%
  readRDS()

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

# show heatmap
heatmapdmn(as.matrix(abun.tbl), dmns[[1]], dmns[[3]], 10)
heatmapdmn(as.matrix(abun.tbl), dmns[[1]], dmns[[4]], 10)
# OTU_1, OTU_2, OTU_3 (intuitively...)

# which genera do these correspond to
named.types <- otu.names %>% 
  filter(OTU %in% c('OTU_1','OTU_2','OTU_3','OTU_6','OTU_13','OTU_9','OTU_5'))

# ####PCoA plotting colored by entertotypes
# dist <- vegdist(raw_values,  method = "bray", na.rm = T)
# all.pcoa <- cmdscale(dist, k = (nrow(raw_values)-1), eig=TRUE)
# in_data <- as.data.frame(vegan::scores(all.pcoa, choices = c(1,2)))
# in_data$Sample <- as.factor(raw$Patient)
# in_data$Time <- as.factor(raw$Case)
# in_data$Group <- as.factor(raw$corona.dash)
# in_data$DMN <- as.character(Dirichlet_multinomial_all$`DMM_k=4`)
# 
# pdf("PCoA_DMN_K4.pdf", width = 8, height = 6)
# ggplot(data = in_data, aes(x = Dim1, y = Dim2, color = DMN)) +
#     geom_point(size = 2, alpha = 0.8) + theme_classic() + 
#     ggtitle("Microbiome (bray)") +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     geom_line(aes(x = Dim1, y = Dim2, group = Sample),
#               size = 0.3) +
#     stat_ellipse(aes(color = DMN))
# dev.off()

### Write enterotype text file ###
### Make sure to check heat map to assign enterotypes!!
### Below is just an example, the enterotype ordering are different in each case
# entero_result <- cbind(raw[ , 1:5], Dirichlet_multinomial_all[ , 4])
# colnames(entero_result)[6] <- "DMM_k=4"
# entero_result$Enterotype <- ifelse(entero_result$`DMM_k=4` == "1",
#                                    yes = "Ruminococcus",
#                                    no = ifelse(entero_result$`DMM_k=4` == "2",
#                                                yes = "Prevotella",
#                                                no = ifelse(entero_result$`DMM_k=4` == "3",
#                                                            yes = "Bacteroides1",
#                                                            no = "Bacteroides2")))
# write.table(entero_result, "Entero_result_assignment.txt", sep = "\t", quote = F)


### CHECKS

four.names <- named.types %>%
  mutate(OTU = str_replace_all(OTU, 'OTU_6','OTU_4'))

three.clusters <- abun.tbl %>% 
  rownames_to_column('SampleID') %>% 
  left_join(all.enterotypes) %>% 
  dplyr::rename(ET='k=3') %>% 
  select(-starts_with('k')) %>% 
  as_tibble() %>% 
  nest(data = -ET) %>%
  mutate(ET = as.character(ET)) %>%
  left_join( separate(named.types, OTU, c('OTU','ET'), sep = '_') ) %>%
  select(-OTU) %>%
  mutate(Prevo = map_dbl(data, ~ mean(.$OTU_1))) %>%
  mutate(Bacter = map_dbl(data, ~ mean(.$OTU_2))) %>%
  mutate(Faecali = map_dbl(data, ~ mean(.$OTU_3))) %>%
  mutate(Blautia = map_dbl(data, ~ mean(.$OTU_6))) %>%
  # mutate(Subdoli = map_dbl(data, ~ mean(.$OTU_13))) %>%
  mutate(Rumino = map_dbl(data, ~ mean(.$OTU_9))) 

three.ranked <- three.clusters %>%
  gather('genus','abundance', -c(ET, data, TaxaID)) %>%
  group_by(TaxaID) %>%
  mutate(abundance = min_rank(-abundance)) %>% 
  ungroup() %>% 
  spread('genus','abundance')


four.clusters <- abun.tbl %>% 
  rownames_to_column('SampleID') %>% 
  left_join(all.enterotypes) %>% 
  dplyr::rename(ET='k=4') %>% 
  select(-starts_with('k')) %>% 
  as_tibble() %>% 
  nest(data = -ET) %>%
  mutate(ET = as.character(ET)) %>%
  left_join( separate(four.names, OTU, c('OTU','ET'), sep = '_') ) %>%
  select(-OTU) %>%
  mutate(Prevo = map_dbl(data, ~ mean(.$OTU_1))) %>%
  mutate(Bacter = map_dbl(data, ~ mean(.$OTU_2))) %>%
  mutate(Faecali = map_dbl(data, ~ mean(.$OTU_3))) %>%
  mutate(Blautia = map_dbl(data, ~ mean(.$OTU_6))) %>%
  # mutate(Subdoli = map_dbl(data, ~ mean(.$OTU_13))) %>%
  mutate(Rumino = map_dbl(data, ~ mean(.$OTU_9)))

four.ranked <- four.clusters %>%
  gather('genus','abundance', -c(ET, data, TaxaID)) %>%
  group_by(TaxaID) %>%
  mutate(abundance = min_rank(-abundance)) %>% 
  ungroup() %>% 
  spread('genus','abundance')
