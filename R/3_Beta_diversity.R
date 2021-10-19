##Project: Ascaris - Pig Microbiome
##Aim: Beta diversity
##Author: Víctor Hugo Jarquín-Díaz
##Root repo setwd("/Ascaris/ascaris/")

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
library("ggrepel")

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

find_hull<- function(df) df[chull(df$v.PCoA1, df$v.PCoA2), ]

venn_phyloseq <-
  function(physeq,
           fact,
           min_nb_seq = 0,
           print_values = TRUE) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be an object of class 'phyloseq'")
    }
    
    moda <-
      as.factor(unlist(unclass(physeq@sam_data[, fact])[fact]))
    if (length(moda) != dim(physeq@otu_table)[1]) {
      data_venn <-
        t(apply(physeq@otu_table, 1, function(x)
          by(x, moda, max)))
    } else if (length(moda) != dim(physeq@otu_table)[2]) {
      data_venn <-
        t(apply(t(physeq@otu_table), 1, function(x)
          by(x, moda, max)))
    } else {
      stop("The factor length and the number of samples must be identical")
    }
    combinations <- data_venn > min_nb_seq
    
    e <- new.env(TRUE, emptyenv())
    cn <- colnames(combinations)
    for (i in seq.int(dim(combinations)[1]))
      if (any(combinations[i, ])) {
        ec <- paste(cn[combinations[i, ]], collapse = "&")
        e[[ec]] <- if (is.null(e[[ec]]))
          1L
        else
          (e[[ec]] + 1L)
      }
    
    en <- ls(e, all.names = TRUE)
    weights <- as.numeric(unlist(lapply(en, get, e)))
    combinations <- as.character(en)
    
    table_value <-
      data.frame(combinations = as.character(combinations),
                 weights = as.double(weights))
    
    venn <- venneuler::venneuler(data_venn > min_nb_seq)
    venn_res <-
      data.frame(
        x = venn$centers[, 1],
        y = venn$centers[, 2],
        radius = venn$diameters / 2
      )
    
    nmod <- nrow(venn_res)
    x1 <- list()
    for (i in 1:nmod) {
      x1[[i]] <- grep(rownames(venn_res)[i], table_value$combinations)
    }
    
    for (i in seq_len(nrow(table_value))) {
      table_value$x[i] <-
        mean(venn$centers[, "x"][unlist(lapply(x1,
                                               function(x)
                                                 sum(x %in% i) > 0))])
      table_value$y[i] <-
        mean(venn$centers[, "y"][unlist(lapply(x1,
                                               function(x)
                                                 sum(x %in% i) > 0))])
    }
    
    df <- venn_res
    df$xlab <- df$x + (df$x - mean(df$x))
    df$ylab <- df$y + (df$y - mean(df$y))
    
    circularise <- function(d, n = 360) {
      angle <- seq(-pi, pi, length = n)
      make_circle <- function(x, y, r, modality) {
        data.frame(x = x + r * cos(angle),
                   y = y + r * sin(angle),
                   modality)
      }
      lmat <- mapply(
        make_circle,
        modality = rownames(d),
        x = d[, 1],
        y = d[, 2],
        r = d[, 3],
        SIMPLIFY = FALSE
      )
      do.call(rbind, lmat)
    }
    
    circles <- circularise(df)
    
    p <-
      ggplot() + geom_polygon(data = circles,
                              aes(x, y, group = modality, fill = modality),
                              alpha = 0.5) + theme_void()
    
    if (print_values) {
      g_legend <- function(agplot) {
        tmp <- ggplot_gtable(ggplot_build(agplot))
        leg <-
          which(sapply(tmp$grobs, function(x)
            x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
      }
      legend <- g_legend(p)
      
      
      grid::grid.newpage()
      vp1 <- viewport(
        width = 0.75,
        height = 1,
        x = 0.375,
        y = .5
      )
      vpleg <-
        viewport(
          width = 0.25,
          height = 0.5,
          x = 0.85,
          y = 0.75
        )
      subvp <- viewport(
        width = 0.3,
        height = 0.3,
        x = 0.85,
        y = 0.25
      )
      print(p + theme(legend.position = "none"), vp = vp1)
      upViewport(0)
      pushViewport(vpleg)
      grid.draw(legend)
      #Make the new viewport active and draw
      upViewport(0)
      pushViewport(subvp)
      grid.draw(gridExtra::tableGrob(table_value[, c(1, 2)], rows = NULL))
    }
    else{
      return(p)
    }
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
                "Pseudoscardovia"= "#BB0021FF", "Aeromonas"= "#FFCD00FF")

tax.palette.SH<- c("Taxa less represented" = "black",  "Unassigned"="lightgray", "Streptococcus"= "#925E9FFF", 
                "Lactobacillus"=  "#631879FF", "Clostridium sensu stricto 1"= "#00468BFF","Bifidobacterium" = "#3C5488FF",
                "Roseburia" = "#0072B5FF", "Escherichia-Shigella"= "#E762D7FF", "Pseudomonas" = "#ED0000FF",
                "Agathobacter" = "#0099B4FF" , "Prevotella" = "#E64B35FF", "Veillonella"= "#3B4992FF", 
                "Turicibacter" = "#AD002AFF", "Terrisporobacter"  = "#00A087FF",  "Romboutsia" = "#F39B7FFF", 
                "Prevotellaceae NK3B31 group"= "#8491B4FF", "Megasphaera"= "#CD534CFF", "Anaerovibrio" = "#FAFD7CFF",
                "Intestinibacter"="#6F99ADFF", "Parasutterella" = "#BC3C29FF", "Helicobacter"  = "#E18727FF",
                "Succinivibrio"= "#7876B1FF", "Prevotellaceae UCG-003" = "#FFDC91FF", 
                "Klebsiella"  = "#EE4C97FF",  "Aliterella"= "#42B540FF", "Succinivibrionaceae UCG-001" = "#925E9FFF", 
                "Erysipelotrichaceae UCG-002"= "#008B45FF", "Alloprevotella"= "#4DBBD5FF",  "Olsenella"= "#B09C85FF", 
                "Pseudoscardovia"= "#BB0021FF", "Aeromonas"= "#FFCD00FF", "Dialister"= "#E69F00", "Phascolarctobacterium"="#56B4E9", "Fusobacterium"= "#009E73", 
                "Bacteroides"= "#0073C2FF", "Aquamonas"="#868686FF", "Enterococcus"="#00FF00")

##Color palette for compartment and system ##

pal.compartment <- c("Ascaris"="#1B9E77","Cecum"= "#D95F02","Colon"= "#7570B3",
                     "Duodenum"= "#E7298A","Ileum"= "#66A61E","Jejunum"="#E6AB02")

pal.system <- c("Pig1"= "#A6761D","Pig2"= "#666666","Pig3"= "#A6CEE3","Pig4"= "#1F78B4",
                "Pig5"= "#B2DF8A","Pig6"= "#33A02C","Pig7"= "#FB9A99","Pig8"="#E31A1C","Pig9"= "#FDBF6F",
                "Pig10"= "#FF7F00","Pig11"= "#CAB2D6","Pig12"= "#6A3D9A","Pig13"= "#FFFF99",  "Pig14"= "#3B3B3BFF", "SH" = "#BB0021FF")

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
################# FIGURE 3 ###########################################
#comparison between infected and non infected pigs (all compartments, merged replicates) 
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

compartment.adonis<- vegan::adonis(bray_dist~ Compartment + InfectionStatus + System,
                                permutations = 999, data = tmp, na.action = F, strata = tmp$Origin)

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

##Just to have an overview!!! (Not included in final plots)
ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, shape= InfectionStatus, fill= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F)+
  labs(tag= "A)", fill  = "Infection status", color= "Compartment")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Compartment), linetype = 2)+
  scale_color_manual(values = pal.compartment)+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))

##Compartments
summary(vegan::anosim(bray_dist, tmp$Compartment, permutations = 999, strata =tmp$Origin))
#ANOSIM statistic R: 0.475
#Significance: 0.001
#permutations = 999

compartment.anosim<- vegan::anosim(bray_dist, tmp$Compartment, permutations = 999, strata =tmp$Origin)
##Conclusion: there is difference between the microbial communities from the different compartments.

anosim.results<- as.data.frame(compartment.anosim$class.vec)
anosim.results$Distance.rank<- as.data.frame(compartment.anosim$dis.rank)
colnames(anosim.results)<- c("Class", "Dis.rank")

##Experiment 1 vs Experiment 2
summary(vegan::anosim(bray_dist, tmp$Origin, permutations = 999, strata =tmp$System))
#ANOSIM statistic R: 0.244
#Significance: 1
#permutations = 999
##Conclusion: there is no difference between the microbial communities of Experiment 1 or Experiment 2 

##Infected vs Non infected
summary(vegan::anosim(bray_dist, tmp$InfectionStatus, permutations = 999, strata =tmp$Origin))
#ANOSIM statistic R: 0.066
#Significance: 0.232
#permutations = 999
##Conclusion: there is no difference between the microbial communities of Infected and non infected pigs

##Individual
summary(vegan::anosim(bray_dist, tmp$System, permutations = 999, strata =tmp$Origin))

#ANOSIM statistic R: 0.169
#Significance: 0.231
#permutations = 999

###Compare distances at PCo1 and PCo2 among groups
##PCo1
seg.data%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(v.PCoA1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus") 

##No differences, don't store!

##PCo2
seg.data%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(v.PCoA2 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")

##No differences, don't store!
## Non-metric multidimensional scaling
##With phyloseq
nmds.ordination<- ordinate(PS.pig.Norm, method="NMDS", distance="bray", 
                           p.adjust.methods= "bonferroni", permutations = 999)

nmds.scores<- as.data.frame(vegan::scores(nmds.ordination))
nmds.scores<- cbind(nmds.scores, tmp)

genus.scores<- as.data.frame(vegan::scores(nmds.ordination, "species"))
genus.data<- as.data.frame(PS.pig.Norm@tax_table)
genus.scores<- cbind(genus.scores, genus.data)
rm(genus.data)

genus.scores%>%
  dplyr::filter(!is.na(NMDS1), !is.na(NMDS2)) %>%
  dplyr::filter(!is.na(Genus))-> genus.scores 
genus.scores %>%
  dplyr::filter(rownames(genus.scores)%in%c("ASV53", "ASV76", "ASV56", "ASV12", "ASV128", "ASV16", 
                                            "ASV69", "ASV126", "ASV26", "ASV28"))-> genus.scores

nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= InfectionStatus, shape= InfectionStatus), size=3) +
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F)+
  labs(tag= "A)", fill  = "Infection status", color= "Compartment")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color= Compartment), linetype = 2)+
  scale_color_manual(values = pal.compartment)+
  geom_segment(data= genus.scores, aes(x = 0, y = 0, xend = (NMDS1)*2, yend = (NMDS2)*2),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_text_repel(data = genus.scores, aes(x = (NMDS1)*2.22, y = (NMDS2)*2.22), label= genus.scores$Genus)+
  annotate("text", x = 1.1, y = 0.9, label= "ANOSIM (compartment) \n")+
  annotate("text", x = 1.1, y = 0.85, label= paste0(label = "R = ", round(compartment.anosim$statistic, digits = 3),
                                                    ", p = ", compartment.anosim$signif), color = "black")-> A1
##Transform dataset to determine contributors
PS.pig.clr <- microbiome::transform(PS.pig, "clr") #Centered log ratio transformation
Ord.pig.clr <- phyloseq::ordinate(PS.pig.clr, "RDA") #principal components analysis

#Examine eigenvalues and % prop. variance explained
head(Ord.pig.clr$CA$eig)
sapply(Ord.pig.clr$CA$eig[1:6], function(x) x / sum(Ord.pig.clr$CA$eig))

##ASVs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.pig.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))
ind_cont_PCA1 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 ASVs contribute for the 18.4% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 ASVs contribute for the 19.8% of the variation in PC2

ind_cont_PCA_top.pig <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge taxonomy
x<- as.data.frame(PS.pig.clr@tax_table)
x<- x[rownames(x)%in%c(ind_cont_PCA_top.pig$ASV),]
x%>%
  rownames_to_column("ASV")->x

ind_cont_PCA_top.pig%>%
  distinct(ASV, .keep_all = T)%>%
  plyr::join(x, by="ASV")%>%
  column_to_rownames("ASV")-> ind_cont_PCA_top.pig

##Taxa explaining variability
write.csv(ind_cont_PCA_top.pig, "Tables/Q1_Principal_Taxa_Infected.csv")

x<- ind.coord[rownames(ind.coord)%in%c(rownames(ind_cont_PCA_top.pig)),]
x%>%
  dplyr::filter(rownames(x)%in%c("ASV53", "ASV76", "ASV56", "ASV12", "ASV128", "ASV16", 
                                 "ASV69", "ASV126", "ASV26", "ASV28"))-> x
  
y<- ind_cont_PCA_top.pig[rownames(ind_cont_PCA_top.pig)%in%c(rownames(x)),]

x<- cbind(x, y)

plot_ordination(PS.pig.clr, ordination = Ord.pig.clr)+ 
  geom_point(size=3, aes(fill= InfectionStatus, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 25), labels = c("Infected", "Non infected"))+
  scale_fill_manual(values = c("#ED0000FF", "#008B45FF"), labels = c("Infected", "Non infected"))+
  labs(tag= "A)", fill  = "Infection status")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_ellipse(aes(color= Compartment), linetype = 2)+
  scale_color_manual(values = pal.compartment)+
  geom_segment(data= x, aes(x = 0, y = 0, xend = (PC1)*30, yend = (PC2)*30),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(24, 25))), shape= F, 
         color= F, arrow= F)+
  annotate("text", x = (x$PC1*37), y = (x$PC2*37), label= x$Genus)+
  xlab(paste0("PC 1 [", round(Ord.pig.clr$CA$eig[1] / sum(Ord.pig.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.pig.clr$CA$eig[2] / sum(Ord.pig.clr$CA$eig)*100, digits = 2), "%]"))

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

##Extract distances to centroid for each group
distances%>%
  cbind(tmp)-> distances

distances%>%
  wilcox_test(mvd.distances ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")

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
nmds.ordination<- ordinate(PS.InfAsc, method="NMDS", distance="bray", trymax= 50,
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
  labs(tag= "B)", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  stat_ellipse(aes(color= Compartment), linetype = 2)+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= genus.scores, aes(x = 0, y = 0, xend = (NMDS1)*2.5, yend = (NMDS2)*2.5),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_text_repel(data = genus.scores, aes(x = (NMDS1)*3, y = (NMDS2)*3), label= genus.scores$Genus)+
  annotate("text", x = 1.5, y = 3, label= "ANOSIM (compartment) \n")+
  annotate("text", x = 1.5, y = 2.9, label= paste0(label = "R = ", round(compartment.anosim$statistic, digits = 3),
                                                    ", p = ", compartment.anosim$signif), color = "black")-> B1

##Transform dataset to determine contributors
tmp<- row.names(PS.PA@sam_data)
tmp<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]

tmp%>%
  dplyr::filter(InfectionStatus!= "Non_infected")%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.InfAsc.clr<- subset_samples(PS.PA, Replicate%in%Inf.Keep)
PS.InfAsc.clr <- microbiome::transform(PS.InfAsc.clr, "clr") #Centered log ratio transformation
Ord.InfAsc.clr <- phyloseq::ordinate(PS.InfAsc.clr, "RDA") #principal components analysis

#Examine eigenvalues and % prop. variance explained
head(Ord.InfAsc.clr$CA$eig)
sapply(Ord.InfAsc.clr$CA$eig[1:6], function(x) x / sum(Ord.InfAsc.clr$CA$eig))

##ASVs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.InfAsc.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))
ind_cont_PCA1 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 ASVs contribute for the 10.1% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 ASVs contribute for the 27.0% of the variation in PC2

ind_cont_PCA_top.InfAsc <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge taxonomy
x<- as.data.frame(PS.InfAsc.clr@tax_table)
x<- x[rownames(x)%in%c(ind_cont_PCA_top.InfAsc$ASV),]
x%>%
  rownames_to_column("ASV")->x

ind_cont_PCA_top.InfAsc%>%
  distinct(ASV, .keep_all = T)%>%
  plyr::join(x, by="ASV")%>%
  column_to_rownames("ASV")-> ind_cont_PCA_top.InfAsc

##Taxa explaining variability
write.csv(ind_cont_PCA_top.InfAsc, "Tables/Q1_Principal_Taxa_Infected_InfAsc.csv")

x<- ind.coord[rownames(ind.coord)%in%c(rownames(ind_cont_PCA_top.InfAsc)),]
x%>%
  dplyr::filter(rownames(x)%in%c("ASV76", "ASV24", "ASV44", "ASV42", "ASV59", "ASV16", 
                                 "ASV29", "ASV28", "ASV35", "ASV50"))-> x

y<- ind_cont_PCA_top.InfAsc[rownames(ind_cont_PCA_top.InfAsc)%in%c(rownames(x)),]

x<- cbind(x, y)

plot_ordination(PS.InfAsc.clr, ordination = Ord.InfAsc.clr)+ 
  geom_point(size=3, aes(fill= Compartment, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 21), labels = c("Infected Pig", "Ascaris"))+
  scale_fill_manual(values = pal.compartment)+
  labs(tag= "B)", fill  = "Compartment", shape= "Host-Parasite")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_ellipse(aes(color= Compartment), linetype = 2)+
  scale_color_manual(values = pal.compartment)+
  geom_segment(data= x, aes(x = 0, y = 0, xend = (PC1)*35, yend = (PC2)*35),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), 
         color=F, arrow= F)+
  annotate("text", x = (x$PC1*40), y = (x$PC2*40), label= x$Genus)+
  xlab(paste0("PC 1 [", round(Ord.InfAsc.clr$CA$eig[1] / sum(Ord.InfAsc.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.InfAsc.clr$CA$eig[2] / sum(Ord.InfAsc.clr$CA$eig)*100, digits = 2), "%]"))

##Save them individually
ggsave(file = "Figures/Q1_NMDS_Composition_Infected_Non_Infected.png", plot = A1, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Composition_Infected_Non_Infected.pdf", plot = A1, width = 10, height = 8, dpi = 450)

ggsave(file = "Figures/Q1_NMDS_Composition_Infected_Compartment.png", plot = B1, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Composition_Infected_Compartment.pdf", plot = B1, width = 10, height = 8, dpi = 450)

##Al together
Plot1<- ggarrange(A1, B1, ncol=1, align = "v")

ggsave(file = "Figures/Figure_3.pdf", plot = Plot1, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Figure_3.png", plot = Plot1, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Figure_3.svg", plot = Plot1, width = 10, height = 12, dpi = 450)

rm(A1, B1, Plot1)

################# FIGURE 4 ###########################################
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

jejunum.ascaris.adonis<- vegan::adonis(bray_dist~ AnimalSpecies + System,
                                   permutations = 999, data = tmp, na.action = F, strata = tmp$Origin)

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

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= System, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", shape= "Host-Parasite", color= "Origin of samples", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))
  
###Create Boxplot to compare distances at PCo1 and PCo2 among groups
## PCo1
seg.data%>%
  wilcox_test(v.PCoA1 ~ Compartment)%>%
  add_significance()%>%
  add_xy_position(x = "Compartment")

### PCo2
seg.data%>%
  wilcox_test(v.PCoA2 ~ Compartment)%>%
  add_significance()%>%
  add_xy_position(x = "Compartment")

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
  geom_point(size=3, aes(fill= System), shape= 21, color= "black")+
  #scale_shape_manual(values = c(23, 22), labels= c("Female", "Male"))+ #if we want to include worm sex now 
  scale_fill_manual(values = pal.system)+
  facet_grid(~Host, scales = "free", space = "free")+
  ylab("Bray-Curtis dissimilarity \n between Ascaris and Jejunum")+
  labs(tag= "B)", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -1000, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> B2

#To test if there is a statistical difference between the microbial communities of two or more groups of samples.
#Null Hypothesis: there is no difference between the microbial communities of your groups of samples.

##Jejunum vs Ascaris
jejunum.ascaris.anosim<- vegan::anosim(bray_dist, tmp$AnimalSpecies, permutations = 999, strata = tmp$Origin)
#ANOSIM statistic R: 0.4932 
#Significance: 0.001
#permutations = 999
##Conclusion: there is difference between the microbial communities of Ascaris microbiome or Jejunum microbiomes 

##System difference
system.anosim<- vegan::anosim(bray_dist, tmp$System, permutations = 999, strata = tmp$Origin)
#ANOSIM statistic R: 0.3927 
#Significance: 0.001
#permutations = 999
##Conclusion: there is difference between the microbial communities of Ascaris microbiome or Jejunum microbiomes 

##Experiment 1 vs Experiment 2
experiment.anosim<- vegan::anosim(bray_dist, tmp$Origin, permutations = 999, strata =tmp$System)
#ANOSIM statistic R: 0.3952 
#Significance: 0.1
#permutations = 999
##Conclusion: there is no difference between the microbial communities of Experiment 1 or Experiment 2 

## Non-metric multidimensional scaling
##With phyloseq
nmds.ordination<- ordinate(PS.JejAsc, method="NMDS", distance="bray", 
                           p.adjust.methods= "bonferroni", permutations = 999)

nmds.scores<- as.data.frame(vegan::scores(nmds.ordination))
nmds.scores<- cbind(nmds.scores, tmp)

genus.scores<- as.data.frame(vegan::scores(nmds.ordination, "species"))
genus.data<- as.data.frame(PS.JejAsc@tax_table)
genus.scores<- cbind(genus.scores, genus.data)
rm(genus.data)

genus.scores%>%
  dplyr::filter(!is.na(NMDS1), !is.na(NMDS2)) %>%
  dplyr::filter(!is.na(Genus))-> genus.scores 
genus.scores %>%
  dplyr::filter(rownames(genus.scores)%in%c("ASV15", "ASV203", "ASV119", "ASV197", "ASV350", "ASV226", 
                                            "ASV156", "ASV171", "ASV84", "ASV118"))-> genus.scores

nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= System, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", shape= "Host-Parasite", color= "Origin of samples", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= genus.scores, aes(x = 0, y = 0, xend = (NMDS1)*2, yend = (NMDS2)*2),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), 
         color=F, arrow= F)+
  geom_text_repel(data = genus.scores, aes(x = (NMDS1)*2.4, y = (NMDS2)*2.4), label= genus.scores$Genus)+
  annotate("text", x = 2.1, y = 1.5, label= "ANOSIM (Host-Parasite) \n")+
  annotate("text", x = 2.1, y = 1.45, label= paste0(label = "R = ", round(jejunum.ascaris.anosim$statistic, digits = 3),
                                                    ", p = ", jejunum.ascaris.anosim$signif), color = "black")-> A2

##Transform dataset to determine contributors
PS.JejAsc.clr<- subset_samples(PS.PA, Replicate%in%Inf.Keep)
PS.JejAsc.clr <- microbiome::transform(PS.JejAsc.clr, "clr") #Centered log ratio transformation
Ord.JejAsc.clr <- phyloseq::ordinate(PS.JejAsc.clr, "RDA") #principal components analysis

#Examine eigenvalues and % prop. variance explained
head(Ord.JejAsc.clr$CA$eig)
sapply(Ord.JejAsc.clr$CA$eig[1:6], function(x) x / sum(Ord.JejAsc.clr$CA$eig))

##ASVs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.JejAsc.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))
ind_cont_PCA1 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 ASVs contribute for the 27.9% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 ASVs contribute for the 20.4% of the variation in PC2

ind_cont_PCA_top.JejAsc <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge taxonomy
x<- as.data.frame(PS.JejAsc.clr@tax_table)
x<- x[rownames(x)%in%c(ind_cont_PCA_top.JejAsc$ASV),]
x%>%
  rownames_to_column("ASV")->x

ind_cont_PCA_top.JejAsc%>%
  distinct(ASV, .keep_all = T)%>%
  plyr::join(x, by="ASV")%>%
  column_to_rownames("ASV")-> ind_cont_PCA_top.JejAsc

##Taxa explaining variability
write.csv(ind_cont_PCA_top.JejAsc, "Tables/Q1_Principal_Taxa_Infected_JejAsc.csv")

x<- ind.coord[rownames(ind.coord)%in%c(rownames(ind_cont_PCA_top.JejAsc)),]
x%>%
  dplyr::filter(rownames(x)%in%c("ASV28", "ASV16", "ASV36", "ASV50", "ASV29", "ASV83", 
                                 "ASV46", "ASV24", "ASV32", "ASV18"))-> x

y<- ind_cont_PCA_top.JejAsc[rownames(ind_cont_PCA_top.JejAsc)%in%c(rownames(x)),]

x<- cbind(x, y)

PS.JejAsc.clr@sam_data$System<- fct_relevel(PS.JejAsc.clr@sam_data$System, 
                            "Pig1","Pig2","Pig3",
                            "Pig10","Pig11", "Pig12", "Pig13")
require(ggrepel)
plot_ordination(PS.JejAsc.clr, ordination = Ord.JejAsc.clr)+ 
  geom_point(size=3, aes(fill= System, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 21), labels = c("Infected Pig", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", fill  = "Individual", shape= "Host-Parasite")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= x, aes(x = 0, y = 0, xend = (PC1)*35, yend = (PC2)*35),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), 
         color=F, arrow= F)+
  geom_text_repel(data = x, aes(x = (PC1)*40, y = (PC2)*40), label= x$Genus)+
  xlab(paste0("PC 1 [", round(Ord.JejAsc.clr$CA$eig[1] / sum(Ord.JejAsc.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.JejAsc.clr$CA$eig[2] / sum(Ord.JejAsc.clr$CA$eig)*100, digits = 2), "%]"))

##Plot abundance of these ASVs
wh0 <- genefilter_sample(PS.JejAsc, filterfun_sample(function(x) x > 5), A=0.01*nsamples(PS.JejAsc))
PS.subset<- prune_taxa(wh0, PS.JejAsc)

##Changes by compartment
phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                    "Jejunum", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3",
                              "Pig10","Pig11", "Pig12", "Pig13"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::filter(!(OTU%in%c("ASV154", "ASV227")))%>% ##Eliminate not bacterial ASVs (plastid derived)
  dplyr::mutate(y.position= c(0.4, 0.12, 0.025, 0.04, 0.35, 0.45))-> stats.test###None of the "driving" ASVs is significantly higher in Jejunum compared to Ascaris

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_ASV_JejAsc.csv")

PS.subset <- subset_taxa(PS.subset, rownames(tax_table(PS.subset)) %in% c(stats.test$OTU))

phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3",
                              "Pig10","Pig11", "Pig12", "Pig13"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  ggplot(data = ., aes(x = Compartment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = System, shape= InfectionStatus), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(24,21), labels = c("Infected Pig", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Relative Abundance (%)", shape = "Infection status") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)+
  facet_wrap(~ OTU, scales = "free")-> Sup3

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

###Merge with NMDS
tmp.Dom%>%
  dplyr::select(c(Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom))%>%
  cbind(nmds.scores)-> nmds.scores

## Non-metric multidimensional scaling with Enterotype
##With phyloseq
##Jejunum vs Ascaris
jejunum.ascaris.anosim<- vegan::anosim(bray_dist, tmp$Gen.Dom, permutations = 999, strata =tmp$System)
#ANOSIM statistic R: 0.4932 
#Significance: 0.001
#permutations = 999
nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= Gen.Dom, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "C)", shape= "Host-Parasite", fill= "Enterotype")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotate("text", x = 0.8, y = 1.6, label= "ANOSIM (Enterotype) \n")+
  annotate("text", x = 0.8, y = 1.5, label= paste0(label = "R = ", round(jejunum.ascaris.anosim$statistic, digits = 3),
                                                    ", p = ", jejunum.ascaris.anosim$signif), color = "black")-> C2
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
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))

jejunum.ascaris.dom.adonis<- vegan::adonis(bray_dist~ AnimalSpecies + System + Gen.Dom,
                                       permutations = 999, data = tmp, na.action = F, strata = tmp$Origin)

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
                              "Pig1","Pig2","Pig3",
                              "Pig10","Pig11", "Pig12", "Pig13"))%>%
  ggplot(aes(x=Replicate, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack", width=.75) + 
  facet_grid(~System, scales = "free", space = "free")+
  scale_fill_manual(values=tax.palette) + 
  theme_bw()+
  labs(tag= "A)")+
  ylab("Relative abundance (%)")+
  xlab("Sample ID")+
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=5))+
  theme(text = element_text(size=16), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())-> D2

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

Enterotype.abund%>%
  group_by(Genus)%>%
  ggplot(aes(x= Gen.Dom, y= Abundance))+
  facet_grid(~Genus, scales = "free", space = "free")+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_jitter(size=3, width = 0.25, aes(fill= System, shape= Compartment), color= "black")+
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  ylab("Relative abundance (%)")+
  labs(tag= "B)", fill= "Individual", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0, hide.ns = T,
                     tip.length = 0)-> E2

##Save them individually
ggsave(file = "Figures/Q1_NMDS_Host_Jejunum_Ascaris.png", plot = A2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Host_Jejunum_Ascaris.pdf", plot = A2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_Distance.png", plot = B2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_Distance.pdf", plot = B2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_ASVs.png", plot = Sup3, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_ASVs.pdf", plot = Sup3, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Enterotype_Jejunum_Ascaris.png", plot = C2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Enterotype_Jejunum_Ascaris.pdf", plot = C2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Composition_Jejunum_Ascaris_Abundance.png", plot = D2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Composition_Jejunum_Ascaris_Abundance.pdf", plot = D2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Jejunum_Ascaris_Abundance.png", plot = E2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Jejunum_Ascaris_Abundance.pdf", plot = E2, width = 10, height = 8, dpi = 450)

###Save figure 4 
##Al together
B2+
  guides(fill = FALSE)-> B2

Plot1<- ggarrange(A2, B2, C2, ncol=1, align = "v")

ggsave(file = "Figures/Figure_4.png", plot = Plot1, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Figure_4.pdf", plot = Plot1, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Figure_4.svg", plot = Plot1, width = 10, height = 12, dpi = 450)

## Supplementary composition

Plot2<- ggarrange(D2, E2, ncol=1)

ggsave(file = "Figures/Supplementary_Figure_4.png", plot = Plot2, width = 15, height = 12, dpi = 450)
ggsave(file = "Figures/Supplementary_Figure_4.pdf", plot = Plot2, width = 15, height = 12, dpi = 450)
ggsave(file = "Figures/Supplementary_Figure_4.svg", plot = Plot2, width = 15, height = 12, dpi = 450)


rm(A2,B2,C2,D2,E2, Plot1, Plot2)

###############################Let's work with Duodenum and Ascaris ####################################
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

duodenum.ascaris.adonis<- vegan::adonis(bray_dist~ AnimalSpecies + System,
                                       permutations = 999, data = tmp, na.action = F, strata = tmp$Origin)

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
  add_xy_position(x = "Compartment") ##ns not plotting 

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
  geom_point(size=3, aes(fill= System), shape=21, color= "black")+
  #scale_shape_manual(values = c(23, 22), labels= c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  facet_grid(~Host, scales = "free", space = "free")+
  ylab("Bray-Curtis dissimilarity \n between Ascaris and Duodenum")+
  labs(tag= "B)", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -1000, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> B

#To test if there is a statistical difference between the microbial communities of two or more groups of samples.
#Null Hypothesis: there is no difference between the microbial communities of your groups of samples.

##Duodenum vs Ascaris
duodenum.ascaris.anosim<- vegan::anosim(bray_dist, tmp$AnimalSpecies, permutations = 999, strata =tmp$Origin)

#ANOSIM statistic R: 0.4732 
#Significance: 0.001
#permutations = 999

##Conclusion: there is difference between the microbial communities of Ascaris and duodenum

##Experiment 1 vs Experiment 2
summary(vegan::anosim(bray_dist, tmp$Origin, permutations = 999, strata =tmp$System))

#ANOSIM statistic R: 0.4515
#Significance: 0.1
#permutations = 999

##Conclusion: there is no difference between the microbial communities of Experiment 1 or Experiment 2 

##Individuals
summary(vegan::anosim(bray_dist, tmp$System, permutations = 999, strata =tmp$Origin))

#ANOSIM statistic R: 0.3651 
#Significance: 0.002 
#permutations = 999

##Conclusion: there is difference between the microbial communities of duodenum and ascaris in different individuals

## Non-metric multidimensional scaling
##With phyloseq
nmds.ordination<- ordinate(PS.DuoAsc, method="NMDS", distance="bray", 
                           p.adjust.methods= "bonferroni", permutations = 999)

nmds.scores<- as.data.frame(vegan::scores(nmds.ordination))
nmds.scores<- cbind(nmds.scores, tmp)

genus.scores<- as.data.frame(vegan::scores(nmds.ordination, "species"))
genus.data<- as.data.frame(PS.DuoAsc@tax_table)
genus.scores<- cbind(genus.scores, genus.data)
rm(genus.data)

genus.scores%>%
  dplyr::filter(!is.na(NMDS1), !is.na(NMDS2)) %>%
  dplyr::filter(!is.na(Genus))-> genus.scores 

genus.scores %>%
  dplyr::filter(rownames(genus.scores)%in%c("ASV36", "ASV50", "ASV29", "ASV28", "ASV26", "ASV64", 
                                            "ASV13", "ASV16", "ASV24", "ASV47", "ASV15", "ASV21"))-> genus.scores

nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= System, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Duodenum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", shape= "Host-Parasite", color= "Origin of samples", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= genus.scores, aes(x = 0, y = 0, xend = (NMDS1)*2, yend = (NMDS2)*2),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), 
         color=F, arrow= F)+
  geom_text_repel(data = genus.scores, aes(x = (NMDS1)*2.2, y = (NMDS2)*2.2), label= genus.scores$Genus)+
  annotate("text", x = 2.1, y = 1, label= "ANOSIM (Host-Parasite) \n")+
  annotate("text", x = 2.1, y = 0.9, label= paste0(label = "R = ", round(duodenum.ascaris.anosim$statistic, digits = 3),
                                                    ", p = ", duodenum.ascaris.anosim$signif), color = "black")-> A3

##Transform dataset to determine contributors
PS.DuoAsc.clr<- subset_samples(PS.PA, Replicate%in%Inf.Keep)

PS.DuoAsc.clr <- microbiome::transform(PS.DuoAsc.clr, "clr") #Centered log ratio transformation
Ord.DuoAsc.clr <- phyloseq::ordinate(PS.DuoAsc.clr, "RDA") #principal components analysis

#Examine eigenvalues and % prop. variance explained
head(Ord.DuoAsc.clr$CA$eig)
sapply(Ord.DuoAsc.clr$CA$eig[1:6], function(x) x / sum(Ord.DuoAsc.clr$CA$eig))

##ASVs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.DuoAsc.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))
ind_cont_PCA1 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 ASVs contribute for the 16.2% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 ASVs contribute for the 14.6% of the variation in PC2

ind_cont_PCA_top.DuoAsc <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge taxonomy
x<- as.data.frame(PS.DuoAsc.clr@tax_table)
x<- x[rownames(x)%in%c(ind_cont_PCA_top.DuoAsc$ASV),]
x%>%
  rownames_to_column("ASV")->x

ind_cont_PCA_top.DuoAsc%>%
  distinct(ASV, .keep_all = T)%>%
  plyr::join(x, by="ASV")%>%
  column_to_rownames("ASV")-> ind_cont_PCA_top.DuoAsc

##Taxa explaining variability
write.csv(ind_cont_PCA_top.DuoAsc, "Tables/Q1_Principal_Taxa_Infected_DuoAsc.csv")

x<- ind.coord[rownames(ind.coord)%in%c(rownames(ind_cont_PCA_top.DuoAsc)),]
x%>%
  dplyr::filter(rownames(x)%in%c("ASV36", "ASV50", "ASV29", "ASV28", "ASV26", "ASV64", 
                                 "ASV13", "ASV16", "ASV24", "ASV47", "ASV15", "ASV21"))-> x

y<- ind_cont_PCA_top.DuoAsc[rownames(ind_cont_PCA_top.DuoAsc)%in%c(rownames(x)),]

x<- cbind(x, y)

PS.DuoAsc.clr@sam_data$System<- fct_relevel(PS.DuoAsc.clr@sam_data$System, 
                                            "Pig1","Pig2","Pig3","Pig4",
                                            "Pig5","Pig6","Pig7","Pig8","Pig9",
                                            "Pig10","Pig11", "Pig12", "Pig13", "Pig14")
require(ggrepel)
plot_ordination(PS.DuoAsc.clr, ordination = Ord.DuoAsc.clr)+ 
  geom_point(size=3, aes(fill= System, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 21), labels = c("Infected Pig (Duodenum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", fill  = "Individual", shape= "Host-Parasite")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= x, aes(x = 0, y = 0, xend = (PC1)*10, yend = (PC2)*10),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), 
         color=F, arrow= F)+
  geom_text_repel(data = x, aes(x = (PC1*12), y = (PC2*12)), label= x$Genus)+
  xlab(paste0("PC 1 [", round(Ord.DuoAsc.clr$CA$eig[1] / sum(Ord.DuoAsc.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.DuoAsc.clr$CA$eig[2] / sum(Ord.DuoAsc.clr$CA$eig)*100, digits = 2), "%]"))-> A2

##Plot abundance of these ASVs
wh0 <- genefilter_sample(PS.DuoAsc, filterfun_sample(function(x) x > 5), A=0.01*nsamples(PS.DuoAsc))
PS.subset<- prune_taxa(wh0, PS.DuoAsc)

##Changes by compartment
phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Ascaris"))%>%
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
  dplyr::mutate(y.position= c(0.25, 5, 2.5, 1.0, 0.35, 1.75, 2.5, 0.17,0.12,0.2,0.075,0.12,0.075,0.10,0.6,1.2,1.3))-> stats.test###None of the "driving" ASVs is significantly higher in Jejunum compared to Ascaris

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_ASV_DuoAsc.csv")

PS.subset <- subset_taxa(PS.subset, rownames(tax_table(PS.subset)) %in% c(stats.test$OTU))

phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  ggplot(data = ., aes(x = Compartment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = System, shape= InfectionStatus), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(24,21), labels = c("Infected Pig", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Relative Abundance (%)", shape = "Infection status") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)+
  facet_wrap(~ OTU, scales = "free")-> B2

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

tmp.Dom%>%
  dplyr::select(c(Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom))%>%
  cbind(nmds.scores)-> nmds.scores

#### Non-metric multidimensional scaling with Enterotype
##With phyloseq
##Duodenum vs Ascaris
duodenum.ascaris.anosim<- vegan::anosim(bray_dist, tmp$Gen.Dom, permutations = 999, strata =tmp$System)
#ANOSIM statistic R: 0.789
#Significance: 0.001
#permutations = 999
nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= Gen.Dom, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Duodenum)", "Ascaris"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "C)", shape= "Host-Parasite", fill= "Enterotype")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotate("text", x = 0.8, y = 1.0, label= "ANOSIM (Enterotype) \n")+
  annotate("text", x = 0.8, y = 0.9, label= paste0(label = "R = ", round(duodenum.ascaris.anosim$statistic, digits = 3),
                                                   ", p = ", duodenum.ascaris.anosim$signif), color = "black")-> C3

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
#supp.labs <- c("Clostridium s. stricto 1","Lactobacillus",
#               "Escherichia-Shigella", "Prevotella", "Romboutsia", "Helicobacter")
#names(supp.labs) <- c("Clostridium sensu stricto 1", "Lactobacillus",
#                      "Escherichia-Shigella", "Prevotella", "Romboutsia", "Helicobacter")

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

##Save them individually
ggsave(file = "Figures/Q1_Host_Duodenum_Ascaris.png", plot = A, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_PCA_Host_Duodenum_Ascaris.png", plot = A2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Host_Duodenum_Ascaris.png", plot = A3, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Duodenum_Ascaris_Distance.png", plot = B, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Duodenum_Ascaris_ASVs.png", plot = B2, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Duodenum_Ascaris.png", plot = C, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Enterotype_Duodenum_Ascaris.png", plot = C3, width = 10, height = 8, dpi = 450)
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
C+
  guides(shape = FALSE)-> C

Plot2<- grid.arrange(C,D)

ggsave(file = "Figures/Q1_Enterotype_Duodenum_Ascaris_All.png", plot = Plot2, width = 12, height = 10, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Duodenum_Ascaris_All.pdf", plot = Plot2, width = 12, height = 10, dpi = 450)

##Kmean
E+
  guides(shape = FALSE)-> E

Plot3<- grid.arrange(E,G)

ggsave(file = "Figures/Q1_Kmean_Duodenum_Ascaris_All.png", plot = Plot3, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Q1_Kmean_Duodenum_Ascaris_All.pdf", plot = Plot3, width = 10, height = 12, dpi = 450)

rm(A,B,C,D,E, G, Plot1, Plot2, Plot3)

################################Figure 5####################################
####Comparisons between Ascaris microbiomes 
##Subset just the same worms from previous comparisons
tmp<- row.names(PS.Asc.Norm@sam_data)
tmp<- alphadiv.Asc[rownames(alphadiv.Asc)%in%tmp, ]

tmp%>%
  dplyr::filter(System!= "Pig14")%>% #No jejunum
  dplyr::filter(System!= "Pig5")%>% #Just one ascaris
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.Asc.Norm<- subset_samples(PS.Asc.Norm, Replicate%in%Inf.Keep)

##Within infected pigs from infected experiment vs Slaughterhouse pigs  
##Location
bray_dist<- phyloseq::distance(PS.Asc.Norm, 
                               method="bray", weighted=F)
ordination<- ordinate(PS.Asc.Norm,
                      method="PCoA", distance="bray")

tmp<- row.names(PS.Asc.Norm@sam_data)
tmp<- alphadiv.Asc[rownames(alphadiv.Asc)%in%tmp, ]

tmp%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                              "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(Location = fct_relevel(Location, 
                                   "FU", "SH"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14", "SH"))-> tmp

worm.adonis<- vegan::adonis(bray_dist~ Location + System + WormSex,
                                   permutations = 999, data = tmp, na.action = F, strata = tmp$Origin)

##Store the result
foo<- as.data.frame(worm.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Worm_all.csv")

####
## Calculate multivariate dispersion (aka distance to the centroid)
mvd<- vegan::betadisper(bray_dist, tmp$Location, type = "centroid")
mvd.perm<- vegan::permutest(mvd, permutations = 999)

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))
vectors<-data.frame(group=mvd$group,data.frame(mvd$vectors))

##Select Axis 1 and 2 
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("Location","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

##Add sample data
tmp%>%
  dplyr::select(!c(Location))%>%
  cbind(seg.data)-> seg.data

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, shape= WormSex, fill= Location), size=3) +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  scale_fill_manual(values = c("#84BD00FF", "#BB0021FF"), labels = c("Local housing", "Slaughterhouse"))+
  guides(fill = guide_legend(override.aes=list(shape=21)))+
  labs(tag= "A)", fill  = "Location of host", shape= "Worm Sex")+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, shape= WormSex, fill= System), size=3) +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  guides(fill = guide_legend(override.aes=list(shape=21)))+
  labs(tag= "B)", fill  = "Host", shape= "Worm Sex")+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))

##Extract distances to centroid for each group
distances<-as.data.frame(data.frame(mvd$distances))

distances%>%
  cbind(tmp)-> distances

distances%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                                               "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(Location = fct_relevel(Location, 
                                       "FU", "SH"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14", "SH"))%>%
  wilcox_test(mvd.distances ~ System)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "System")##NS among systems

##Sex difference within location
distances%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                                     "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(Location = fct_relevel(Location, 
                                       "FU", "SH"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14", "SH"))%>%
  group_by(Location)%>%
  wilcox_test(mvd.distances ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Centroid_Distances_Ascaris_Sex_Location.csv")

##Sex difference within origin
distances%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                                     "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(Location = fct_relevel(Location, 
                                       "FU", "SH"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14", "SH"))%>%
  group_by(Origin)%>%
  wilcox_test(mvd.distances ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Centroid_Distances_Ascaris_Sex_Origin.csv")

##Sex difference within system
distances%>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                                     "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(Location = fct_relevel(Location, 
                                       "FU", "SH"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig5","Pig10","Pig11", 
                                     "Pig12", "Pig13", "Pig14", "SH"))%>%
  dplyr::filter(System!= "Pig5")%>% #Just one ascaris
  group_by(System)%>%
  wilcox_test(mvd.distances ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Centroid_Distances_Ascaris_Sex_System.csv")

#To test if there is a statistical difference between the microbial communities of two or more groups of samples.
#Null Hypothesis: there is no difference between the microbial communities of your groups of samples.

##Location
summary(vegan::anosim(bray_dist, tmp$Location, permutations = 999, strata =tmp$Origin))

#ANOSIM statistic R: 0.4419
#Significance: 1
#permutations = 999

##Conclusion: there is no difference between the microbial communities of Ascaris from FU and SH

##Experiment 1 vs Experiment 2
summary(vegan::anosim(bray_dist, tmp$Origin, permutations = 999, strata =tmp$System))

#ANOSIM statistic R: 0.4746
#Significance: 0.1
#permutations = 999

##Conclusion: there is no difference between the microbial communities of Ascaris from Experiment 1, Experiment 2 and SH

##Individuals
Individual.Asc.adonis<- vegan::anosim(bray_dist, tmp$System, permutations = 999, strata =tmp$Origin)

#ANOSIM statistic R: 0.3228
#Significance: 0.003 
#permutations = 999

##Conclusion: there is difference between the microbial communities of duodenum and ascaris in different individuals

## Non-metric multidimensional scaling
##With phyloseq
nmds.ordination<- ordinate(PS.Asc.Norm, method="NMDS", distance="bray", trymax= 50,
                           p.adjust.methods= "bonferroni", permutations = 999)

nmds.scores<- as.data.frame(vegan::scores(nmds.ordination))
nmds.scores<- cbind(nmds.scores, tmp)

genus.scores<- as.data.frame(vegan::scores(nmds.ordination, "species"))
genus.data<- as.data.frame(PS.Asc.Norm@tax_table)
genus.scores<- cbind(genus.scores, genus.data)
rm(genus.data)

genus.scores%>%
  dplyr::filter(!is.na(NMDS1), !is.na(NMDS2)) %>%
  dplyr::filter(!is.na(Genus))-> genus.scores 

genus.scores %>%
  dplyr::filter(rownames(genus.scores)%in%c("ASV28", "ASV29", "ASV36", "ASV50", "ASV16", "ASV46", 
                                            "ASV51", "ASV151", "ASV24", "ASV182", "ASV68"))-> genus.scores

nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= System, shape= WormSex), size=3) +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", shape= "Sex", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= genus.scores, aes(x = 0, y = 0, xend = (NMDS1)*2, yend = (NMDS2)*2),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), 
         color=F, arrow= F)+
  geom_text_repel(data = genus.scores, aes(x = (NMDS1)*2.2, y = (NMDS2)*2.2), label= genus.scores$Genus)+
  annotate("text", x = 1.5, y = 1.5, label= "ANOSIM (Individual) \n")+
  annotate("text", x = 1.5, y = 1.4, label= paste0(label = "R = ", round(Individual.Asc.adonis$statistic, digits = 3),
                                                   ", p = ", Individual.Asc.adonis$signif), color = "black")-> A3

##Transform dataset to determine contributors
#PS.Asc.clr<- subset_samples(PS.Asc, Replicate%in%Inf.Keep)
Inf.Keep<- row.names(PS.Asc@sam_data)
Inf.Keep<- alphadiv.Asc[rownames(alphadiv.Asc)%in%Inf.Keep, ]

Inf.Keep%>%
  dplyr::filter(System!= "Pig14")%>% #No jejunum
  dplyr::filter(System!= "Pig5")%>% #Just one ascaris
  dplyr::select(Replicate)%>%
  ungroup()%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.Asc.clr<- subset_samples(PS.Asc, Replicate%in%Inf.Keep)

PS.Asc.clr <- microbiome::transform(PS.Asc.clr, "clr") #Centered log ratio transformation
Ord.Asc.clr <- phyloseq::ordinate(PS.Asc.clr, "RDA") #principal components analysis

#Examine eigenvalues and % prop. variance explained
head(Ord.Asc.clr$CA$eig)
sapply(Ord.Asc.clr$CA$eig[1:6], function(x) x / sum(Ord.Asc.clr$CA$eig))

##ASVs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.Asc.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))
ind_cont_PCA1 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 ASVs contribute for the 20.9% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 ASVs contribute for the 17.4% of the variation in PC2

ind_cont_PCA_top.Asc <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge taxonomy
x<- as.data.frame(PS.Asc.clr@tax_table)
x<- x[rownames(x)%in%c(ind_cont_PCA_top.Asc$ASV),]
x%>%
  rownames_to_column("ASV")->x

ind_cont_PCA_top.Asc%>%
  distinct(ASV, .keep_all = T)%>%
  plyr::join(x, by="ASV")%>%
  column_to_rownames("ASV")-> ind_cont_PCA_top.Asc

##Taxa explaining variability
write.csv(ind_cont_PCA_top.Asc, "Tables/Q1_Principal_Taxa_Infected_Asc.csv")

x<- ind.coord[rownames(ind.coord)%in%c(rownames(ind_cont_PCA_top.Asc)),]
x%>%
  dplyr::filter(rownames(x)%in%c("ASV28", "ASV29", "ASV36", "ASV50", "ASV16", "ASV46", 
                                 "ASV51", "ASV151", "ASV24", "ASV182", "ASV68"))-> x

y<- ind_cont_PCA_top.Asc[rownames(ind_cont_PCA_top.Asc)%in%c(rownames(x)),]

x<- cbind(x, y)

PS.Asc.clr@sam_data$System<- fct_relevel(PS.Asc.clr@sam_data$System, 
                                            "Pig1","Pig2","Pig3","Pig4",
                                            "Pig5","Pig6","Pig7","Pig8","Pig9",
                                            "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH")
require(ggrepel)
plot_ordination(PS.Asc.clr, ordination = Ord.Asc.clr)+ 
  geom_point(size=3, aes(fill= System, shape= WormSex), color= "black")+
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", fill  = "Individual", shape= "Sex")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= x, aes(x = 0, y = 0, xend = (PC1)*10, yend = (PC2)*10),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), 
         color=F, arrow= F)+
  geom_text_repel(data = x, aes(x = (PC1*12), y = (PC2*12)), label= x$Genus)+
  xlab(paste0("PC 1 [", round(Ord.Asc.clr$CA$eig[1] / sum(Ord.Asc.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.Asc.clr$CA$eig[2] / sum(Ord.Asc.clr$CA$eig)*100, digits = 2), "%]"))

##Plot abundance of these ASVs
wh0 <- genefilter_sample(PS.Asc, filterfun_sample(function(x) x > 5), A=0.01*nsamples(PS.Asc))
PS.subset<- prune_taxa(wh0, PS.Asc)

sample_data(PS.subset)<- sample_data(tmp)

PS.subset <- subset_taxa(PS.subset, rownames(tax_table(PS.subset)) %in% c("ASV14", "ASV165", "ASV2","ASV20", "ASV22", "ASV293", "ASV304", "ASV333" ,"ASV35","ASV352",
                                                                          "ASV43","ASV454" ,"ASV55","ASV62","ASV68","ASV7","ASV87","ASV9","ASV93"))
##Changes by Location
phyloseq::psmelt(PS.subset) %>%
  dplyr::mutate(Origin = fct_relevel(Origin, 
                                     "Experiment_1", "Experiment_2", "Slaughterhouse"))%>%
  dplyr::mutate(Location = fct_relevel(Location, 
                                       "FU", "SH"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                                     "Pig1","Pig2","Pig3",
                                     "Pig10","Pig11", 
                                     "Pig12", "Pig13", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ Location)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Location")%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
 dplyr::mutate(y.position= c(0.075, 0.025, 0.7, 1.7, 0.06, 
                             0.012, 0.005, 0.006,0.035,0.007,
                             0.65,0.003,0.33,0.035,0.09,
                             0.9, 0.22, 0.6, 0.25))-> stats.test

###Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_ASV_Asc_Location.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3",
                              "Pig10","Pig11", "Pig12", "Pig13"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  ggplot(data = ., aes(x = Location, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = System, shape= WormSex), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Relative Abundance (%)", shape = "Sex") +
  facet_wrap(~ OTU, scales = "free")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> Supplementary_Figure_5A

##Check for dominant taxa for worms
###Effect of dominant taxa (Enterotype)
###Summarize to Genus
PS.Asc.Norm.Gen<-  tax_glom(PS.Asc.Norm, "Genus", NArm = F)

##Dominant taxa per sample 
dom.gen.Asc<- find.top.asv(PS.Asc.Norm.Gen, "Genus", 1) #--> Genus

dom.gen.Asc%>%
  dplyr::select(c(1,3,9))%>%
  dplyr::mutate(Gen.Abund= Abundance)%>%
  dplyr::mutate(Gen.Dom= Genus)%>%
  dplyr::select(c(Replicate, Gen.Abund, Gen.Dom))-> dom.gen.Asc 

tmp%>%
  left_join(dom.gen.Asc, by="Replicate")-> tmp

dom.phy.Asc<- find.top.asv(PS.Asc.Norm, "Phylum", 1) #--> Phylum

dom.phy.Asc%>%
  dplyr::select(c(1,3,5))%>%
  dplyr::mutate(Phy.Abund= Abundance)%>%
  dplyr::mutate(Phy.Dom= Phylum)%>%
  dplyr::select(c(Replicate, Phy.Abund, Phy.Dom))-> dom.phy.Asc 

tmp%>%
  left_join(dom.phy.Asc, by="Replicate")-> tmp 

tmp%>%
  dplyr::select(c(Replicate, Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom)) -> tmp.Dom

rownames(tmp.Dom)<- tmp.Dom$Replicate 

tmp.Dom%>%
  dplyr::select(c(Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom))%>%
  cbind(seg.data)-> seg.data 

tmp.Dom%>%
  dplyr::select(c(Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom))%>%
  cbind(nmds.scores)-> nmds.scores

#### Non-metric multidimensional scaling with Enterotype
##With phyloseq
##Location Ascaris
location.ascaris.anosim<- vegan::anosim(bray_dist, tmp$Gen.Dom, permutations = 999, strata =tmp$System)
#ANOSIM statistic R: 0.8459
#Significance: 0.001
#permutations = 999
nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= Gen.Dom, shape= WormSex), size=3) +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "B)", shape= "Sex", fill= "Enterotype")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotate("text", x = 1.1, y = 1.5, label= "ANOSIM (Enterotype) \n")+
  annotate("text", x = 1.1, y = 1.4, label= paste0(label = "R = ", round(location.ascaris.anosim$statistic, digits = 3),
                                                   ", p = ", location.ascaris.anosim$signif), color = "black")-> B3

###Estimate centroids for enterotypes
###Higher amount of variance is explain by experiment of origin 
mvd<- vegan::betadisper(bray_dist, tmp$Gen.Dom, type = "centroid")
anova(mvd)

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))

ggplot() + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Gen.Dom, shape= WormSex), size=3) +
  scale_shape_manual(values = c(23, 22), labels = c("Female", "Male"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "C)", shape= "Worm Sex", fill= "Enterotype")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))

worm.adonis.dom<- vegan::adonis(bray_dist~ Location + System + WormSex + Gen.Dom,
                            permutations = 999, data = tmp, na.action = F, strata = tmp$Origin)

##Store the result
foo<- as.data.frame(worm.adonis.dom$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Worm_dominant.csv")

###Compare abundance of dominants 
gen.Asc<- count.high.genus(x = PS.Asc.Norm.Gen, num = 0) ##

tmp%>%
  left_join(gen.Asc, by="Replicate")-> gen.Asc

###We have 3 big enterotypes:
###Enterotype 1: Clostridium sensu stricto 1
###Enterotype 2: Lactobacillus, Romboutsia, Streptococcus
###Enterotype 3: Escherichia-Shigella, Aeromonas

##Select for enterotype taxa 
gen.Asc%>%
  dplyr::filter(Genus%in%c("Clostridium sensu stricto 1", 
                           "Lactobacillus",  "Romboutsia", "Escherichia-Shigella", 
                           "Prevotella", "Streptococcus", "Aeromonas"))%>%
  dplyr::select(10:20, 23, 29)%>%
  mutate(Gen.Dom = fct_relevel(Gen.Dom,"Clostridium sensu stricto 1", 
                               "Lactobacillus",  "Escherichia-Shigella",
                               "Prevotella", "Romboutsia", "Streptococcus",
                                "Aeromonas"))-> Enterotype.abund

Enterotype.abund %>%
  group_by(Genus)%>%
  wilcox_test(Abundance ~ Gen.Dom)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "Gen.Dom")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Asc_Enterotype_abundances.csv")

stats.test%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(y.position= c(110, 120, 110, 120, 90, 100))-> stats.test 

# New facet label names for supp variable
supp.labs <- c("Clostridium s. stricto 1","Lactobacillus",
               "Escherichia-Shigella")
names(supp.labs) <- c("Clostridium sensu stricto 1", "Lactobacillus",
                      "Escherichia-Shigella")

Enterotype.abund%>%
  dplyr::filter(Genus%in%c("Clostridium sensu stricto 1", "Lactobacillus",
                           "Escherichia-Shigella"))%>%
  group_by(Genus)%>%
  ggplot(aes(x= Gen.Dom, y= Abundance))+
  facet_grid(~Genus, scales = "free")+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_jitter(size=3, width = 0.25, aes(fill= System, shape= WormSex), color= "black")+
  scale_shape_manual(values = c(22, 23), labels= c("Female", "Male"))+
  scale_fill_manual(values = pal.system)+
  ylab("Relative abundance (%)")+
  labs(tag= "B)", fill= "Individual", shape= "Worm Sex")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0, hide.ns = T,
                     tip.length = 0)-> Supplementary_Figure_5B

##Barplot by sample 
gen.Asc<- count.high.genus(x = PS.Asc.Norm.Gen, num = 1) ##Taxa less represented had less than 1% Relative abundance

tmp%>%
  left_join(gen.Asc, by="Replicate")-> gen.Asc

#set color palette to accommodate the number of genera
length(unique(gen.Asc$Genus))

#plot
gen.Asc%>%
  dplyr::filter(System== "SH")%>%
  mutate(Replicate = fct_relevel(Replicate, "SH.Ascaris.2", "SH.Ascaris.3", "SH.Ascaris.4", "SH.Ascaris.5",   
                                 "SH.Ascaris.7", "SH.Ascaris.8","SH.Ascaris.9", "SH.Ascaris.10", 
                                 "SH.Ascaris.11", "SH.Ascaris.12", "SH.Ascaris.13", "SH.Ascaris.14",  "SH.Ascaris.15",
                                 "SH.Ascaris.16",  "SH.Ascaris.17", "SH.Ascaris.18", "SH.Ascaris.19", "SH.Ascaris.20"))%>%
  ggplot(aes(x=Replicate, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack", width=.75) + 
  facet_grid(~WormSex, scales = "free", space = "free")+
  scale_fill_manual(values=c(tax.palette.SH)) + 
  theme_bw()+
  labs(tag= "C)")+
  ylab("Relative abundance (%)")+
  xlab("Sample ID")+
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=4))+
  theme(text = element_text(size=16), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())-> barplot.SH

##Save them individually
ggsave(file = "Figures/Q1_NMDS_Host_Ascaris.png", plot = A3, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Host_Ascaris.pdf", plot = A3, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Ascaris_ASVs.png", plot =Supplementary_Figure_5A, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Ascaris_ASVs.pdf", plot =Supplementary_Figure_5A, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Enterotype_Ascaris.png", plot = B3, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Enterotype_Ascaris.pdf", plot = B3, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Ascaris_Abundance.png", plot = Supplementary_Figure_5B, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Enterotype_Ascaris_Abundance.pdf", plot = Supplementary_Figure_5B, width = 10, height = 8, dpi = 450)

###Save figures 
##Al together
##Figure 5
Plot1<- cowplot::plot_grid(A3, B3, align = "hv", ncol= 1)

ggsave(file = "Figures/Figure_5.png", plot = Plot1, width = 10, height = 10, dpi = 450)
ggsave(file = "Figures/Figure_5.pdf", plot = Plot1, width = 10, height = 10, dpi = 450)
ggsave(file = "Figures/Figure_5.svg", plot = Plot1, width = 10, height = 10, dpi = 450)

##Supplementary Figure 5
Supplementary_Figure_5A+
  guides(shape = FALSE, fill= F)->Supplementary_Figure_5A

Plot2<- ggarrange(Supplementary_Figure_5A,Supplementary_Figure_5B, ncol =  2)

Plot2<- cowplot::plot_grid(Plot2, barplot.SH, align = "hv", ncol= 1)

ggsave(file = "Figures/Supplementary_Figure_5.png", plot = Plot2, width = 14, height = 12, dpi = 450)
ggsave(file = "Figures/Supplementary_Figure_5.pdf", plot = Plot2, width = 14, height = 12, dpi = 450)
ggsave(file = "Figures/Supplementary_Figure_5.svg", plot = Plot2, width = 14, height = 12, dpi = 450)

rm(A3,B3,Supplementary_Figure_5A, Supplementary_Figure_5B, Plot1, Plot2)

#############################Core microbiome analysis##################################################
##For Pigs
###############Infected
##Subset just the infected pigs 
tmp<- row.names(PS.pig.Norm@sam_data)
tmp<- alphadiv.pig[rownames(alphadiv.pig)%in%tmp, ]

tmp%>%
  dplyr::filter(InfectionStatus!= "Non_infected")%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.pig.Inf<- subset_samples(PS.pig.Norm, Replicate%in%Inf.Keep)

PS.rel <- microbiome::transform(PS.pig.Inf, "compositional")
core.taxa <- core_members(PS.rel, detection = 0.0001, prevalence = 50/100)

PS.core <- core(PS.rel, detection = 0.0001, prevalence = .5)

# Compute prevalence of each feature, store as data.frame
PrevAll <- apply(X = otu_table(PS.core),
                 MARGIN = ifelse(taxa_are_rows(PS.core), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
core.taxa.piginf <- data.frame(Prevalence = PrevAll,
                      TotalAbundance = taxa_sums(PS.core),
                      tax_table(PS.core))

##Transform prevalence into a porcentage 
core.taxa.piginf%>%
  dplyr::mutate(Prevalence= round((Prevalence/46)*100, digits = 2))-> core.taxa.piginf

core.taxa.piginf%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> core.taxa.piginf
  
###Site of infection
tmp%>%
  dplyr::filter(InfectionStatus!= "Non_infected")%>%
  dplyr::filter(Compartment== "Jejunum")%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.pig.JeInf<- subset_samples(PS.pig.Norm, Replicate%in%Inf.Keep)

PS.rel <- microbiome::transform(PS.pig.JeInf, "compositional")
core.taxa <- core_members(PS.rel, detection = 0.0001, prevalence = 50/100)

PS.core <- core(PS.rel, detection = 0.0001, prevalence = .5)

# Compute prevalence of each feature, store as data.frame
PrevAll <- apply(X = otu_table(PS.core),
                 MARGIN = ifelse(taxa_are_rows(PS.core), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
core.taxa.Jeinf <- data.frame(Prevalence = PrevAll,
                               TotalAbundance = taxa_sums(PS.core),
                               tax_table(PS.core))

core.taxa.Jeinf%>%
  dplyr::mutate(Prevalence= round((Prevalence/9)*100, digits = 2))->core.taxa.Jeinf

core.taxa.Jeinf%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> core.taxa.Jeinf
  
###############Subset just the uninfected pigs 
tmp%>%
  dplyr::filter(InfectionStatus!= "Infected")%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.pig.Unnf<- subset_samples(PS.pig.Norm, Replicate%in%Inf.Keep)

PS.rel <- microbiome::transform(PS.pig.Unnf, "compositional")
core.taxa <- core_members(PS.rel, detection = 0.0001, prevalence = 50/100)

PS.core <- core(PS.rel, detection = 0.0001, prevalence = .5)

# Compute prevalence of each feature, store as data.frame
PrevAll <- apply(X = otu_table(PS.core),
                 MARGIN = ifelse(taxa_are_rows(PS.core), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
core.taxa.piguninf <- data.frame(Prevalence = PrevAll,
                              TotalAbundance = taxa_sums(PS.core),
                              tax_table(PS.core))

core.taxa.piguninf%>%
  dplyr::mutate(Prevalence= round((Prevalence/16)*100, digits = 2))->core.taxa.piguninf

core.taxa.piguninf%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> core.taxa.piguninf

###Site of infection
tmp%>%
  dplyr::filter(InfectionStatus== "Non_infected")%>%
  dplyr::filter(Compartment== "Jejunum")%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.pig.Jeunf<- subset_samples(PS.pig.Norm, Replicate%in%Inf.Keep)

PS.rel <- microbiome::transform(PS.pig.Jeunf, "compositional")
core.taxa <- core_members(PS.rel, detection = 0.0001, prevalence = 50/100)

PS.core <- core(PS.rel, detection = 0.0001, prevalence = .5)

# Compute prevalence of each feature, store as data.frame
PrevAll <- apply(X = otu_table(PS.core),
                 MARGIN = ifelse(taxa_are_rows(PS.core), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
core.taxa.Jeunf <- data.frame(Prevalence = PrevAll,
                                 TotalAbundance = taxa_sums(PS.core),
                                 tax_table(PS.core))

core.taxa.Jeunf%>%
  dplyr::mutate(Prevalence= round((Prevalence/3)*100, digits = 2))->core.taxa.Jeunf

core.taxa.Jeunf%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> core.taxa.Jeunf

####################For Ascaris
##From Infected Pigs
tmp<- row.names(PS.Asc.Norm@sam_data)
tmp<- alphadiv.Asc[rownames(alphadiv.Asc)%in%tmp, ]

tmp%>%
  dplyr::filter(Origin!= "Slaughterhouse")%>%
  dplyr::filter(System!= "Pig14")%>% #No jejunum
  dplyr::filter(System!= "Pig5")%>% #Just one ascaris
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.Asc.Inf<- subset_samples(PS.Asc.Norm, Replicate%in%Inf.Keep)

PS.rel <- microbiome::transform(PS.Asc.Inf, "compositional")
core.taxa <- core_members(PS.rel, detection = 0.0001, prevalence = 50/100)

PS.core <- core(PS.rel, detection = 0.0001, prevalence = .5)

# Compute prevalence of each feature, store as data.frame
PrevAll <- apply(X = otu_table(PS.core),
                 MARGIN = ifelse(taxa_are_rows(PS.core), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
core.taxa.Ascinf <- data.frame(Prevalence = PrevAll,
                              TotalAbundance = taxa_sums(PS.core),
                              tax_table(PS.core))

core.taxa.Ascinf%>%
  dplyr::mutate(Prevalence= round((Prevalence/41)*100, digits = 2))->core.taxa.Ascinf

core.taxa.Ascinf%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))-> core.taxa.Ascinf

##From Infected Pigs SH
tmp%>%
  dplyr::filter(Origin== "Slaughterhouse")%>%
  dplyr::select(Replicate)-> Inf.Keep

Inf.Keep<- Inf.Keep$Replicate

PS.Asc.SH<- subset_samples(PS.Asc.Norm, Replicate%in%Inf.Keep)

PS.rel <- microbiome::transform(PS.Asc.SH, "compositional")
core.taxa <- core_members(PS.rel, detection = 0.0001, prevalence = 50/100)

PS.core <- core(PS.rel, detection = 0.0001, prevalence = .5)

# Compute prevalence of each feature, store as data.frame
PrevAll <- apply(X = otu_table(PS.core),
                 MARGIN = ifelse(taxa_are_rows(PS.core), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
core.taxa.AscSH  <- data.frame(Prevalence = PrevAll,
                               TotalAbundance = taxa_sums(PS.core),
                               tax_table(PS.core))

core.taxa.AscSH%>%
  dplyr::mutate(Prevalence= round((Prevalence/18)*100, digits = 2))->core.taxa.AscSH 

core.taxa.AscSH%>%
  unite(Bacteria_name, c("Genus", "Species"), remove = F)%>%
  dplyr::filter(Bacteria_name!= "NA_NA")%>%
  dplyr::mutate(Bacteria_name = gsub("_NA", " sp.", basename(Bacteria_name)))%>%
  dplyr::mutate(Bacteria_name = gsub("_", " ", basename(Bacteria_name)))->core.taxa.AscSH

###Determine which ASVs are distinctive for each microbiome
##Infected vs non infected (Jejunum)
diff.piginf<- setdiff(rownames(core.taxa.Jeinf), rownames(core.taxa.Jeunf)) ##Just in infected
diff.piguninf<- setdiff(rownames(core.taxa.Jeunf), rownames(core.taxa.Jeinf)) ## Just in non infected

int.inf<- intersect(rownames(core.taxa.Jeinf), rownames(core.taxa.Jeunf))

core.inf<- c(diff.piginf, diff.piguninf, int.inf)

PS.pig.diff<- microbiome::transform(PS.pig.Norm, "compositional")
PS.pig.diff<- subset_samples(PS.pig.diff, Compartment== "Jejunum")
PS.pig.diff<-subset_taxa(PS.pig.diff, rownames(t(PS.pig.diff@otu_table))%in%core.inf)

##Make a table
x<- as.data.frame(PS.pig.diff@tax_table[rownames(PS.pig.diff@tax_table)%in%c(diff.piginf),])
x$Origin<- as.factor("Infected")

y<- as.data.frame(PS.pig.Jeunf@tax_table[rownames(PS.pig.Jeunf@tax_table)%in%c(diff.piguninf),])
y$Origin<- as.factor("Non_Infected")

##What is different between infected pigs and ascaris 
diff.Ascinf<- setdiff(core.taxa.Ascinf$OTU, core.taxa.piginf$OTU) ##just in Ascaris from infected pigs 

diff.inf<- unique(c(diff.piginf, diff.piguninf, diff.Ascinf))
PS.PA.diff<- microbiome::transform(PS.PA.Norm, "compositional")
PS.PA.diff<- subset_samples(PS.PA.diff, Compartment%in%c("Jejunum", "Ascaris"))
PS.PA.diff<- subset_samples(PS.PA.diff, Origin!= "Slaughterhouse")
PS.PA.diff<-subset_taxa(PS.PA.diff, rownames(t(PS.PA.diff@otu_table))%in%diff.inf)
PS.PA.diff<-subset_taxa(PS.PA.diff, Genus!= "NA")

z<- as.data.frame(PS.Asc.Inf@tax_table[rownames(PS.Asc.Inf@tax_table)%in%c(diff.Ascinf),])
z$Origin<- as.factor("Ascaris")

inf.uninf.asv<- rbind(x,y,z)

write.csv(inf.uninf.asv, "Tables/Q2_Core_ASVs_Inf_Asc.csv")

##FU vs SH
diff.AscFU<- setdiff(core.taxa.Ascinf$OTU, core.taxa.AscSH$OTU) ##just in Ascaris from FU 
diff.AscSH<- setdiff(core.taxa.AscSH$OTU, core.taxa.Ascinf$OTU) ##just in Ascaris from SH 

##Make a table
x<- as.data.frame(PS.Asc.Inf@tax_table[rownames(PS.Asc.Inf@tax_table)%in%c(diff.AscFU),])
x$Origin<- as.factor("Ascaris_FU")

y<- as.data.frame(PS.Asc.SH@tax_table[rownames(PS.Asc.SH@tax_table)%in%c(diff.AscSH),])
y$Origin<- as.factor("Ascaris_SH")

fu.sh.asv<- rbind(x,y)

write.csv(fu.sh.asv, "Tables/Q2_Core_ASVs_FU_SH.csv")

#####################Alternative Figure 4################################################
#Merge PA with SH worms 
PS.JejWorms<- merge_phyloseq(PS.JejAsc, PS.Asc.SH)
bray_dist<- phyloseq::distance(PS.JejWorms, 
                               method="bray", weighted=F)

ordination<- ordinate(PS.JejWorms,
                      method="PCoA", distance="bray")

tmp<- row.names(PS.JejWorms@sam_data)
tmp.1<- alphadiv.PA[rownames(alphadiv.PA)%in%tmp, ]
tmp.1%>%
  dplyr::filter(Compartment=="Jejunum")%>%
  dplyr::mutate(Location = case_when(Origin %in% c("Experiment_1", "Experiment_2")  ~ "FU",
                                     Origin == "Slaughterhouse" ~ "SH"))-> tmp.1

tmp.2<- alphadiv.Asc[rownames(alphadiv.Asc)%in%tmp, ]

tmp<- rbind(tmp.1, tmp.2)

tmp%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))-> tmp

jejunum.ascaris.all.adonis<- vegan::adonis(bray_dist~ AnimalSpecies + System,
                                        permutations = 999, data = tmp, na.action = F, strata = tmp$Origin)

##Store data
foo<- as.data.frame(jejunum.ascaris.all.adonis$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Jejunum_All_Ascaris.csv")

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
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", shape= "Host-Parasite", color= "Origin of samples", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))

###Distance to host Jejunum
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

BC.JejAscSH<- as.data.frame(x[c(Inf.asc.Keep), c(Inf.pig.Keep)])

BC.JejAscSH %>%
  rownames_to_column(var = "Replicate")%>%
  gather("Pig1.Jejunum", "Pig10.Jejunum", "Pig11.Jejunum", "Pig12.Jejunum", "Pig13.Jejunum", "Pig2.Jejunum",
         "Pig3.Jejunum" , key = Host, value = BC_dist)%>%
  left_join(tmp, by= "Replicate")-> BC.JejAscSH

#Are the Ascaris from SH closer than Ascaris from FU to FU Pigs jejunum
BC.JejAscSH%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Host = fct_relevel(Host, 
                            "Pig1.Jejunum", "Pig2.Jejunum",
                            "Pig3.Jejunum", "Pig10.Jejunum", "Pig11.Jejunum", "Pig12.Jejunum", "Pig13.Jejunum"))%>%
  dplyr::group_by(Host)%>%
  wilcox_test(BC_dist ~ System)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "System")%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(xmax= 8)%>%
  dplyr::mutate(xmin = case_when(xmin== 1  ~ 1, xmin== 2  ~ 2, xmin== 3  ~ 3,
                                 xmin== 10  ~ 4, xmin== 11  ~ 5, xmin== 12  ~ 6,
                                 xmin== 13  ~ 7))%>%
  dplyr::mutate(y.position= c(1.02, 1.05, 1.08,
                              1.02, 1.05, 1.08, 1.12, 1.16,
                              1.02, 1.05, 1.08, 1.12, 1.16, 1.19,
                              1.02, 1.05, 1.08, 1.12,
                              1.02, 1.05, 1.08, 1.12, 1.16,
                              1.02, 1.05, 1.08, 1.12, 1.16, 1.19,
                              1.02, 1.05, 1.08, 1.12, 1.16))-> stats.test
BC.JejAscSH%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3","Pig4",
                              "Pig5","Pig6","Pig7","Pig8","Pig9",
                              "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH"))%>%
  mutate(Host = fct_relevel(Host, 
                            "Pig1.Jejunum", "Pig2.Jejunum",
                            "Pig3.Jejunum", "Pig10.Jejunum", "Pig11.Jejunum", "Pig12.Jejunum", "Pig13.Jejunum"))%>%
  dplyr::group_by(Host)%>%
  ggplot(aes(x= System, y= BC_dist))+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_point(size=3, aes(fill= System), shape=21, color= "black")+
  scale_fill_manual(values = pal.system)+
  facet_grid(~Host, scales = "free", space = "free")+
  ylab("Bray-Curtis dissimilarity \n between Ascaris and Jejunum")+
  labs(tag= "B)", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  stat_pvalue_manual(stats.test, hide.ns = T,
                     tip.length = 0)-> BC.JejAscSH.plot

#To test if there is a statistical difference between the microbial communities of two or more groups of samples.
#Null Hypothesis: there is no difference between the microbial communities of your groups of samples.

##Jejunum vs Ascaris
summary(vegan::anosim(bray_dist, tmp$AnimalSpecies, permutations = 999, strata =tmp$Origin))
#ANOSIM statistic R: -0.083
#Significance: 0.559
#permutations = 999

##Conclusion: there is no difference between the microbial communities of all Ascaris and Jejunum

##FU vs SH
summary(vegan::anosim(bray_dist, tmp$Location, permutations = 999, strata =tmp$System))
#ANOSIM statistic R: 0.4253
#Significance: 0.1
#permutations = 999

##Conclusion: there is no difference between the microbial communities of both locations

##Individuals
Jejunum.ascaris.all.anosim<- vegan::anosim(bray_dist, tmp$System, permutations = 999, strata =tmp$Origin)
#ANOSIM statistic R: 0.2911 
#Significance: 0.003 
#permutations = 999

##Conclusion: there is difference between the microbial communities of jejunum and Ascaris in different individuals

## Non-metric multidimensional scaling
##With phyloseq
nmds.ordination<- ordinate(PS.JejWorms, method="NMDS", distance="bray",  trymax= 100,
                           p.adjust.methods= "bonferroni", permutations = 999)

nmds.scores<- as.data.frame(vegan::scores(nmds.ordination))
nmds.scores<- cbind(nmds.scores, tmp)

genus.scores<- as.data.frame(vegan::scores(nmds.ordination, "species"))
genus.data<- as.data.frame(PS.JejWorms@tax_table)
genus.scores<- cbind(genus.scores, genus.data)
rm(genus.data)

genus.scores%>%
  dplyr::filter(!is.na(NMDS1), !is.na(NMDS2)) %>%
  dplyr::filter(!is.na(Genus))-> genus.scores 

genus.scores %>%
  dplyr::filter(rownames(genus.scores)%in%c("ASV36", "ASV50", "ASV29", "ASV28", "ASV26", "ASV64", 
                                            "ASV13", "ASV16", "ASV24", "ASV47", "ASV15", "ASV21", "ASV68"))-> genus.scores

nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= System, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", shape= "Host-Parasite", color= "Origin of samples", fill= "Individual")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= genus.scores, aes(x = 0, y = 0, xend = (NMDS1)*2, yend = (NMDS2)*2),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), 
         color=F, arrow= F)+
  geom_text_repel(data = genus.scores, aes(x = (NMDS1)*2.2, y = (NMDS2)*2.2), label= genus.scores$Genus)+
  annotate("text", x = 2.5, y =  1.5, label= "ANOSIM (Individual) \n")+
  annotate("text", x = 2.5, y =  1.4, label= paste0(label = "R = ", round(Jejunum.ascaris.all.anosim$statistic, digits = 3),
                                                   ", p = ", Jejunum.ascaris.all.anosim$signif), color = "black")-> nmds.JejAsc

##Transform dataset to determine contributors
PS.JejWorms.clr <- microbiome::transform(PS.JejWorms, "clr") #Centered log ratio transformation
Ord.JejWorms.clr <- phyloseq::ordinate(PS.JejWorms.clr, "RDA") #principal components analysis

#Examine eigenvalues and % prop. variance explained
head(Ord.JejWorms.clr$CA$eig)
sapply(Ord.JejWorms.clr$CA$eig[1:6], function(x) x / sum(Ord.JejWorms.clr$CA$eig))

##ASVs contributing into PC1 and PC2
ind.coord <- data.frame(Ord.JejWorms.clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))
ind_cont_PCA1 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA1")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA1_top

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##25 ASVs contribute for the 21.9% of the variation in PC1

ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2 %>% 
  rownames_to_column("ASV") %>% 
  mutate(Component= "PCoA2")%>%
  arrange(desc(PCA))%>%
  slice_head(n = 25)-> ind_cont_PCA2_top

sum(ind_cont_PCA2_top$PCA) / sum(ind_cont_PCA2$PCA)
##25 ASVs contribute for the 55.6% of the variation in PC2

ind_cont_PCA_top.JejAscAll <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top)

##Merge taxonomy
x<- as.data.frame(PS.JejWorms.clr@tax_table)
x<- x[rownames(x)%in%c(ind_cont_PCA_top.JejAscAll$ASV),]
x%>%
  rownames_to_column("ASV")->x

ind_cont_PCA_top.JejAscAll%>%
  distinct(ASV, .keep_all = T)%>%
  plyr::join(x, by="ASV")%>%
  column_to_rownames("ASV")-> ind_cont_PCA_top.JejAscAll

##Taxa explaining variability
write.csv(ind_cont_PCA_top.JejAscAll, "Tables/Q1_Principal_Taxa_Infected_JejAscAll.csv")

x<- ind.coord[rownames(ind.coord)%in%c(rownames(ind_cont_PCA_top.JejAscAll)),]
x%>%
  dplyr::filter(rownames(x)%in%c("ASV36", "ASV50", "ASV29", "ASV28", "ASV22", "ASV47", 
                                 "ASV20", "ASV2", "ASV7", "ASV9"))-> x

y<- ind_cont_PCA_top.JejAscAll[rownames(ind_cont_PCA_top.JejAscAll)%in%c(rownames(x)),]

x<- cbind(x, y)

PS.JejWorms.clr@sam_data$System<- fct_relevel(PS.JejWorms.clr@sam_data$System, 
                                            "Pig1","Pig2","Pig3","Pig4",
                                            "Pig5","Pig6","Pig7","Pig8","Pig9",
                                            "Pig10","Pig11", "Pig12", "Pig13", "Pig14", "SH")
require(ggrepel)
plot_ordination(PS.JejWorms.clr, ordination = Ord.JejWorms.clr)+ 
  geom_point(size=3, aes(fill= System, shape= InfectionStatus), color= "black")+
  scale_shape_manual(values = c(24, 21), labels = c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "A)", fill  = "Individual", shape= "Host-Parasite")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_segment(data= x, aes(x = 0, y = 0, xend = (PC1)*10, yend = (PC2)*10),
               arrow = arrow(length = unit(0.2, "cm")))+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), 
         color=F, arrow= F)+
  geom_text_repel(data = x, aes(x = (PC1*12), y = (PC2*12)), label= x$Genus)+
  xlab(paste0("PC 1 [", round(Ord.JejWorms.clr$CA$eig[1] / sum(Ord.JejWorms.clr$CA$eig)*100, digits = 2), "%]"))+
  ylab(paste0("PC 2 [", round(Ord.JejWorms.clr$CA$eig[2] / sum(Ord.JejWorms.clr$CA$eig)*100, digits = 2), "%]"))-> PCA.JejAsc.All

##Plot abundance of these ASVs
wh0 <- genefilter_sample(PS.JejWorms, filterfun_sample(function(x) x > 5), A=0.01*nsamples(PS.JejWorms))
PS.subset<- prune_taxa(wh0, PS.JejWorms)

##
require("indicspecies")

##Changes by compartment
phyloseq::psmelt(PS.subset) %>%
  dplyr::mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  dplyr::mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3",
                              "Pig10","Pig11", "Pig12", "Pig13", "SH"))%>%
  dplyr::mutate(Location = case_when(Origin %in% c("Experiment_1", "Experiment_2")  ~ "FU",
                                     Origin == "Slaughterhouse" ~ "SH"))%>%
  tidyr::unite(col = "CL", c(Compartment, Location), sep = "_", remove = F)%>%
  dplyr::mutate(CL = fct_relevel(CL, "Jejunum_FU", "Ascaris_FU", "Ascaris_SH"))%>%
  dplyr::mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  dplyr::group_by(OTU)%>%
  wilcox_test(Abundance ~ CL)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "CL")%>%
  dplyr::filter(p.adj.signif!= "ns")-> stats.test
  
##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Abundance_ASV_JejAscAll.csv")

phyloseq::psmelt(PS.subset) %>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Jejunum", "Ascaris"))%>%
  mutate(System = fct_relevel(System, 
                              "Pig1","Pig2","Pig3",
                              "Pig10","Pig11", "Pig12", "Pig13", "SH"))%>%
  mutate(Abundance = (Abundance/1E6)*100)%>% ##Transform to relative abundance 
  ggplot(data = ., aes(x = Compartment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(fill = System, shape= InfectionStatus), height = 0, width = .2, size= 3, color= "black") +
  scale_shape_manual(values = c(24,21), labels = c("Infected Pig", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.title.x=element_blank())+
  labs(x = "", y = "Relative Abundance (%)", shape = "Infection status") +
  facet_wrap(~ OTU, scales = "free")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0.005, hide.ns = T,
                     tip.length = 0)-> ASV.JejAsc.All

###Effect of dominant taxa (Enterotype)
###Summarize to Genus
PS.JejWorms.Gen<-  tax_glom(PS.JejWorms, "Genus", NArm = F)

##Dominant taxa per sample 
dom.gen.JejWorms<- find.top.asv(PS.JejWorms, "Genus", 1) #--> Genus

dom.gen.JejWorms%>%
  dplyr::select(c(1,3,9))%>%
  dplyr::mutate(Gen.Abund= Abundance)%>%
  dplyr::mutate(Gen.Dom= Genus)%>%
  dplyr::select(c(Replicate, Gen.Abund, Gen.Dom))-> dom.gen.JejWorms 

tmp%>%
  left_join(dom.gen.JejWorms, by="Replicate")-> tmp

dom.phy.JejWorms<- find.top.asv(PS.JejWorms, "Phylum", 1) #--> Phylum

dom.phy.JejWorms%>%
  dplyr::select(c(1,3,5))%>%
  dplyr::mutate(Phy.Abund= Abundance)%>%
  dplyr::mutate(Phy.Dom= Phylum)%>%
  dplyr::select(c(Replicate, Phy.Abund, Phy.Dom))-> dom.phy.JejWorms 

tmp%>%
  left_join(dom.phy.JejWorms, by="Replicate")-> tmp 

tmp%>%
  dplyr::select(c(Replicate, Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom)) -> tmp.Dom

rownames(tmp.Dom)<- tmp.Dom$Replicate 

tmp.Dom%>%
  dplyr::select(c(Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom))%>%
  cbind(seg.data)-> seg.data 

tmp.Dom%>%
  dplyr::select(c(Gen.Abund, Gen.Dom, Phy.Abund, Phy.Dom))%>%
  cbind(nmds.scores)-> nmds.scores

#### Non-metric multidimensional scaling with Enterotype
##With phyloseq
##Jejunum vs  all Ascaris
jejunum.all.ascaris.anosim<- vegan::anosim(bray_dist, tmp$Gen.Dom, permutations = 999, strata =tmp$System)
#ANOSIM statistic R: 0.386
#Significance: 0.001
#permutations = 999
nmds.scores%>%
  ggplot(aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(fill= Gen.Dom, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "C)", shape= "Host-Parasite", fill= "Enterotype")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))))+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotate("text", x = 0.5, y = -1.0, label= "ANOSIM (Enterotype) \n")+
  annotate("text", x = 0.5, y = -0.6, label= paste0(label = "R = ", round(jejunum.all.ascaris.anosim$statistic, digits = 3),
                                                   ", p = ", jejunum.all.ascaris.anosim$signif), color = "black")-> nmds.JejAsc.Entero

###Estimate centroids for enterotypes
###Higher amount of variance is explain by experiment of origin 
mvd<- vegan::betadisper(bray_dist, tmp$Gen.Dom, type = "centroid")
anova(mvd)

##Extract centroids and vectors 
centroids<-data.frame(grps=rownames(mvd$centroids),data.frame(mvd$centroids))

ggplot() + 
  #geom_point(data=centroids[c(1,4),1:4], aes(x=PCoA1,y=PCoA2, color= grps, group=grps), size=4, shape= 4) +
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, fill= Gen.Dom, shape= Compartment), size=3) +
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = tax.palette)+
  labs(tag= "C)", shape= "Host-Parasite", fill= "Enterotype")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= F)+
  #stat_ellipse(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, color= Gen.Dom), linetype = 2)+
  scale_color_manual(values = tax.palette)+
  theme_bw()+
  theme(text = element_text(size=16))+
  xlab(paste0("PCo 1 [", round(ordination$values[1,2]*100, digits = 2), "%]"))+
  ylab(paste0("PCo 2 [", round(ordination$values[2,2]*100, digits = 2), "%]"))

jej.ascaris.dom.adonis.all<- vegan::adonis(bray_dist~ AnimalSpecies + System + Gen.Dom,
                                            permutations = 999, data = tmp, na.action = F, strata = tmp$Origin)

##Store data
foo<- as.data.frame(jej.ascaris.dom.adonis.all$aov.tab)
#write.csv(foo, file = "Tables/Q1_Adonis_Jejunum_Ascaris_All_Dom.csv")

##Barplot by sample 
gen.JejWorm<- count.high.genus(x = PS.JejWorms.Gen, num = 1) ##Taxa less represented had less than 1% Relative abundance

tmp%>%
  left_join(gen.JejWorm, by="Replicate")-> gen.JejWorm

#set color palette to accommodate the number of genera
length(unique(gen.JejWorm$Genus))

#plot
gen.JejWorm%>%
  dplyr::filter(System== "SH")%>%
  mutate(Replicate = fct_relevel(Replicate, "SH.Ascaris.2", "SH.Ascaris.3", "SH.Ascaris.4", "SH.Ascaris.5",   
                                 "SH.Ascaris.7", "SH.Ascaris.8","SH.Ascaris.9", "SH.Ascaris.10", 
                                 "SH.Ascaris.11", "SH.Ascaris.12", "SH.Ascaris.13", "SH.Ascaris.14",  "SH.Ascaris.15",
                                 "SH.Ascaris.16",  "SH.Ascaris.17", "SH.Ascaris.18", "SH.Ascaris.19", "SH.Ascaris.20"))%>%
  ggplot(aes(x=Replicate, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack", width=.75) + 
  facet_grid(~WormSex, scales = "free", space = "free")+
  scale_fill_manual(values=c(tax.palette.SH)) + 
  theme_bw()+
  labs(tag= "A)")+
  ylab("Relative abundance (%)")+
  xlab("Sample ID")+
  theme(legend.position="bottom")+
  guides(fill=guide_legend(nrow=4))+
  theme(text = element_text(size=16), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())-> barplot.SH

###Compare abundance of dominants 
gen.JejWorm<- count.high.genus(x = PS.JejWorms.Gen, num = 0) ##Taxa less represented had less than 1% Relative abundance

tmp%>%
  left_join(gen.JejWorm, by="Replicate")-> gen.JejWorm

###Based on PCoA we have 7 Enterotypes:
###Lactobacillus, Escherichia-Shigella, Prevotella, Streptococcus, Clostridium sensu stricto 1, Aeromonas, Romboutsia

gen.JejWorm%>%
  dplyr::filter(Genus%in%c("Clostridium sensu stricto 1", 
                           "Lactobacillus",  "Romboutsia", "Escherichia-Shigella", 
                           "Prevotella", "Streptococcus", "Aeromonas"))%>%
  dplyr::select(10:25, 29)%>%
  mutate(Gen.Dom = fct_relevel(Gen.Dom,"Clostridium sensu stricto 1", "Lactobacillus",
                               "Escherichia-Shigella", "Prevotella", "Romboutsia", "Streptococcus", "Aeromonas"))-> Enterotype.abund

Enterotype.abund %>%
  group_by(Genus)%>%
  wilcox_test(Abundance ~ Gen.Dom)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "Gen.Dom")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_JejAscAll_Enterotype_abundances.csv")

stats.test%>%
  dplyr::filter(p.adj.signif!= "ns")%>%
  dplyr::mutate(y.position= c(100, 105,110, 100, 105, 90, 95, 100))-> stats.test

Enterotype.abund%>%
  group_by(Genus)%>%
  ggplot(aes(x= Gen.Dom, y= Abundance))+
  facet_grid(~Genus, scales = "free", space = "free")+
  geom_boxplot(color= "black", alpha= 0.5, outlier.shape=NA)+
  geom_jitter(size=3, width = 0.25, aes(fill= System, shape= Compartment), color= "black")+
  scale_shape_manual(values = c(24, 21), labels= c("Infected Pig (Jejunum)", "Ascaris"))+
  scale_fill_manual(values = pal.system)+
  ylab("Relative abundance (%)")+
  labs(tag= "B)", fill= "Individual", shape= "Host-Parasite")+
  guides(fill = guide_legend(override.aes=list(shape=c(21))), color= FALSE)+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0, step.increase = 0, hide.ns = T,
                     tip.length = 0)-> E2

##Save them individually
ggsave(file = "Figures/Q1_NMDS_Host_Jejunum_Ascaris_All.png", plot = nmds.JejAsc, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Host_Jejunum_Ascaris_All.pdf", plot = nmds.JejAsc, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_Distance_All.png", plot = BC.JejAscSH.plot, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_Distance_All.pdf", plot = BC.JejAscSH.plot, width = 10, height = 8, dpi = 450)
#ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_ASVs_All.png", plot = Sup3, width = 10, height = 8, dpi = 450)
#ggsave(file = "Figures/Q1_Host_Jejunum_Ascaris_ASVs_All.pdf", plot = Sup3, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Enterotype_Jejunum_Ascaris_All.png", plot = nmds.JejAsc.Entero, width = 10, height = 8, dpi = 450)
ggsave(file = "Figures/Q1_NMDS_Enterotype_Jejunum_Ascaris_All.pdf", plot = nmds.JejAsc.Entero, width = 10, height = 8, dpi = 450)
#ggsave(file = "Figures/Q1_Composition_Jejunum_Ascaris_Abundance_All.png", plot = barplot.SH, width = 10, height = 8, dpi = 450)
#ggsave(file = "Figures/Q1_Composition_Jejunum_Ascaris_Abundance_All.pdf", plot = barplot.SH, width = 10, height = 8, dpi = 450)
#ggsave(file = "Figures/Q1_Enterotype_Jejunum_Ascaris_Abundance_All.png", plot = E2, width = 10, height = 8, dpi = 450)
#ggsave(file = "Figures/Q1_Enterotype_Jejunum_Ascaris_Abundance_All.pdf", plot = E2, width = 10, height = 8, dpi = 450)

###Save figure 4 
##Al together
BC.JejAscSH.plot+
  guides(fill = FALSE)-> BC.JejAscSH.plot

Plot1<- ggarrange(nmds.JejAsc, BC.JejAscSH.plot, nmds.JejAsc.Entero, ncol=1, align = "v")

ggsave(file = "Figures/Figure_4.2.png", plot = Plot1, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Figure_4.2.pdf", plot = Plot1, width = 10, height = 12, dpi = 450)
ggsave(file = "Figures/Figure_4.2.svg", plot = Plot1, width = 10, height = 12, dpi = 450)

## Supplementary composition

#Plot2<- ggarrange(barplot.SH, E2, ncol=1)

#ggsave(file = "Figures/Supplementary_Figure_4.2.png", plot = Plot2, width = 15, height = 12, dpi = 450)
#ggsave(file = "Figures/Supplementary_Figure_4.2.pdf", plot = Plot2, width = 15, height = 12, dpi = 450)
#ggsave(file = "Figures/Supplementary_Figure_4.2.svg", plot = Plot2, width = 15, height = 12, dpi = 450)
