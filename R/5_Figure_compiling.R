##Project: Ascaris - Pig Microbiome
##Aim: Compile figures
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
library("ggrepel")

##For Figure 3
Beta.div.JA<- readRDS("Figures/Beta_Jej_Asc_part1.RDS")
B<- readRDS("Figures/Q1_Diff_Abundance_JejAsc.RDS")

Figure.3<- grid.arrange(Beta.div.JA, B, widths = c(5, 5),
                           layout_matrix = rbind(c(1, 1),
                                                 c(2, 2)))

ggsave(file = "Figures/Figure.3.Final.pdf", plot = Figure.3, width = 16, height = 20, dpi = 600)
ggsave(file = "Figures/Figure.3.Final.png", plot = Figure.3, width = 16, height = 20, dpi = 600)
ggsave(file = "Figures/Figure.3.Final.svg", plot = Figure.3, width = 16, height = 20, dpi = 600)

##For Figure 4
WormsA.V2<- readRDS("Figures/Alpha_Sex_Worms_A_V2.RDS") ##with sex coloring 
WormsA.V2+
  guides(fill=F)-> WormsA.V2

A3<- readRDS("Figures/Beta_Sex_Worms_A_V2.RDS")
C3<- readRDS("Figures/Diff_Phylum_Sex_Worms_A_V2.RDS") ##with sex coloring 
C3+
  guides(fill=F)-> C3

D<- readRDS("Figures/Q1_Diff_Abundance_AscFM.RDS") ##to compile with other
barplot.Asc<- readRDS("Figures/RA_Sex_Worms_A_V2.RDS") ##with sex coloring 


Figure.4.Top<- grid.arrange(WormsA.V2, A3, C3, D, widths = c(5, 8),
                        layout_matrix = rbind(c(1, 2),
                                              c(3, 4)))

Figure.4<-grid.arrange(Figure.4.Top, barplot.Asc)

ggsave(file = "Figures/Figure.4.Final.pdf", plot = Figure.4, width = 19, height = 22, dpi = 600)
ggsave(file = "Figures/Figure.4.Final.png", plot = Figure.4, width = 16, height = 20, dpi = 600)
ggsave(file = "Figures/Figure.4.Final.svg", plot = Figure.4, width = 19, height = 22, dpi = 600)