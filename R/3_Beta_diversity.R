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