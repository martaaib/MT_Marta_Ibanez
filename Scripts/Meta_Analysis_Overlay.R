## R version 4.0.4 (2021-02-15)
##
## R script for overlay of DEGs and DMGs, meta-analysis and graphs
##
##
##  Copyright: Marta Ibañez Lligoña (marta.ib@live.com)
##

message("+-------------------------------------------------------------------------------+")
message("+       Load the packages and basic settings                                    +")
message("+-------------------------------------------------------------------------------+")
## ----- Packages ----- ##
library(readxl)
library(ggplot2)
library(forestplot)
library(tidyverse)
library(meta)
library(dplyr)
library("ggrepel")
library(RColorBrewer)
require(plyr)
library(reshape2)
library("scales")
library(clusterProfiler)
library(mygene)
library("ggpubr")
library(gplots)
library(rpsychi)
library(fishmethods)
library(paswr)

## ------ Set working directory ----- ##
setwd("~/Documents/GitHub/MT_Marta_Ibanez/")

message("+-------------------------------------------------------------------------------+")
message("+              TRANSCRIPTOME-WIDE STUDIES: DEGs                                 +")
message("+-------------------------------------------------------------------------------+")

data_studies <- read_excel("Data/DEGs_Studies_full_list.xlsx")

## ----- Overlay of DEGs ------ ##
Study1 <- data_studies[1:405,1]
colnames(Study1) <- "Genes"
Study2 <- data_studies[1:281,2]
colnames(Study2) <- "Genes"
Study3 <- data_studies[1:243,3]
colnames(Study3) <- "Genes"
Study4 <- data_studies[,4]
colnames(Study4) <- "Genes"

## Create matrix to define presence of gene in a study
all_genes <- join_all(list(Study1, Study2, Study3, Study4),  by = 'Genes', type = 'full')
all_genes <- na.omit(unique(all_genes))
key <- unlist(unique(c(all_genes)))
studies <- data.frame(
  key = key,
  do.call(
    cbind,
    lapply(
      data_studies[,1:4],
      function(x) ifelse(key %in% x == TRUE, 1, 0))))
## Convert to numeric matrix
studies <- unique(studies)
new_data <- as.matrix(sapply( studies[,2:ncol(studies)], as.numeric))
rownames(new_data) <- studies$key

## -- Getting most matched genes, those that are present in more than 3 studies -- ##
int <- rowSums(new_data) >= 3
most_matched_genes <- new_data[int,]
nrow(most_matched_genes)

## ------- Forestplot ------- ##
participants_BMI <- read_excel("Data/Transcriptome_BMI_CI_full_studies.xlsx")
tabletext <- cbind(c("Study",participants_BMI$Study_name), c("Number of control samples", participants_BMI$Nc), c("Number of GDM samples", participants_BMI$Ne))
forestplot(labeltext = tabletext, 
           fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
           #legend_args = fpLegend(pos = list(x = .2, y = 0.25), 
           #gp = gpar(col = "#CCCCCC", fill = "#F9F9F9")),
           legend = c("Control", "GDM"),
           lineheight = "auto",
           ci.vertices = TRUE,
           lwd.zero = FALSE,
           title = "TRANSCRIPTOME STUDIES",
           mean  = cbind(c(NA,as.numeric(participants_BMI$Mc)), c(NA,as.numeric(participants_BMI$Me))), 
           lower = cbind(c(NA,as.numeric(participants_BMI$lower_Control)), c(NA,as.numeric(participants_BMI$lower_GDM))),
           upper = cbind(c(NA,as.numeric(participants_BMI$upper_Control)), c(NA,as.numeric(participants_BMI$upper_GDM))), 
           new_page = TRUE,
           clip = c(15,45), 
           xlab = "                             BMI (kg/m2)     ",
           boxsize = .1,
           graphwidth = unit(20, 'cm'),
           txt_gp = fpTxtGp(ticks=gpar(cex=1), legend = gpar(cex=1), xlab = gpar(cex=1) ),
           col = fpColors(box = c("blue", "darkred")))

## --- Participant characteristics statistical test between studies --- ##
chars <- read_excel("Data/Transcriptome_Participant_characteristics.xlsx")
chars <- t(chars)
chars <- data.frame(chars)
colnames(chars) <- chars[1,] 
chars <- chars[-1,]
chars$group <- c(1,2,3,4)
chars$group2 <- c("Binder AM, 2015",	"Radaelli T, 2003",	"Ding R, 2018",	"Zhao YH, 2011")

boxplot(as.numeric(chars$BMI..kg.m2.) ~ chars$group2,main='BMI according to study',xlab='Study', ylab='BMI')

## -- ONE-WAY ANOVA TESTS -- ##
res.anova1 <- ind.oneway.second(m = as.numeric(chars$`BMI (kg/m2)`), sd = as.numeric(chars$BMI.sd), n = as.numeric(chars$N))
res.anova1
res.anova2.1 <-  ind.oneway.second(m = as.numeric(chars$Maternal.age[c(1,3,4)]), sd = as.numeric(chars$Maternal.age.sd[c(1,3,4)]), n = as.numeric(chars$N[c(1,3,4)]))
res.anova2.1
res.anova3.1 <- tsum.test(mean.x = as.numeric(chars$Maternal.age[1]), s.x = as.numeric(chars$Maternal.age.sd[1]), n.x = as.numeric(chars$N[1]),
                          mean.y = as.numeric(chars$Maternal.age[3]), s.y = as.numeric(chars$Maternal.age.sd[3]), n.y = as.numeric(chars$N[3]) )
res.anova3.1
res.anova4 <- ind.oneway.second(m = as.numeric(chars$`Gestational age (weeks)`), sd = as.numeric(chars$GA.sd), n = as.numeric(chars$N))
res.anova4


message("+-------------------------------------------------------------------------------+")
message("+              METHYLOME-WIDE STUDIES: DMGs                                 +")
message("+-------------------------------------------------------------------------------+")

## --- Overlay of DMGs --- ##

data_studies2 <- read_excel("Data/DMGs_Studies_full_list.xlsx")

Study1 <- data_studies2[1:1066,1]
colnames(Study1) <- "Genes"
Study2 <- data_studies2[,2]
colnames(Study2) <- "Genes"
Study3 <- data_studies2[1:429,3]
colnames(Study3) <- "Genes"

## Create matrix to define presence of gene in a study
all_genes2 <- join_all(list(Study1, Study2, Study3),  by = 'Genes', type = 'full')
all_genes2 <- na.omit(unique(all_genes2))
key2 <- unlist(unique(c(all_genes2)))
studies2 <- data.frame(
  key = key,
  do.call(
    cbind,
    lapply(
      data_studies2[,1:3],
      function(x) ifelse(key %in% x == TRUE, 1, 0))))
## Convert to numeric matrix
studies2 <- unique(studies2)
new_data2 <- as.matrix(sapply( studies2[,2:ncol(studies2)], as.numeric))
rownames(new_data2) <- studies2$key

## -- Getting most matched DMGs, those that are present in the 3 studies -- ##
int <- rowSums(new_data2) == 3
most_matched_genes2 <- new_data2[int,]
nrow(most_matched_genes2)

## ------- Forestplot ------- ##
participants_BMI <- read_excel("Data/BMI_full_studies_Methylation.xlsx")
tabletext2 <- cbind(c("Study",participants_BMI$Study), c("Number of control samples", participants_BMI$Nc), c("Number of GDM samples", participants_BMI$Ne))
forestplot(labeltext = tabletext2, 
           fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI),
           #legend_args = fpLegend(pos = list(x = .2, y = 0.25), 
           #gp = gpar(col = "#CCCCCC", fill = "#F9F9F9")),
           legend = c("Control", "GDM"),
           lineheight = "auto",
           #lwd.xaxis = 26.8,
           ci.vertices = TRUE,
           lwd.zero = FALSE,
           title = "METHYLOME STUDIES",
           mean  = cbind(c(NA,as.numeric(participants_BMI$Mc)), c(NA,as.numeric(participants_BMI$Me))), 
           lower = cbind(c(NA,as.numeric(participants_BMI$lower_Control)), c(NA,as.numeric(participants_BMI$lower_GDM))),
           upper = cbind(c(NA,as.numeric(participants_BMI$upper_Control)), c(NA,as.numeric(participants_BMI$upper_GDM))), 
           new_page = TRUE,
           clip = c(15,45), 
           xlab = "                             BMI (kg/m2)     ",
           txt_gp = fpTxtGp(ticks=gpar(cex=1), legend = gpar(cex=1), xlab = gpar(cex=1) ),
           boxsize = .1,
           graphwidth = unit(20, 'cm'),
           col = fpColors(box = c("blue", "darkred")))

## --- Reactome results --- ##
reactome <- read.csv("Results/Reactome_results_DMGs.csv")
significant_results.reactome <- reactome[reactome$Entities.pValue < 0.05,]
ggplot(data = significant_results.reactome, aes(x = -log10(Entities.pValue), y = Pathway.name)) + geom_bar(position = "dodge", fill = "#3C5488B2", stat = "identity") + theme_minimal(base_size = 2) + xlab("-log10(P-value)") + ylab("Reactome pathways") +
  geom_label(aes(label = Submitted.entities.found), nudge_y = 0, nudge_x = 0.005, angle =-90, size = 6) +   theme(text  = element_text(size=20) )                                                                                                                                                          

## ------ Participant characteristics statistical test between studies ------ ##
chars <- read_excel("Data/BMI_full_studies_Transcriptome.xlsx")
chars <- t(chars)
chars <- data.frame(chars)
colnames(chars) <- chars[1,] 
chars <- chars[-1,]
chars$group <- c(1,2,3)
chars$group2 <- c("Binder AM, 2015", "Finer S, 2015",	"Awamleh Z, 2021"	)
colnames(chars)[1] <- "N"

## -- ONE-WAY ANOVA TESTS -- ##
res.anova1 <- ind.oneway.second(m = as.numeric(chars$BMI..kg.m2.), sd = as.numeric(chars$BMI.sd), n = as.numeric(chars$N))
res.anova1
res.anova2 <- ind.oneway.second(m = as.numeric(chars$Maternal.age), sd = as.numeric(chars$Maternal.age.sd), n = as.numeric(chars$N))
res.anova2
res.anova3 <- ind.oneway.second(m = as.numeric(chars$Birth.weight..kg.), sd = as.numeric(chars$BMI.sd), n = as.numeric(chars$N))
res.anova3
res.anova4 <- ind.oneway.second(m = as.numeric(chars$Gestational.age..weeks.), sd = as.numeric(chars$Gestational.age..weeks.sd), n = as.numeric(chars$N))
res.anova4


## --- Overlap with single cell --- ##
## Upload data
sc <- read_excel("Data/scSeq-Results.xlsx")
meth_genes <- read_excel("Results/Common_DMGs.xlsx")
colnames(meth_genes) <- "genes"
## Separate single-cell types from data frame 
CYT1 <- sc[,1]
colnames(CYT1) <- "genes"
CYT2 <-sc[,2]
colnames(CYT2) <- "genes"
CYT3 <- sc[,3]
colnames(CYT3) <- "genes"
DC <-sc[,4]
colnames(DC) <- "genes"
EVT <- sc[,5]
colnames(EVT) <- "genes"
SYN <- sc[,6]
colnames(SYN) <- "genes"
ESF <- sc[,7]
colnames(ESF) <- "genes"
DEC <-  sc[,8]
colnames(DEC) <- "genes"

## Genes found in each cell type
c1 <- merge(CYT1, meth_genes, by = "genes")
c2 <- merge(CYT2, meth_genes, by = "genes")
c3 <- merge(CYT3, meth_genes, by = "genes")
c4 <- merge(DC, meth_genes, by = "genes")
c5 <- merge(DEC, meth_genes, by = "genes")
c6 <- merge(ESF, meth_genes, by = "genes")
c7 <- merge(EVT, meth_genes, by = "genes")
c8 <- merge(SYN, meth_genes, by = "genes")



message("+-------------------------------------------------------------------------------+")
message("+             OVERLAY OF DEGs/DMGs with in-house dataset                        +")
message("+-------------------------------------------------------------------------------+")

