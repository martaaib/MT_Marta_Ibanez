## R version 4.0.4 (2021-02-15)
##
## Omics data (Gene expression and DNA methylation) integration and differential analysis
##
##
##  Copyright: Marta Ibañez Lligoña (marta.ib@live.com)
##

message("+-------------------------------------------------------------------------------+")
message("+       Load the packages and basic settings                                    +")
message("+-------------------------------------------------------------------------------+")
## ----- Packages ----- ##
library(hta20transcriptcluster.db)
library(limma)
library(oligo)
library(biomaRt)
library(lmm2met)
library(ggplot2)
library("genefilter")
library("multtest")
library("annotate")
library("xtable")
library("gplots")
require(GEOquery)
library(parallel)
library("AnnotationDbi")
library(affycoretools)
library(geneExpressionFromGEO)
library(plyr)
library("edgeR")
library("limma")
library("beadarray")
library(ggplot2)
library(sva)
library(ggpubr)
library(factoextra)
library(pheatmap)
library(readxl)
library(ropls)
library(mixOmics)
library(methylumi)
library(clusterProfiler)
library(illuminaHumanv2.db)
library(affy)
library(methylumi)
library(illuminaio)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library("wateRmelon")
require(DESeq2)
library(RColorBrewer)
library(cowplot)
library(umap)
library(stringr)
library(ggVennDiagram)
library(splitstackshape)

## ---- Directory ---- ##
setwd("~/Documents/GitHub/MT_Marta_Ibanez/")
## --- Additional functions/expressions used in the code --- ##
`%nin%` = Negate(`%in%`)


message("+-------------------------------------------------------------------------------+")
message("+           INTEGRATION OF TRANSCRIPTOME-WIDE STUDIES                           +")
message("+-------------------------------------------------------------------------------+")

## ----------------------------- DATASET 1: GSE70493 ------------------------------- ##
## ----- Samples table ----- ##
sample_names <- read.delim("Data/samples_GSE70493.txt")
sample_names.S1 <- sample_names[,-1]


##  ----- Downloaded files from GEO: CEL format ----- ##
CELfiles <- list.celfiles("Datasets/")
pathdir <- "Datasets/"
rawData <- read.celfiles(file.path(pathdir,CELfiles))
sampleNames <- c(CELfiles)
sampleNames <-  data.frame(sapply(sampleNames, substr, start=1, stop=10))
colnames(sampleNames) <-"Sample"
sample_names.S1 <- merge(sample_names.S1, sampleNames, by = "Sample")
sample_names.S1$color <- "#99CCFF"
sample_names.S1$color[sample_names.S1$Description == "GDM"] <- "#FFCC99"
## Intensity of raw data
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, names = sample_names.S1$matName, col = sample_names.S1$color)

## --- Converting into txt file to compute correlation --- ##
rawData1.3 <- oligo::rma(rawData)
boxplot(rawData1.3, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, names=sample_names.S1$matName, col = sample_names.S1$color )

# Inspect the eset object
head(exprs(rawData1.3))
head(featureNames(rawData1.3))
probeids <- featureNames(rawData1.3)
nrow(exprs(rawData1.3))

## -- Filtering -- ##
data_to_filter <- rawData1.3
annotation(data_to_filter) <- "org.Hs.eg.db"
eset_filtered <- nsFilter(data_to_filter, var.func=IQR,
                          var.cutoff=0.25, var.filter=TRUE,
                          filterByQuantile=TRUE)$eset
#NUMBER OF GENES OUT
print(nrow(exprs(data_to_filter)))
print(nrow(exprs(eset_filtered)))

# -- Annotation: BiomaRt, probe ids to ENSEMBL -- ##
annotation(eset_filtered) <- "pd.hta.2.0"
eset.main <- getMainProbes(input = eset_filtered, level = "core")
probeids <- featureNames(eset.main)
eset <- exprs(eset.main)
affy_id <- sub(x = probeids, pattern =  "[.]\\d$", replacement = "")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl",mart=mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hta_2_0",
    "ensembl_gene_id",
    "external_gene_name"),
  filter = "affy_hta_2_0",
  values = affy_id, uniqueRows=TRUE)

annotLookup$affy_hta_2_0 <- paste0( annotLookup$affy_hta_2_0, ".1")
indicesLookup <- match(rownames(eset), annotLookup$affy_hta_2_0) 
new_row.names <- data.frame(rownames(eset), annotLookup[indicesLookup,])
table(rownames(eset) == annotLookup[indicesLookup,'affy_hta_2_0'])
rownames(eset) <- new_row.names$ensembl_gene_id
head(eset)
Non.annotated1 <- which(is.na(rownames(eset)) == TRUE)

## ----- PCA ----- ##
## function to compute Principal Component Analysis (PCA) and plot it
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10),ylim=c(min(pcX$x[,2])-10, max(pcX$x[,2])+10))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(eset, labels= sample_names.S1$matName, dataDesc="Normalized data",
        formapunts=c(rep(16,63)), myCex=0.8, colors = sample_names.S1$color)


## ----- PLS-DA ----- ##
## Function to select the desired number of most variable genes
select_top_variable.genes <- function(counts, TOPNUM){
  sel <- order(apply(counts, 1, var), decreasing=TRUE)[1:TOPNUM]
  counts2 <- counts[sel,]
  counts_t <- as.matrix(t(counts2))
}

counts.dataset1 <- select_top_variable.genes(eset, 2000 )
rownames(counts.dataset1) <- sample_names.S1$Sample

## --- PLS-DA WITH ALL SAMPLES --- ##
group.dataset1 <- sample_names.S1$Description
placenta_gene_expression.plsda1 <- opls(counts.dataset1, factor(group.dataset1), predI = 2) 
plot(placenta_gene_expression.plsda1,
     typeVc = "x-score",
     parAsColFcVn = group.dataset1,
     parLabVc = as.character(paste0(sample_names$matName)),
     parPaletteVc = c("darkblue", "darkgreen"))

## --- PLS-DA WITH SELECTED SAMPLES --- ##
## select and remove unwanted samples
colnames_toberemoved <- sample_names.S1[sample_names.S1$matName %in%  c("normal_R53", "normal_R53_rep1", "normal_R53_rep2", "normal_R3", "normal_R36", "GDM_R4", "normal_R5", "normal_R33",  "GDM_R6", "normal_R1", "normal_R1_rep1", "normal_R1_rep2", "GDM_R2", "GDM_R6", "GDM_R13", "normal_R7", "GDM_R26"),1]
samples1.removed1 <- sample_names.S1[sample_names.S1$matName %nin% c("normal_R53", "normal_R53_rep1", "normal_R53_rep2", "normal_R3", "normal_R36", "GDM_R4", "normal_R5", "normal_R33",  "GDM_R6", "normal_R1", "normal_R1_rep1", "normal_R1_rep2", "GDM_R2", "GDM_R6", "GDM_R13", "normal_R7", "GDM_R26"),]
counts.dataset1.2 <- counts.dataset1[rownames(counts.dataset1) %nin% colnames_toberemoved,]
group.dataset1.2 <- samples1.removed1$Description

## perform PLS-DA with ropls package
placenta_gene_expression.plsda2 <- opls(counts.dataset1.2, factor(group.dataset1.2), predI = 2) 
plot(placenta_gene_expression.plsda2,
     typeVc = "x-score",
     parAsColFcVn = group.dataset1.2,
     parLabVc = as.character(paste0(samples1.removed1$matName)),
     parPaletteVc = c("darkblue", "darkgreen"))

## PCA with selected samples 
plotPCA(t(counts.dataset1.2), labels= sample_names.S1$matName, dataDesc="Normalized data",
        formapunts=c(rep(16,63)), myCex=0.8, colors = sample_names.S1$color)

## --------------------------- DATASET 2: GSE19649 --------------------------- ##
## ----- Input matrix from GEO ----- ##
x<-getGEO("GSE19649")
rawData2.expression <- x[[1]]
ids.data2 <- featureData(rawData2.expression)@data
## Collect ids from dataframe
pobeids2 <- featureNames(rawData2.expression)
## Upload data as a normal matrix (needs to be downloaded from GEO-NCBI)
rawData2 <- as.matrix(read.delim("Datasets/GSE19649_non-normalized.txt"))

## Samples
samples <- read.delim("Data/Samples_GSE19649.txt")

## Select samples of GEO matrix from placenta
rawData2.1 <- as.matrix(data.frame(rawData2[,1], rawData2[,8], rawData2[,9],rawData2[,10], rawData2[,11]))
colnames(rawData2.1) <- c("ID_REF", "GSM490136","Detection.Pval.3", "GSM490137" ,"Detection.Pval.4")
rownames(rawData2.1) <- rawData2.1[,1]
rownames.rawData21 <- rownames(rawData2.1)
rawData2.1 <- rawData2.1[,-1]
rawData2.1 <- apply(rawData2.1,2, function(x) as.numeric(as.character(x)))
data2 <- as.matrix(data.frame(log2(rawData2.1)[,1], rawData2.1[,2], log2(rawData2.1)[,3], rawData2.1[,4]))
colnames(data2) <- c("GSM490136","Detection.Pval.3", "GSM490137" ,"Detection.Pval.4")
rownames(data2) <- rownames(rawData2.1)

## intensity of log2 intensity data
boxplot(data2[,c(1,3)], which="all",las=2, main="Intensity distribution of log2 intensity data- not normalized", 
        cex.axis=0.6, names = c("GSM490136", "GSM490137"), col = c("#99CCCC", "#66CC99"))

## Quantile normalization ##
data2.intensities <- as.matrix(data.frame(data2[,1], data2[,3]))
colnames(data2.intensities) <- c("GSM490136", "GSM490137")
p_values <- as.matrix(data2[,c(2,4)])
data2.raw <- as.matrix(data.frame(rawData2.1[,1], rawData2.1[,3]))
data2.raw <- apply(data2.raw,2, function(x) as.numeric(as.character(x)))
colnames(data2.raw) <- c("GSM490136", "GSM490137")
data2.NORM <- neqc(data2.raw, detection.p = p_values)
## boxplot 
boxplot(data2.NORM, which="all",las=2, main="Intensity distribution of normalized data", 
        cex.axis=0.6, names = c("GSM490136", "GSM490137"), col = c("#99CCCC", "#66CC99"))

## Filtering ##
data2.norm <- as.matrix(data.frame(data2.NORM[,1], p_values[,1], data2.NORM[,2], p_values[,2]))
colnames(data2.norm) <- c( "GSM490136","Detection.Pval.3", "GSM490137" ,"Detection.Pval.4")
rownames(data2.norm) <- rownames.rawData21
## Remove genes with threshold
data2_filtered <- data2.norm[data2.norm[,2] < 0.05 & data2.norm[,4] < 0.05,]
data2_filtered <- data2_filtered[sum(data2_filtered[,1], data2_filtered[,3]) >= 1,]

## Annotation ##
## PART 1: Change to platform id ##
probesid.filtered <- data.frame(rownames(data2_filtered))
colnames(probesid.filtered) <- "Search_Key"
platform <- read.delim("GPL6102-11574.txt")
illid2 <- merge(platform, probesid.filtered, by = "Search_Key")
illid2Lookup <- match(rownames(data2_filtered), illid2$Search_Key) 
new_row.names2 <- data.frame(rownames(data2_filtered), illid2[illid2Lookup,])
table(rownames(data2_filtered) == illid2[illid2Lookup,'Search_Key'])
illid3 <- new_row.names2$ID
rownames(data2_filtered) <- illid3

## PART 2: BIOMART TO ENSEMBL ##
annotLookup2 <- getBM(
  mart=mart,
  attributes=c(
    "illumina_humanwg_6_v2",
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene_id"),
  filter = "illumina_humanwg_6_v2",
  values = illid3, uniqueRows=TRUE)

illid3Lookup <- match(rownames(data2_filtered), annotLookup2$illumina_humanwg_6_v2) 
new_row.names3 <- data.frame(rownames(data2_filtered), annotLookup2[illid3Lookup,])
table(rownames(data2_filtered) == annotLookup2[illid3Lookup,'illumina_humanwg_6_v2'])
rownames(data2_filtered) <- new_row.names3$ensembl_gene_id

## -------- Correlation between common transcripts from the two datasets (before cross-platform normalization) ------- ##
transcripts.data1 <- data.frame(rownames(eset))
head(transcripts.data1)
colnames(transcripts.data1) <- "ID"
transcripts.data2 <- data.frame(rownames(data2_filtered))
colnames(transcripts.data2) <- "ID"
head(transcripts.data2)
common.transcripts <- merge(transcripts.data1, transcripts.data2, by = "ID")
common.transcripts <- na.omit(common.transcripts)

## Select data1 and data 2 transcripts transcripts
data2_filtered.common <- data2_filtered[common.transcripts$ID,]
data2_filtered.common <- data2_filtered.common[,-c(2,4)]
data1_filtered.common <- eset[common.transcripts$ID,]

## Remove possible duplicated rows
data1_filtered.common <- unique(data1_filtered.common)
all(duplicated(rownames(data1_filtered.common)))
data2_filtered.common <- unique(data2_filtered.common)
all(duplicated(rownames(data2_filtered.common)))

## Check order of transcripts of both datasets
all(rownames(data2_filtered.common) == rownames(data1_filtered.common))

## Save filtered datasets
write.csv(data2_filtered.common, "data2_filtered.common.csv")
write.csv(data1_filtered.common, "data1_filtered.common.csv")

## ------------ Combine Data -------------- ##
## --- Load data (if necessary) --- ##
data2_filtered.common <- read.csv("Data/data2_filtered.common_TRANSCRIPTOME.csv")
rownames(data2_filtered.common) <- data2_filtered.common$X
data2_filtered.common <- data2_filtered.common[,-1]
data1_filtered.common <- read.csv("Data/data1_filtered.common_TRANSCRIPTOME.csv")
rownames(data1_filtered.common) <- data1_filtered.common$X
data1_filtered.common <- data1_filtered.common[,-1]

## --- Combine samples data --- ##
sample_names.S1 <- sample_names.S1[sample_names.S1$Sample %in% colnames(data1_filtered.common),]
colnames(data1_filtered.common) <- sample_names.S1$Sample
data1_filtered.common <- data1_filtered.common[,colnames(data1_filtered.common) %in% sample_names.S1$Sample]
combined.data <- as.matrix(data.frame(data1_filtered.common,data2_filtered.common ))
#combined_samples <- sample_names.S1
combined_samples <-sample_names.S1
combined_samples <- combined_samples[,-4]
#samples$color <- c("#FFCC99", "#99CCFF")
samples$matName <- c( "GDM_R01", "Normal_R02")
colnames(samples) <- c("Sample", "Description",    "matName" )
combined_samples <- rbind(combined_samples, samples)
combined_samples$batch <- 1
combined_samples$batch[nrow(combined_samples)] <- 2
combined_samples$batch[nrow(combined_samples)-1] <- 2
combined_samples$color <- "#FFCC99"
combined_samples$color[combined_samples$Description == "Normal"] <- "#99CCFF"

## Remove outliers
combined_samples <- combined_samples[-which((combined_samples$matName %in% c("normal_R37", "normal_R32", "GDM_R25", "GDM_R21", "GDM_R45", "GDM_R43", "GDM_R30", "GDM_R19", "GDM_R10", "GDM_R24", "GDM_R39", "normal_R42", "normal_R23", "GDM_R43", "GDM_R11", "GDM_R46", "GDM_R34", "GDM_R41", "normal_R52_rep1", "normal_R52_rep2", "GDM_55_rep1", "GDM_R55_rep2", "GDM_R49", "normal_R38", "normal_R40", "GDM_R54"))),]
combined.data <- combined.data[,(colnames(combined.data) %in% combined_samples$Sample)]

## Save combined samples and matrix of log2 intensity values
write.csv2(combined.data, "Combined_Data1_Data2_matrix_transcriptome.csv")
write.csv2(combined_samples, "Final_Samples_Data1&2.csv")

## Cross-platform normalization and PCA
# PCA of combined data
plotPCA(combined.data, labels= combined_samples$matName, dataDesc="Combined data",
        formapunts=c(rep(16,63)), myCex=0.8, color = combined_samples$color)
groups <- as.factor(combined_samples$Description)
groups <- relevel(groups, "Normal")
mod = model.matrix(~groups, data=combined_samples)
# Cross-platform normalization
combat_edata3 = ComBat(dat=combined.data, batch = combined_samples$batch, mod = mod )
## PCA after normalization
plotPCA(combat_edata3, labels= combined_samples$matName, dataDesc="Combined data",
        formapunts=c(rep(16,63)), myCex=0.8, color = combined_samples$color)

## Compute and plot PCA through ggplot ##
pca1 <- prcomp(t(combat_edata3))
pca1.plot <- data.frame(pca1$x)
pca1.plot$label <- combined_samples$matName
pca1.plot$group <- combined_samples$Description
pca1.plot$dataset <- combined_samples$batch
eigs <- pca1$sdev^2
percentVar <- c(round(100*(eigs[1] / sum(eigs))), round(100*(eigs[2] / sum(eigs))))
pca1.plot$dataset[pca1.plot$dataset == 1] <- "Binder AM, 2015"
pca1.plot$dataset[pca1.plot$dataset == 2] <- "Zhao YH, 2011"
ggplot(pca1.plot, aes(x = PC1, y = PC2)) + geom_point(aes(col = group, shape = as.factor(dataset), size = 0)) + theme_minimal() + labs(col = "Group", shape = "Study") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + theme_bw() +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme(legend.text=element_text(size=15), axis.title = element_text(size=15))

## Pearson correlation of combined data
Pear.corr_crossnorm <- cor(combat_edata3, combat_edata3, method = "pearson")
summary(Pear.corr_crossnorm)

## ---- PLS-DA of combined data ---- ##
combat_edata3.T <- select_top_variable.genes(combat_edata3, nrow(combat_edata3))
group.dataset.combined <- groups
placenta_gene_expression.plsda1 <- opls(combat_edata3.T, factor(group.dataset.combined), predI = 2) 
plot(placenta_gene_expression.plsda1,
     typeVc = "x-score",
     parAsColFcVn = group.dataset.combined,
     parLabVc = as.character(paste0(combined_samples$matName)),
     parPaletteVc = c("red", "blue"))
## Save VIP scores of genes
key_genes <- data.frame(getVipVn(placenta_gene_expression.plsda1))

### ------------ Differential expression analysis ------------- ###

groups <- as.factor(combined_samples$Description)
groups <- relevel(groups, "Normal")
batch <- factor(c(rep("batch1", nrow(combined_samples)-2), rep("batch2", 2)))
design <-  model.matrix(object = ~ groups + batch)
colnames(design) <- c("Normal", "GDM", "batch")
rownames(design) <- combined_samples$matName
fit <- eBayes(lmFit(combined.data, design))
y <- topTable(fit, coef=2, number=Inf)
DEgenes <- y[y$adj.P.Val<0.05,]

## ----- PLS-DA ------ ##
vip <- data.frame(key_genes[match(rownames(key_genes), rownames(y)),])
colnames(vip) <- "VIP"
resLFC.VIP <- cbind(y, vip)
pls.DEGs <- resLFC.VIP[resLFC.VIP$adj.P.Val < 0.05 & resLFC.VIP$VIP> 1, ]

## Annotation through ENSEMBL ##
message("+-------------------------------------------------------------------------------+")
message("+    Connect to ENSEMBL database through Biomart package                        +")
message("+-------------------------------------------------------------------------------+")

ensembl        <- useMart("ensembl") 
ensembl        <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

ensEMBL2id     <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'chromosome_name'), 
                        mart = ensembl, useCache = FALSE)  

chr.rmInd      <- c(grep("CHR*", ensEMBL2id$chromosome), 
                    grep("GL*", ensEMBL2id$chromosome),
                    grep("KI*", ensEMBL2id$chromosome))
ensEMBL2id     <- ensEMBL2id[-chr.rmInd,]

head(ensEMBL2id)

## --- DEGs: Annotation of gene name and description --- ##
id <-  match( rownames(pls.DEGs) ,ensEMBL2id$ensembl_gene_id)
na_vals <- rownames(pls.DEGs)
null_pos <- which(is.na(id))
pls.DEGs$external_geneid <- ensEMBL2id$external_gene_name[id]
pls.DEGs$external_geneid[null_pos] <- na_vals[null_pos]
pls.DEGs$description <- ensEMBL2id$description[id]
no_name <-pls.DEGs[null_pos,]

## Save DEGs in file
write.csv2(pls.DEGs,"DEGs_data_integration.csv")
write.csv2()

## Heatmap
degs.data_int <- read.csv2("Results/DEGs_data_integration.csv")
comb_data <- read.csv2("Data/Combined_Data1_Data2_matrix_transcriptome.csv")
samples <- read.csv2("Final_Samples_Data1&2_transcriptome.csv")
mod = model.matrix(~relevel(as.factor(Description), "Normal"), data=samples)
colnames(mod) <- c("Normal", "GDM")
combat_edata3 = ComBat(dat=comb_data[,-1], batch = samples$batch, mod = mod )
rownames(combat_edata3) <- comb_data$X
heatmap.data <- combat_edata3[degs.data_int$X,]
test <- data.frame(colnames(heatmap.data),samples$Sample) ## same order
colnames(heatmap.data) <- samples$matName
heatmap.data <- heatmap.data[,order(colnames(heatmap.data))] 
logCPM <- t(scale(t(heatmap.data)))
pheatmap(logCPM, scale = "row", cluster_cols = FALSE, labels_row = rep(" ", nrow(heatmap.data)))

## Integration of single-cell, UMAP dimensionality reduction and generation of heatmap
sc_mat <- read_excel("Data/scSeq-Results_mat.xlsx")
deg <- read.csv2("Results/DEGs_data_integration.csv")
colnames(deg)[9] <- "Gene"
head(sc_mat)
colnames(sc_mat) <- sc_mat[3,]
sc_mat <- sc_mat[-c(1:3),]
found <- merge(sc_mat, deg, by = "Gene")
found <- found[,colnames(found) != "SE"]
rownames(found) <- found$Gene
found <- found[,-1]
data <- found[,1:8] 
data<- apply(data, 2, function(x) as.numeric(as.character(x)))
rownames(data) <- rownames(found)
data <- data[rowSums(data) > 0,]
try <- umap(data)
df <- try$layout
colnames(df) <-  c("UMAP_1", "UMAP_2")
groups <- data.frame(apply(data, 1, max))
groups$cell <- colnames(data)[apply(data,1,which.max)]
df <- cbind(df, data.frame(rownames(df)), type = "UMAP",celltype =  groups$cell)

ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.7, alpha = 0.5, aes(col = celltype))  + theme_minimal() + guides(col = guide_legend("Cell type"))

degs2 <- data.frame(rownames(df))
pheatmap(data, scale = "row",  labels_row = rep(" ", nrow(data)), cluster_cols = FALSE, color = hcl.colors(20, "BlueRed2", rev = FALSE))

## ------------- Plots for enrichment analysis --------------- ##
## BIOLOGICAL PROCESS
BP <- read.delim("Results/GO_Biological_Process_2021_table_TRANSCRIPTOME.txt")
BP <- BP[BP$P.value < 0.05,]
BP <- BP[order(str_count(string = BP$Genes, pattern = ";"), decreasing = TRUE),]
BP$label <- paste0("n = ", str_count(string = BP$Genes, pattern = ";")+1)
p <- ggplot(data = BP[1:20,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "#99CC99") + theme_minimal()+ xlab("-log10(P-value)") + ylab("Biological process") +
  theme(text =  element_text(size=20)) + geom_label(aes(label = label),hjust = -0.5, size = 4,
                                                    position = position_dodge(width = 1),
                                                    inherit.aes = TRUE)

## MOLECULAR FUNCTION
MF <- read.delim("Results/GO_Molecular_Function_2021_table_TRANSCRIPTOME.txt")
MF <- MF[MF$P.value < 0.05,]
MF <- MF[order(str_count(string = MF$Genes, pattern = ";"), decreasing = TRUE),]
MF$label <- paste0("n = ", str_count(string = MF$Genes, pattern = ";")+1)
ggplot(data = MF[1:20,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "#6699CC") + theme_minimal()+ xlab("-log10(P-value)") + ylab("Molecular function") +
  theme(text =  element_text(size=20)) + geom_label(aes(label = label),hjust = -0.5, size = 4,
                                                    position = position_dodge(width = 1),
                                                    inherit.aes = TRUE)
## CELLULAR COMPONENT
CC <- read.delim("Results/GO_Cellular_Component_2021_table_TRANSCRIPTOME.txt")
CC <- CC[CC$P.value < 0.05,]
CC <- CC[order(str_count(string = CC$Genes, pattern = ";"), decreasing = TRUE),]
CC$label <- paste0("n = ", str_count(string = CC$Genes, pattern = ";")+1)

ggplot(data = CC[1:14,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "salmon2") + theme_minimal()+ xlab("-log10(P-value)") + ylab("Cellular component") +
  theme(text =  element_text(size=30)) + geom_label(aes(label = label),hjust = -0.5, size = 6,
                                                    position = position_dodge(width = 1),
                                                    inherit.aes = TRUE)
## KEGG 
kegg <- read.delim("Results/KEGG_2021_Human_table_TRANSCRIPTOME.txt")
kegg <- kegg[kegg$P.value < 0.05,]
kegg <- kegg[order(str_count(string = kegg$Genes, pattern = ";"), decreasing = TRUE),]
kegg$label <- paste0("n = ", str_count(string = kegg$Genes, pattern = ";")+1)

ggplot(data = kegg[1:7,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "pink3") + theme_minimal()+ xlab("-log10(P-value)") + ylab("KEGG pathways") +
  theme(text =  element_text(size=20)) + geom_label(aes(label = label),hjust = -0.5, size = 4,
                                                    position = position_dodge(width = 1),
                                                    inherit.aes = TRUE)



message("+-------------------------------------------------------------------------------+")
message("+              INTEGRATION OF METHYLOMEE-WIDE STUDIES                           +")
message("+-------------------------------------------------------------------------------+")

## ---------------------- DATASET 1: GSE70453 ---------------------- ##
## LOAD DATA, STORE INFORMATION AND P-VALUES
## Intensity matrix
untar("GSE70453_RAW.tar", exdir = "GSE70453_RAW/idat")
Study1 <- read.csv("GSE70453_Matrix_signal_intensities.csv")
Study1 <- as.matrix(Study1)
rownames(Study1) <- Study1[,1]
Study1.ids <- Study1[,1]
Study1  <- Study1[,-1]
rownames(Study1)
head(Study1)

## Separate matrix
p.values.stored <- seq(from = 3, to = ncol(Study1), by = 3)
p.values <- as.matrix(Study1[, p.values.stored])
probes <- rownames(p.values)
mn.stored <- seq(from = 2, to = ncol(Study1), by = 3)
un.stored <- seq(from = 1, to = ncol(Study1), by = 3)
mn <- as.matrix(Study1[, mn.stored])
un <- as.matrix(Study1[,un.stored ])

## Load into minfi object: MethylSet
data1 <- readGEORawFile(filename ="GSE70453_Matrix_signal_intensities.csv", Uname = "Unmethylated.Signal", Mname = "Methylated.signal",  )

## Samples
GSE70453_samples <- read.delim("Data/GSE70453_samples.txt")
files <- sampleNames(data1)
GSE70453_samples <- cbind(GSE70453_samples, files)

## Load samples in minfi object
GSE70453_samples$Group <- "GDM"
normal <- grep("Normal", GSE70453_samples$MatName)
GSE70453_samples$Group[normal] <- "Normal"

# visualise what the data looks like before and after normalisation
densityPlot(getBeta(data1), sampGroups=GSE70453_samples$Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(GSE70453_samples$Group)), 
       text.col=brewer.pal(8,"Dark2"))


## Filtering of samples
p.values <- apply(p.values, 2, function(x) as.numeric(as.character(x)))
rownames(p.values) <- probes
head(p.values)
keep <- colMeans(p.values) < 0.05
data1.samples_filtered <- data1[,keep] ## all significant samples, none removed
data1.samples_filtered

## Normalization
data1.samples_filtered.norm <- preprocessQuantile(data1.samples_filtered)
densityPlot(getBeta(data1.samples_filtered.norm), sampGroups=GSE70453_samples$Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(GSE70453_samples$Group)), 
       text.col=brewer.pal(8,"Dark2"))

# ensure probes are in the same order in the mSetSq and detP objects
p.values <- p.values[match(featureNames(data1.samples_filtered.norm),rownames(p.values)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(p.values < 0.05) == ncol(data1.samples_filtered.norm) 
table(keep)
data1.samples_filtered.norm <- data1.samples_filtered.norm[keep,]
data1.samples_filtered.norm ## as this dataset only contains women, we do not need to remove any sex-specific probes

# remove probes with SNPs at CpG site
data1.samples_filtered.norm <- dropLociWithSnps(data1.samples_filtered.norm)
data1.samples_filtered.norm

# exclude cross reactive probes 
xReactiveProbes <- read.csv(file="Data/48639-non-specific-probes-Illumina450k.csv")
head(xReactiveProbes)
keep <- !(featureNames(data1.samples_filtered.norm) %in% xReactiveProbes$TargetID)
table(keep)
data1.filtered.norm <- data1.samples_filtered.norm[keep,] 
data1.filtered.norm

## Save normalized  M values 
norm.Mvals <- getM(data1.filtered.norm)
head(norm.Mvals)

## Save normalized  BETA values 
norm.Bvals <- getBeta(data1.filtered.norm)
head(norm.Bvals)
colnames(norm.Bvals) <- GSE70453_samples$MatName

## Principal Component Analysis (PCA) of Beta values
## set colors
GSE70453_samples$color <- "#0072B2"
GSE70453_samples$color[grep("Normal", GSE70453_samples$Group)] <- "#D55E00"
## with labels
pca1 <- prcomp(t(getBeta(data1.filtered.norm)))
pca1 <- pca1$x
rownames(pca1) <- GSE70453_samples$MatName
pca1 <- pca1[,1:2]
batch1 <- pca1[pca1[,2] < 0,]
batch1 <- rownames(batch1)
batch2 <-pca1[pca1[,2] > 0,]
batch2 <- rownames(batch2)
GSE70453_samples$batch <- NA
GSE70453_samples$batch[which(GSE70453_samples$MatName %in% batch1)] <- 1
GSE70453_samples$batch[which(GSE70453_samples$MatName %in% batch2)] <- 2

## Remove samples
toremove <- c("R22", "R19", "R20", "R13", "R66", "R59", "R26")
GSE70453_samples_r <- filter(GSE70453_samples, grepl(paste(toremove, collapse="|"), GSE70453_samples$MatName))
GSE70453_samples <- GSE70453_samples[GSE70453_samples$MatName %nin% GSE70453_samples_r$MatName,]
norm.Bvals <- norm.Bvals[,which(colnames(norm.Bvals) %nin% GSE70453_samples_r$MatName)]

## PCA plot
plotPCA(norm.Bvals, labels= GSE70453_samples$MatName, dataDesc="Normalized data",
        formapunts=c(rep(16,82)), myCex=0.8, colors = GSE70453_samples$color)

## M-values
## PCA
plotPCA(getM(data1.filtered.norm),  dataDesc="Normalized data",
        formapunts=c(rep(16,82)), myCex=0.8, colors = GSE70453_samples$color)



## ---------------------- DATASET 2:GSE153220  ---------------------- ##
untar("GSE153220_RAW.tar", exdir = "GSE153220/idat")
## RUN: list idat files and unzip
idatFiles <- list.files("GSE153220/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

## Annotation
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

# Samples
GSE153220_samples <- read.delim("Data/GSE153220_samples.txt")
GSE153220_samples <- cbind(GSE153220_samples, sampleNames(data2))
GSE153220_samples$Group <- "GDM"
normal <- grep("Normal", GSE153220_samples$MatName)
GSE153220_samples$Group[normal] <- "Normal"

## Load idat files
data2 <- read.metharray.exp("GSE153220/idat")
data2
sampleNames(data2) <- GSE153220_samples$MatName

## Quality control
detP <- detectionP(data2)
head(detP)
pal <- brewer.pal(8,"Dark2")
barplot(colMeans(detP), col=pal[factor(GSE153220_samples$Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(GSE153220_samples$Group)), fill=pal,
       bg="white")
qcReport(data2, sampGroups=GSE153220_samples$Group, 
         pdf="qcReport.pdf")

# remove poor quality samples
keep <- colMeans(detP) < 0.05
table(keep) ## No sample needs to be removed


## Normalization
mSetSq2 <- preprocessFunnorm(data2)

# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(data2, sampGroups=GSE153220_samples$Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(GSE153220_samples$Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq2), sampGroups=GSE153220_samples$Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(GSE153220_samples$Group)), 
       text.col=brewer.pal(8,"Dark2"))


## Filtering probes
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq2),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.05) == ncol(mSetSq2) 
table(keep)
mSetSqFlt2 <- mSetSq2[keep,]
mSetSqFlt2 

# remove probes with common SNPs which may affect  CpG sites/ we will be removing the default ones, but you can also remove all of those with minor allele frequencies
mSetSqFlt2 <- dropLociWithSnps(mSetSqFlt2)
mSetSqFlt2


## We can also filter all of those  of those probes that are know to be cross-reactive, probes that have been demonstrated to map to multiple places in the genome
## Chen et al. 2013 
keep <- !(featureNames(mSetSqFlt2) %in% xReactiveProbes$TargetID)
table(keep)
mSetSqFlt2 <- mSetSqFlt2[keep,] 
mSetSqFlt2

## MDS
pal <- brewer.pal(8,"Dark2")
plotMDS(getM(mSetSqFlt2), top=1000, gene.selection="common", 
        col=pal[factor(GSE153220_samples$Group)], cex=0.8)
legend("right", legend=levels(factor(GSE153220_samples$Group)), text.col=pal,
       cex=0.65, bg="white")

## PCA
GSE153220_samples$color <- "#0072B2"
GSE153220_samples$color[normal] <- "#D55E00"

## Without batch effects removed
## Beta values
plotPCA(getBeta(mSetSqFlt2),  dataDesc="Normalized data",
        formapunts=c(rep(16,36)), myCex=0.8, colors = GSE153220_samples$color, labels = GSE153220_samples$MatName)

## Save normalized M values
norm.Mvals2 <- getM(mSetSqFlt2)
hist(norm.Mvals2)
head(norm.Mvals2)
norm.Mvals2.2=norm.Mvals2[apply(norm.Mvals2,1,var)>0,]


## Save normalized BETA values
norm.Bvals2 <- getBeta(mSetSqFlt2)
head(norm.Bvals2)

## NORMALIZED INTENSITY VALUES
mSetSqFlt2

## --- Common M and Beta values --- ##
Mvals.data1 <- data.frame(rownames(norm.Mvals))
head(Mvals.data1)
colnames(Mvals.data1) <- "ID"
Mvals.data2 <- data.frame(rownames(norm.Mvals2))
head(Mvals.data2)
colnames(Mvals.data2) <- "ID"
head(Mvals.data2)
common.Mvals <- merge(Mvals.data1, Mvals.data2, by = "ID")
common.Mvals <- na.omit(common.Mvals)

## Select data1 and data 2 common M and B values
Mvals1.common <- norm.Mvals[common.Mvals$ID,]
Mvals2.common <- as.matrix(norm.Mvals2[common.Mvals$ID,])

Bvals1.common <- norm.Bvals[common.Mvals$ID,]
Bvals2.common <- norm.Bvals2[common.Mvals$ID,]

## Check order 
all(rownames(Mvals1.common) == rownames(Mvals2.common))
all(rownames(Bvals1.common) == rownames(Bvals2.common))

## ------ Pearson correlarion ------- ##
## Using normalized Beta  values
joined_Bvals <- cbind(Bvals1.common, Bvals2.common)
joined_samples <- rbind(GSE153220_samples[,c(1:3,6)], GSE70453_samples[,c(1:3,7)])
joined_samples$color<- NA

## Color by condition
joined_samples$color[grep("GDM", joined_samples$MatName)] <- "#0072B2"
joined_samples$color[grep("Normal", joined_samples$MatName)] <- "#D55E00"
joined_samples$color2 <- NA

## Color by study
joined_samples$color2[grep("Genomic", joined_samples$Description)] <- "#0072B2"
joined_samples$color2[grep("Placenta", joined_samples$Description)] <- "#D55E00"

## Without batch effects removed
Meth.corr2 <- cor(Bvals1.common, Bvals2.common, method = "pearson")
summary(Meth.corr2) 

## With batch effects removed

## PCA
plotPCA(joined_Bvals, labels= joined_samples$SampleFile, dataDesc="Normalized data",
        formapunts=c(rep(16,63)), myCex=0.8, colors = joined_samples$color)

## PCA plot with ggplot
joined_samples$Group <- "GDM"
joined_samples$Group[grep("Normal", joined_samples$MatName)] <- "Normal"
pca2 <- prcomp(t(joined_Bvals))
pca.plot <- as.data.frame(pca2$x)
eigs <- pca2$sdev^2
percentVar <- c(round(100*(eigs[1] / sum(eigs))), round(100*(eigs[2] / sum(eigs))))
pca.plot$group <- joined_samples$Group
pca.plot$study <- joined_samples$batch
pca.plot$study[pca.plot$study %in% c(1,2)] <- "Binder AM, 2015"
pca.plot$study[pca.plot$study %in% c(3,4)] <- "Awamleh Z, 2021"

ggplot(pca.plot, aes(PC1, PC2), show.legend = TRUE) +
  geom_point(size=2, alpha=0.75, aes(shape = study, col = group)) +
  #egeom_encircle(aes(group = group, fill = NULL)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(col = "Group", shape = "Study", size = 3) +
  theme_bw() + theme(legend.text=element_text(size=15), axis.title = element_text(size=15), legend.title = element_text(size=15))

## Remove selected unwanred samples
colnames(joined_Bvals) <- joined_samples$SampleFile

## Samoples that will be removed
tobrem <- c("R58", "R75", "R51", "R50", "R80", "R68", "R47", "R61", "R62", "R64", "R49", "R70", "R60", "R73", "R63", "R77", "R67", "R65", "R56",
            "R48", "R53", "R71", "R72", "R66", "R69", "R78", "R82", "R81", "R79", "R52", "R76", "R55", "R54", "R74", "R57", "R59", "R038", "R033", "R030", "R070", "R029", "R45", "R46")
toberem2 <- c( "R8", "R020", "R11", "R7", "R10", "R16", "R013", "R023", "R42", "R15", "R17", "R002", "R031", "R035","R045",
               "R025", "R037", "R028", "R44", "R024", "R27", "R043", "R041", "R044", "R018", "R042", "R048", "R37", "R048", "R1")
tobrem <- c(tobrem, toberem2)

## Samples to keep
tokeep <- c("GSM4636037",	"GSM4636038",	"GSM4636039",	"GSM4636040",	"GSM4636041",	"GSM4636043",	"GSM4636044",	"GSM4636047",	"GSM4636048"	,"GSM4636052",	"GSM4636058",	"GSM4636062",	"GSM4636068",	"GSM4636069",	"GSM4636071",	"GSM1755028",	"GSM1755029",	"GSM1755030",	"GSM1755031",	"GSM1755032",	"GSM1755035",	
            "GSM1755047",	"GSM1755049",	"GSM1755050",	"GSM1755051",	"GSM1755054",	"GSM1755055",	"GSM1755056",	"GSM1755057",	"GSM1755058",	"GSM1755059",	"GSM1755060",	"GSM1755061"	,"GSM1755062",	"GSM1755064",	"GSM1755065",	"GSM1755066",	"GSM1755067")

joined_samples2 <- filter(joined_samples, grepl(paste(tobrem, collapse="|"), joined_samples$MatName))
joined_samples.flt <- joined_samples[joined_samples$MatName %nin% joined_samples2$MatName,]
joined_samples.flt <- joined_samples[joined_samples$SampleFile %in% tokeep,]
colnames(joined_Bvals) <- joined_samples$SampleFile

### LOAD DATA FROM B-VALUES OF DESIRED SAMPLES ##
joined_Bvals <- read.csv2("Data/Bvals_final.csv")
head(joined_Bvals)
joined_Bvals.flt <- joined_Bvals[,-1]
rownames(joined_Bvals.flt) <- joined_Bvals$X
joined_samples.flt <- read.csv2("Data/final_samples_methylation.csv")

## Beta
joined_Bvals.flt <- joined_Bvals[,colnames(joined_Bvals) %in%  joined_samples.flt$SampleFile]

## PCA OF COMBINED SAMPLES
plotPCA(joined_Bvals.flt, labels= joined_samples.flt$MatName, dataDesc="Normalized data",
        formapunts=c(rep(16,63)), myCex=0.8, colors = joined_samples.flt$color)

joined_samples$Group <- "GDM"
joined_samples$Group[grep("Normal", joined_samples$MatName)] <- "Normal"

## Filtered
joined_samples.flt$Group <- "GDM"
joined_samples.flt$Group[grep("Normal", joined_samples.flt$MatName)] <- "Normal"

## final PCA plot with ggplot ##
pca3 <- prcomp(t(joined_Bvals.flt))
pca3.plot <- as.data.frame(pca3$x)
eigs <- pca3$sdev^2
percentVar <- c(round(100*(eigs[1] / sum(eigs))), round(100*(eigs[2] / sum(eigs))))
pca3.plot$group <- joined_samples.flt$Group
pca3.plot$study <- joined_samples.flt$batch
pca3.plot$study[pca3.plot$study %in% c(1,2)] <- "Binder AM, 2015"
pca3.plot$study[pca3.plot$study %in% c(3,4)] <- "Awamleh Z, 2021"

ggplot(pca3.plot, aes(PC1, PC2), show.legend = TRUE) +
  geom_point(size=2, alpha=0.75, aes(shape = study, col = group)) +
  #geom_encircle(aes(group = group, col = group)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(col = "Group", shape = "Study", size = 3) +
  theme_bw() + theme(legend.text=element_text(size=15), axis.title = element_text(size=15), legend.title = element_text(size=15))


## --------- Differential methylation analysis --------- ##
groups <- as.factor(joined_samples.flt$Group)
groups <- relevel(groups, "Normal")
design <-  model.matrix(object = ~ groups)
colnames(design) <- c("Normal", "GDM")
rownames(design) <- joined_samples.flt$MatName
fit <- eBayes(lmFit(joined_Bvals.flt, design))
y <- topTable(fit, coef=2, number=Inf)
DEgenes <- y[y$adj.P.Val<0.05,]

## Annotation of regions ##
ann450kSub <- ann450k[match(rownames(joined_Bvals.flt),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit, num= 7023, coef=2, genelist = ann450kSub)
head(DMPs)
write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

DMPs.B <- match(rownames(DMPs), rownames(joined_Bvals.flt))
DMPs.B <- joined_Bvals.flt[DMPs.B,]
all(rownames(DMPs.B) == rownames(DMPs))

DMPs.Bvals <- cbind(DMPs, DMPs.B)
write.table(DMPs.Bvals, file="DMPs_Bvals.csv", sep=",", row.names=FALSE)


## Keep names of genes and save into file ##
genes <- data.frame(DMPs$UCSC_RefGene_Name)
genes <- genes[- which(genes$DMPs.UCSC_RefGene_Name == ""),]
genes <- data.frame(unlist(strsplit(genes, ";")))
genes <- unique(genes)
write.csv(genes, "methylated_genes.csv")
write.csv2(joined_samples.flt,"final_samples_meth.csv")
write.csv2(joined_Bvals.flt, "Bvals_final.csv")

## Heatmap ##
DMGs <- read.csv("DMPs_Bvals.csv")
samples <- read.csv2("Data/final_samples_meth.csv")
heatmap.data <- DMGs[, 29:ncol(DMGs)]
test <- data.frame(colnames(heatmap.data),samples$SampleFile) ## same order
colnames(heatmap.data) <- samples$MatName
heatmap.data <- heatmap.data[,order(colnames(heatmap.data))] 
logCPM <- t(scale(t(heatmap.data)))
pheatmap(heatmap.data, scale = "row", cluster_cols = TRUE, labels_row = rep(" ", nrow(heatmap.data)), color = hcl.colors(20, "Cork", rev = FALSE))
pheatmap(logCPM, scale = "row", cluster_cols = FALSE, labels_row = rep(" ", nrow(heatmap.data)))

## ENRICHMENT ANALYSIS
## BIOLOGICAL PROCESS
bp <- read.delim("Results/GO_Biological_Process_2021_table_METHYLATION.txt")
bp <- bp[bp$P.value < 0.05,]
bp <- bp[order(str_count(string = bp$Genes, pattern = ";"), decreasing = TRUE),]
bp$label <- paste0("n = ", str_count(bp$Genes, ";")+1)
p <-ggplot(data = bp[1:20,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "#99CC99") + theme_minimal()+ xlab("-log10(P-value)") + ylab("Biological process") +
  theme(text =  element_text(size=20)) + geom_label(aes(label = label),hjust = -0.5, size = 6,
                                                    position = position_dodge(width = 1),
                                                    inherit.aes = TRUE)

## MOLECULAR FUNCTION
mf <- read.delim("Results/GO_Molecular_Function_2021_table_METHYLATION.txt")
mf <- mf[mf$P.value < 0.05,]
mf <- mf[order(str_count(string = mf$Genes, pattern = ";"), decreasing = TRUE),]
mf$label <- paste0("n = ", str_count(mf$Genes, ";")+1)
p <- ggplot(data = mf[1:20,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "#6699CC") + theme_minimal()+ xlab("-log10(P-value)") + ylab("Molecular function") + 
  geom_label(aes(label = label),hjust = -0.5, size = 6,
             position = position_dodge(width = 1),
             inherit.aes = TRUE) + theme(text =  element_text(size=20))

## CELLULAR COMPONENT
cc <- read.delim("Results/GO_Cellular_Component_2021_table_METHYLATION.txt")
cc <- cc[cc$P.value < 0.05,]
cc <- cc[order(str_count(string = cc$Genes, pattern = ";"), decreasing = TRUE),]
cc$label <- paste0("n = ", str_count(cc$Genes, ";")+1)
p <- ggplot(data = cc, aes(x = -log10(P.value), y = Term)) + geom_col(fill = "salmon2") + theme_minimal()+ xlab("-log10(P-value)") + 
  ylab("Cellular component") + theme(text =  element_text(size=30)) + geom_label(aes(label = label),hjust = -0.5, size = 6,
                                                                                 position = position_dodge(width = 1),
                                                                                 inherit.aes = TRUE)

## KEGG
KEGG <- read.delim("Results/KEGG_2021_Human_table_METHYLATION.txt")
KEGG <- KEGG[KEGG$P.value < 0.05,]
KEGG <- KEGG[order(str_count(string = KEGG$Genes, pattern = ";"), decreasing = TRUE),]
KEGG$label <- paste0("n = ", str_count(KEGG$Genes, ";")+1)
ggplot(data = KEGG, aes(x = -log10(P.value), y = Term)) + geom_col(fill = "pink3") + theme_minimal()+ xlab("-log10(P-value)") + 
  ylab("KEGG pathways") + theme(text =  element_text(size=20)) + geom_label(aes(label = label),hjust = -0.5, size = 4,
                                                                            position = position_dodge(width = 1),
                                                                            inherit.aes = TRUE)

## Integration with single-cell placental transcriptome + UMAP dimensionality reduction
sc_mat <- read_excel("Data/scSeq-Results_mat.xlsx")
dmg <- genes
colnames(dmg) <- "Gene"
head(dmg)
colnames(sc_mat) <- sc_mat[3,]
sc_mat <- sc_mat[-c(1:3),]
found <- merge(sc_mat, dmg, by = "Gene")
found <- found[,colnames(found) != "SE"]
rownames(found) <- found$Gene
found <- found[,-1]
data <- found[,1:8] 
data<- apply(data, 2, function(x) as.numeric(as.character(x)))
rownames(data) <- rownames(found)
data <- data[rowSums(data) > 0,]
try <- umap(data)
df <- try$layout
colnames(df) <-  c("UMAP_1", "UMAP_2")
groups <- data.frame(apply(data, 1, max))
groups$cell <- colnames(data)[apply(data,1,which.max)]
df <- cbind(df, data.frame(rownames(df)), type = "UMAP",celltype =  groups$cell)

## plot UMAP
ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.7, alpha = 0.5, aes(col = celltype))  + theme_minimal() + labs(col = "Cell type") + theme(legend.text=element_text(size=10), legend.title  = element_text(size=10))
degs2 <- data.frame(rownames(df))

## PLOT HEATMAP OF SINGLE-CELL MATRIX CORRESPONDING GENES
cols <- colorRampPalette(brewer.pal(5, "Blues"))
pos_df = data.frame("Pos" = df_used$Pos) ## add source
rownames(pos_df) = rownames(df_num)
pheatmap(data, scale = "row",  labels_row = rep(" ", nrow(data)), cluster_cols = FALSE, color =  hcl.colors(20, "BlueRed2", rev = FALSE)) 


## OVERLAP WITH QUANTSEQ DATA + PUBLIC DATA AND CREATE VENN DIAGRAMS ##
public_data.degs <- read_excel("Data/DEGs_Studies_full_list.xlsx")
public_data.dmgs <-read_excel("Data/DMGs_full_studies.xlsx")
DEGs.DMGs_int <- read_excel("Results/DEGs_DMGs_integration.xlsx") ## DEGs and DMGs from integration of transcriptome and methylome-wide datasets
QS.DEGs <- rread.csv2("pls_DEGs_GDMvsLEAN.csv") ## in-house DEGs
qs <- data.frame(QS.DEGs$external_geneid)
qs <- qs[qs$QS.DEGs.external_geneid != "",]
colnames(public_data.degs) <- c("genes", "genes", "genes", "genes")
colnames(public_data.dmgs) <- c("genes", "genes", "genes")
public.DEGs <- na.omit(unique(rbind(data.frame(genes = public_data.degs[,1]), data.frame(genes = public_data.degs[,2]), data.frame(genes = public_data.degs[,3]), data.frame(genes = public_data.degs[,4]))))
public.DEGs <- apply(public.DEGs,2, function(x) as.character(x))
public.DMGs <- apply(public.DMGs,2, function(x) as.character(x))
public.DMGs <- na.omit(unique(rbind(data.frame(genes = public_data.dmgs[,1]), data.frame(genes = public_data.dmgs[,2]), data.frame(genes = public_data.dmgs[,3]))))
DEGs_int <- na.omit(DEGs.DMGs_int[,1])
DEGs_int <- apply(DEGs_int,2, function(x) as.character(x))
DMGs_int <- na.omit(DEGs.DMGs_int[,2])
DMGs_int <- apply(DMGs_int,2, function(x) as.character(x))
data.venn <- list(
  public.DEGs = c(public.DEGs$Gene),
  public.DMGs = c(public.DMGs),
  DEGs_int =c(DEGs_int ),
  DMGs_int = c(DMGs_int),
  QS = c(qs))
ggVennDiagram(data.venn)
qs <- data.frame(qs)
DEGs.venn <- list(
  public.DEGs = c(public.DEGs$Gene),
  DEGs_int =c(DEGs_int$Gene ),
  QS = c(qs$Gene))
t <- ggVennDiagram(DEGs.venn)
colnames(public.DEGs) <- "Gene"
colnames(DEGs_int) <- "Gene"
colnames(qs) <- "Gene"

public.DEGs <- data.frame(public.DEGs)
DEGs_int <- data.frame(DEGs_int)
m <- join_all(list(public.DEGs,DEGs_int, qs), type = "inner", by = "Gene")

## Venn diagram of public DMGs and DMGs from integration of methylome-wide results
DMGs.venn <- list(
  public.DMGs = c(public.DMGs),
  DMGs_int = c(DMGs_int))
ggVennDiagram(DMGs.venn)
write.table(public.DEGs,"t")
data.studies <-matrix(NA, ncol = 5)

## Venn diagram of integration results: common DEGs and DMGs 
int.venn <- list(
  DEGs_int =c(DEGs_int$Gene ),
  DMGs_int = c(DMGs_int))
ggVennDiagram(int.venn)

## ------ Correlate gene expression with methylation ------- ##
expression.vals <- read.csv2("Data/Combined_Data1_Data2_matrix_transcriptome.csv")
degs.data_int <- read.csv2("Results/DEGs_data_integration.csv")
samples <- read.csv2("Data/Final_Samples_Data1&2_transcriptome.csv")
mod = model.matrix(~relevel(as.factor(Description), "Normal"), data=samples)
colnames(mod) <- c("Normal", "GDM")
combat_edata3 = ComBat(dat=expression.vals[,-1], batch = samples$batch, mod = mod )
rownames(combat_edata3) <- expression.vals$X
cor.data1 <- combat_edata3[which(rownames(combat_edata3) %in% degs.data_int$X),]
DMGs1 <- data.frame(Gene = DMGs_int)
colnames(DMGs1) <- "Gene"
DEGS1 <- data.frame(Gene = degs.data_int$external_geneid)
colnames(DEGS1) <- "Gene"
ovlp <- merge(DMGs1,DEGS1 , by = "Gene") 
ens <- degs.data_int$X[degs.data_int$external_geneid %in% c(ovlp$Gene)]
cor.data1 <- data.frame(combat_edata3[ens,])
rownames(cor.data1) <- ovlp$Gene
cor.data1 <- cor.data1[order(rownames(cor.data1)),]
DMGs <- read.csv("Results/DMPs_Bvals.csv")
out <- cSplit(DMGs, "UCSC_RefGene_Name", sep=";", "long")
DMGs <- out[out$UCSC_RefGene_Name  %in% c(ovlp$Gene) ,]
DMGs$mean <- rowMeans(DMGs[,29:ncol(DMGs)])
cor.data2 <- aggregate(mean ~ UCSC_RefGene_Name, DMGs, mean)

cor.data1.2 <- data.frame(rowMeans(cor.data1))
cor.data2.2 <- data.frame(cor.data2$mean)
rownames(cor.data2.2) <- cor.data2$UCSC_RefGene_Name
cor.data2.2 <- cor.data2.2[order(rownames(cor.data2.2)),]
t <- cbind(cor.data1.2, cor.data2.2)

correlation <- cor.test(matrix(t(cor.data1.2)), matrix(t(cor.data2.2)), type = "pearson")





