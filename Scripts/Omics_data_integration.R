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

## ---- Directory ---- ##

## --- Additional functions/expressions used in the code --- ##
`%nin%` = Negate(`%in%`)


message("+-------------------------------------------------------------------------------+")
message("+           INTEGRATION OF TRANSCRIPTOME-WIDE STUDIES                           +")
message("+-------------------------------------------------------------------------------+")

## ----------------------------- DATASET 1: GSE70493 ------------------------------- ##
## ----- Samples table ----- ##
sample_names <- read.delim("text.txt")
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
## PLS-DA WITH ALL SAMPLES
group.dataset1 <- sample_names.S1$Description
placenta_gene_expression.plsda1 <- opls(counts.dataset1, factor(group.dataset1), predI = 2) 
plot(placenta_gene_expression.plsda1,
     typeVc = "x-score",
     parAsColFcVn = group.dataset1,
     parLabVc = as.character(paste0(sample_names$matName)),
     parPaletteVc = c("darkblue", "darkgreen"))

## PLS-DA WITH SELECTED SAMPLES
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
## Upload data as a normal matrix
rawData2 <- as.matrix(read.delim("Datasets/GSE19649_non-normalized.txt"))

## Samples
samples <- read.delim("Samples_GSE19649.txt")

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
data2_filtered.common <- read.csv("data2_filtered.common.csv")
rownames(data2_filtered.common) <- data2_filtered.common$X
data2_filtered.common <- data2_filtered.common[,-1]
data1_filtered.common <- read.csv("data1_filtered.common.csv")
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
write.csv2(combined.data, "Combined_Data1_Data2_matrix.csv")
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
degs.data_int <- read.csv2("DEGs_data_integration.csv")
comb_data <- read.csv2("Combined_Data1_Data2_matrix.csv")
samples <- read.csv2("Final_Samples_Data1&2.csv")
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
sc_mat <- read_excel("scSeq-Results.xlsx")
markers <- read_excel("Markers.xlsx")
deg <- read.csv2("DEGs_data_integration.csv")
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
BP <- read.delim("GO_Biological_Process_2021_table_final.txt")
BP <- BP[BP$P.value < 0.05,]
BP <- BP[order(str_count(string = BP$Genes, pattern = ";"), decreasing = TRUE),]
BP$label <- paste0("n = ", str_count(string = BP$Genes, pattern = ";")+1)
p <- ggplot(data = BP[1:20,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "#99CC99") + theme_minimal()+ xlab("-log10(P-value)") + ylab("Biological process") +
  theme(text =  element_text(size=20)) + geom_label(aes(label = label),hjust = -0.5, size = 4,
                                                    position = position_dodge(width = 1),
                                                    inherit.aes = TRUE)

## MOLECULAR FUNCTION
MF <- read.delim("GO_Molecular_Function_2021_table-2_final.txt")
MF <- MF[MF$P.value < 0.05,]
MF <- MF[order(str_count(string = MF$Genes, pattern = ";"), decreasing = TRUE),]
MF$label <- paste0("n = ", str_count(string = MF$Genes, pattern = ";")+1)
ggplot(data = MF[1:20,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "#6699CC") + theme_minimal()+ xlab("-log10(P-value)") + ylab("Molecular function") +
  theme(text =  element_text(size=20)) + geom_label(aes(label = label),hjust = -0.5, size = 4,
                                                    position = position_dodge(width = 1),
                                                    inherit.aes = TRUE)
## CELLULAR COMPONENT
CC <- read.delim("GO_Cellular_Component_2021_table-2_final.txt")
CC <- CC[CC$P.value < 0.05,]
CC <- CC[order(str_count(string = CC$Genes, pattern = ";"), decreasing = TRUE),]
CC$label <- paste0("n = ", str_count(string = CC$Genes, pattern = ";")+1)

ggplot(data = CC[1:14,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "salmon2") + theme_minimal()+ xlab("-log10(P-value)") + ylab("Cellular component") +
  theme(text =  element_text(size=30)) + geom_label(aes(label = label),hjust = -0.5, size = 6,
                                                    position = position_dodge(width = 1),
                                                    inherit.aes = TRUE)
## KEGG 
kegg <- read.delim("KEGG_2021_Human_table-2.txt")
kegg <- kegg[kegg$P.value < 0.05,]
kegg <- kegg[order(str_count(string = kegg$Genes, pattern = ";"), decreasing = TRUE),]
kegg$label <- paste0("n = ", str_count(string = kegg$Genes, pattern = ";")+1)

ggplot(data = kegg[1:7,], aes(x = -log10(P.value), y = Term)) + geom_col(fill = "pink3") + theme_minimal()+ xlab("-log10(P-value)") + ylab("KEGG pathways") +
  theme(text =  element_text(size=20)) + geom_label(aes(label = label),hjust = -0.5, size = 4,
                                                    position = position_dodge(width = 1),
                                                    inherit.aes = TRUE)
