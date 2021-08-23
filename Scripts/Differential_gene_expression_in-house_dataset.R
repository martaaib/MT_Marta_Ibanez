## R version 4.0.4 (2021-02-15)
##
## Rscript for differential expression analysis of in-house dataset
##
##
##  Copyright: Marta Ibañez Lligoña (marta.ib@live.com)
##

message("+-------------------------------------------------------------------------------+")
message("+       Load the packages and basic settings                                    +")
message("+-------------------------------------------------------------------------------+")

## Packages
library(ropls)
library(edgeR)
library(readxl)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(mygene)
library("org.Hs.eg.db")
library(pheatmap)
library(mygene)
library("genefilter")
require(plyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggalt)
library(reshape2)
library("ggVennDiagram")
`%nin%` = Negate(`%in%`)

## Working directory
base.dir <- "~/Documents/GitHub/MT_Marta_Ibanez/"
setwd(base.dir)
TOPNUM <- 2000

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

message("+-------------------------------------------------------------------------------+")
message("+                    Load samples data and counts matrix                        +")
message("+-------------------------------------------------------------------------------+")

counts.matrix <- read.table("CTR_ans48_0009-QuantSeq_rmUMIs_PCRdup_count_matrix.txt", header = TRUE, sep = " ")
samples <- read.csv2("samples_information.csv")

## --- collapse replicates --- ##
colnames(counts.matrix)
colnames(counts.matrix) <- sapply(colnames(counts.matrix), gsub, pattern = "_L001", replacement = "")
colnames(counts.matrix) <- sapply(colnames(counts.matrix), gsub, pattern = "_L002", replacement = "")
counts.matrix <- as.matrix(counts.matrix)
counts.matrix_collapse <- sumTechReps(counts.matrix)

## --- Sorting and re-naming matrices --- ##
samples <- samples[order(samples$Sample),]
test <- data.frame(colnames(counts.matrix_collapse), samples$Sample) ## same order
colnames(counts.matrix_collapse) <- samples$Sample


message("+-------------------------------------------------------------------------------+")
message("+                 PCA: Top 2000 variable genes                                  +")
message("+-------------------------------------------------------------------------------+")


## ----- GDM vs LEAN SAMPLES ----- ##
samples.GDMvsLEAN <- samples[-grep("obese", samples$matName),]
## remove outlier
samples.GDMvsLEAN <- samples.GDMvsLEAN[-which(samples.GDMvsLEAN$Sample == "40633_110"),]

## DEFINITIVE SAMPLES: Samples which need to be removed
rmv <- c("lean_R6", "lean_R11", "lean_R10", "GDM_R8", "GDM_R7", "GDM_R4", "GDM_R6")
samples.GDMvsLEAN <- samples.GDMvsLEAN[-which(samples.GDMvsLEAN$matName %in% rmv),]

## Remove samples from counts and filtering of low count genes and normalizarion for PCA and PLS-DA
counts.matrix_collapseGDMvsLEAN <- counts.matrix_collapse[, colnames(counts.matrix_collapse) %in% samples.GDMvsLEAN$Sample]
dgList <- DGEList(counts = counts.matrix_collapseGDMvsLEAN, samples = samples.GDMvsLEAN, group = samples.GDMvsLEAN$Group)
apply(dgList$counts, 2, sum) # total gene counts per sample
dgList <- calcNormFactors(dgList, method = "TMM")
keep <- rowSums(cpm(dgList) >=5) >= 8
dgList <- dgList[keep,]

#counts.matrix_collapseGDMvsLEAN_cpm <- cpm(dgList$counts, log=TRUE, prior.count=5) 
counts.matrix_collapseGDMvsLEAN_cpm <- cpm(dgList)

## --- selecting top 2000 genes with most variance --- ##
sel <- order(apply(counts.matrix_collapseGDMvsLEAN_cpm, 1, var), decreasing=TRUE)[1:2000]

## ----- PCA ----- ##
pca <- prcomp(t(counts.matrix_collapseGDMvsLEAN_cpm[sel,]), scale. = TRUE)
eigs <- pca$sdev^2
percentVar.star <- c(round(100*(eigs[1] / sum(eigs))), round(100*(eigs[2] / sum(eigs))))
pca_data <- as.data.frame(pca$x)
pca_data[,ncol(pca_data)+1] <- samples.GDMvsLEAN$Group
colnames(pca_data)[ncol(pca_data)] <- "Group"
pca_data[,ncol(pca_data)+1] <- samples.GDMvsLEAN$Sample
colnames(pca_data)[ncol(pca_data)] <- "Label"
pca_data <- cbind(pca_data, samples.GDMvsLEAN$RIN)
pca_data$RIN <- NA
pca_data$RIN[pca_data$`samples.GDMvsLEAN$RIN` < 3] <- "< 3"
pca_data$RIN[pca_data$`samples.GDMvsLEAN$RIN` > 3] <- "> 3"
pca_data$Group[pca_data$Group == "lean"] <- "Normal"
## ----- PLOT PCA ----- ##
pdf(file = "~/Documents/GitHub/QuantSeq-data-analysis/Results/PCA/PCA_GDMvsLEAN_finalsamples.pdf")
ggplot(pca_data, aes(PC1, PC2, color=Group, fill =Group, label=Label), show.legend = TRUE) +
  geom_point(size=2, alpha=0.75, aes(shape = RIN), show.legend = TRUE ) +
  geom_text_repel(aes(label=Label), show.legend = FALSE, size=2, max.overlaps = 40) +
  geom_encircle(aes(group = Group, fill = NULL)) +
  xlab(paste0("PC1: ",percentVar.star[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.star[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes)") ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='top', aspect.ratio=1) 
dev.off()


message("+-------------------------------------------------------------------------------+")
message("+                                 PLS-DA                                        +")
message("+-------------------------------------------------------------------------------+")

## Function to select top variable genes/ order genes by variance + transposing of matrix for PLS-DA
select_top_variable.genes <- function(counts, TOPNUM){
  sel <- order(apply(counts, 1, var), decreasing=TRUE)[1:TOPNUM]
  counts2 <- counts[sel,]
  counts_t <- as.matrix(t(counts2))
}


## ----- GDM vs LEAN SAMPLES ----- ##
## Select group
samples.GDMvsLEAN$Group[samples.GDMvsLEAN$Group == "lean"] <- "Normal"
group.GDMvsLean <- samples.GDMvsLEAN$Group

## Transpose matrix
counts.matrix_collapseGDMvsLEAN.T <- select_top_variable.genes(counts.matrix_collapseGDMvsLEAN_cpm, nrow(counts.matrix_collapseGDMvsLEAN_cpm))

## --- PLS-DA function from ropls package --- ##
placenta_gene_expression.plsda <- opls(counts.matrix_collapseGDMvsLEAN.T, factor(group.GDMvsLean), predI = 2) 
## predI = 2 --> It means that the first predictive  component was not meaningful and checks for any additional ones.

## plot pls-da 
pdf(file = "~/Documents/GitHub/QuantSeq-data-analysis/Results/PLS-DA/PLSDA_GDMvsLEAN.pdf")
plot(placenta_gene_expression.plsda,
     typeVc = "x-score",
     parAsColFcVn = group.GDMvsLean,
     parLabVc = as.character(paste0(samples.GDMvsLEAN$Sample)),
     parPaletteVc = c("darkblue", "darkgreen"))
dev.off()

## Save VIPs from PLS ##
key_genes <- data.frame(getVipVn(placenta_gene_expression.plsda))


message("+-------------------------------------------------------------------------------+")
message("+               edgeR: Differential expression analysis                         +")
message("+-------------------------------------------------------------------------------+")


## --- Design matrix --- ##
groups <- factor(group.GDMvsLean)
groups <- relevel(groups, "lean")
design <- model.matrix(object = ~ groups)
colnames(design)<-c("LEAN", "GDM")
rownames(design)<- samples.GDMvsLEAN$sampleNames
print(design); 

## --- DISPERSION --- ##
dgList <- estimateCommonDisp(dgList, design=design)
dgList <- estimateTagwiseDisp(dgList, design=design)

## --- DEGs --- ##
fit <- glmQLFit(dgList,design = design)
lrt <- glmLRT(fit, coef = 2)
summary(de <- decideTestsDGE(lrt, p.value = 0.05, lfc = 1))
resLFC <- topTags(lrt, n = Inf)$table

## Select genes with logFC cut-off of 1 and FDR < 0.05
edgeR.significant.GDMvsLEAN <- resLFC[which((resLFC$logFC < -1 | resLFC$logFC > 1) & resLFC$FDR < 0.05),]

## --- Annotation of gene name and description --- ##
id <-  match( rownames(edgeR.significant.GDMvsLEAN) ,ensEMBL2id$ensembl_gene_id)
na_vals <- rownames(edgeR.significant.GDMvsLEAN)
null_pos <- which(is.na(id))
edgeR.significant.GDMvsLEAN$external_geneid <- ensEMBL2id$external_gene_name[id]
edgeR.significant.GDMvsLEAN$external_geneid[null_pos] <- na_vals[null_pos]
edgeR.significant.GDMvsLEAN$description <- ensEMBL2id$description[id]
no_name <-edgeR.significant.GDMvsLEAN[null_pos,]

## --- Add normalized counts to table --- ##
logCPM.edgeR <- counts.matrix_collapseGDMvsLEAN_cpm[match(rownames(edgeR.significant.GDMvsLEAN), rownames(counts.matrix_collapseGDMvsLEAN_cpm)),]
edgeR.significant.GDMvsLEAN <- cbind(edgeR.significant.GDMvsLEAN, logCPM.edgeR)

## Save document
write.csv2(edgeR.significant.GDMvsLEAN, file = "Results/edgeR_DEGs_GDMvsLEAN.csv")


message("+-------------------------------------------------------------------------------+")
message("+               edgeR + PLS: Differential expression analysis                   +")
message("+-------------------------------------------------------------------------------+")

## match genes with DEGs
all.genes_VIP <- data.frame(key_genes[match(rownames(key_genes), rownames(resLFC)),])
colnames(all.genes_VIP) <- "VIP"
resLFC.VIP <- cbind(resLFC, all.genes_VIP)
pls.DEGs <-  resLFC.VIP[which((resLFC.VIP$logFC < -1 | resLFC.VIP$logFC > 1) & resLFC.VIP$FDR < 0.05 & resLFC.VIP$VIP > 1),]

## --- DEGs: Annotation of gene name and description --- ##
id <-  match( rownames(pls.DEGs) ,ensEMBL2id$ensembl_gene_id)
na_vals <- rownames(pls.DEGs)
null_pos <- which(is.na(id))
pls.DEGs$external_geneid <- ensEMBL2id$external_gene_name[id]
pls.DEGs$external_geneid[null_pos] <- na_vals[null_pos]
pls.DEGs$description <- ensEMBL2id$description[id]
no_name <-pls.DEGs[null_pos,]

## --- Add normalized counts to table --- ##
logCPM.PLS <- counts.matrix_collapseGDMvsLEAN_cpm[match(rownames(pls.DEGs), rownames(counts.matrix_collapseGDMvsLEAN_cpm)),]
pls.DEGs <- cbind(pls.DEGs, logCPM.PLS)

## --- All genes: Annotation and save csv file with all counts from all genes --- ##
id <-  match( rownames(resLFC.VIP) ,ensEMBL2id$ensembl_gene_id)
na_vals <- rownames(resLFC.VIP)
null_pos <- which(is.na(id))
resLFC.VIP$external_geneid <- ensEMBL2id$external_gene_name[id]
resLFC.VIP$external_geneid[null_pos] <- na_vals[null_pos]
resLFC.VIP$description <- ensEMBL2id$description[id]
no_name <-resLFC.VIP[null_pos,]

## Save document
write.csv2(pls.DEGs, file = "Results/pls_DEGs_GDMvsLEAN.csv")
write.csv2(resLFC.VIP, file = "Results/pls_allgenes_GDMvsLEAN.csv")
WriteXLS::WriteXLS(pls.DEGs, ExcelFileName = "Results/pls_DEGs_GDMvsLEAN.xlsx")


message("+-------------------------------------------------------------------------------+")
message("+       Volcano plot and heatmap of differentially expressed genes              +")
message("+-------------------------------------------------------------------------------+")

## --- Load data (if necessary) --- ##
DEGs <- read.csv2("~/Documents/GitHub/QuantSeq-data-analysis/Results/DEGs_GDMvsLEAN/pls_DEGs_GDMvsLEAN.csv")
all.genes <- read.csv2("~/Documents/GitHub/QuantSeq-data-analysis/Results/DEGs_GDMvsLEAN/pls_allgenes_GDMvsLEAN.csv")

## ---  Classical volcano plot --- ##
## Set variables for plot: Regulation and label
all.genesV <- all.genes[,1:9]
all.genesV$Regulation1 <- NA
all.genesV$Regulation1[which(all.genesV$logFC < -1 & all.genesV$FDR < 0.05)] <- "Down"
all.genesV$Regulation1[which(all.genesV$logFC > 1 & all.genesV$FDR < 0.05)] <- "Up"
all.genesV$label1 <- NA
all.genesV$label1[which(is.na(all.genesV$Regulation1) == FALSE)] <- all.genesV$external_geneid[which(is.na(all.genesV$Regulation1) == FALSE)]

## plot
pdf(file = "~/Documents/GitHub/QuantSeq-data-analysis/Results/Classical_Volcano_plot_EdgeR_GDMvsLEAN.pdf", width = 8, height = 6)
ggplot(data=all.genesV, aes(x=logFC, y=-log10(FDR), col=Regulation1, label=label1)) + 
  geom_point() +  geom_vline(xintercept = c(-1,1), col = "black", linetype = "dashed") + geom_hline(yintercept = 1.3, col = "black", linetype = "dashed") +
  theme_minimal()  + geom_label_repel(aes(label=label1), show.legend = FALSE, size=2, max.overlaps = 25) + 
  labs(title = "edgeR:  GDM vs LEAN ")  + scale_colour_manual( values = c("#3399FF","#FF3333" ), name = "Regulation")
dev.off()

## --- PLS: Adding size as VIP in volcano plot --- ##
## -- Set direction of regulation according to PLS and edgeR -- ##
all.genesV$Regulation2 <- "NA"
all.genesV$Regulation2[which(all.genesV$logFC < -1 & all.genesV$VIP > 1 & all.genesV$FDR < 0.05)] <- "Down"
all.genesV$Regulation2[which(all.genesV$logFC > 1 & all.genesV$VIP > 1 & all.genesV$FDR < 0.05)] <- "Up"

## Set size of VIP ##
all.genesV$size <- as.factor(all.genesV$VIP)
all.genesV$size <- as.numeric(as.character(all.genesV$size))
all.genesV$size[all.genesV$VIP < 1] <- 0
all.genesV$size <- as.factor(round(all.genesV$size))

## -- Set labels 2 -- ##
all.genesV$label2 <- NA
all.genesV$label2[which(all.genesV$Regulation2 != "NA")] <- all.genesV$external_geneid[which(all.genesV$Regulation2 != "NA")]

## plot with labels
pdf(file = "~/Documents/GitHub/QuantSeq-data-analysis/Results/Volcanoplot_PLS_VIP_GDMvsLEAN.pdf")
ggplot(data=all.genesV, aes(x=logFC, y=-log10(FDR), col=Regulation2, label=label2)) + 
  geom_point(aes(size = size)) + 
  theme_minimal() +  geom_vline(xintercept = c(-1,1), col = "black", linetype = "dashed" ) + geom_hline(yintercept = 1.3, col = "black", linetype = "dashed") +
  geom_label_repel(aes(label=label2), show.legend = FALSE, size=2, max.overlaps = 10) + labs(title = "edgeR + PLS: DEGs from GDM vs LEAN samples", size = "VIP")  + 
  scale_colour_manual( values = c("#3399FF","#CCCCCC","#FF3333" ), name = "Regulation") + guides(fill=guide_legend(title="VIP")) + ylim(c(-1, 25))
dev.off()

## plot without labels
pdf(file = "~/Documents/GitHub/QuantSeq-data-analysis/Results/Volcanoplot_PLS_VIP_GDMvsLEAN_without_labels.pdf")
ggplot(data=all.genesV, aes(x=logFC, y=-log10(FDR), col=Regulation2, label=label2)) + 
  geom_point(aes(size = size)) + 
  theme_minimal() +  geom_vline(xintercept = c(-1,1), col = "black", linetype = "dashed" ) + geom_hline(yintercept = 1.3, col = "black", linetype = "dashed") +
  labs(title = "edgeR + PLS: DEGs from GDM vs LEAN samples", size = "VIP")  + 
  scale_colour_manual( values = c("#3399FF","#CCCCCC","#FF3333" ), name = "Regulation") + guides(fill=guide_legend(title="VIP"))+ ylim(c(-1, 40)) + xlim(c(-5.5, 5.5))
dev.off()
## ----- HEATMAP ----- ##
DEGs.heatmap <- DEGs[,10:ncol(DEGs)]
## Avoid duplicated gene names
DEGs$external_geneid[which(DEGs$external_geneid == "Y_RNA")] <- DEGs$X[which(DEGs$external_geneid == "Y_RNA")]
DEGs$external_geneid[which(DEGs$external_geneid == '')] <- DEGs$X[which(DEGs$external_geneid == '')]
rownames(DEGs.heatmap) <- DEGs$external_geneid
DEGs.heatmap <- DEGs.heatmap[order(DEGs$FDR),]
## Change sample names
colnames(DEGs.heatmap) <- samples.GDMvsLEAN$Sample
colnames <- match(colnames(DEGs.heatmap), samples.GDMvsLEAN$Sample)
colnames(DEGs.heatmap) <- samples.GDMvsLEAN$matName[colnames]
DEGs.heatmap <- DEGs.heatmap[,order(colnames(DEGs.heatmap))] 
## calculate z-scores
logCPM <- t(scale(t(DEGs.heatmap)))
logCPM <- logCPM[order(DEGs$FDR),]
colnames(logCPM) <- c("GDM_R1",  "GDM_R10", "GDM_R11", "GDM_R12", "GDM_R13", "GDM_R14", "GDM_R3",  "GDM_R5",  "GDM_R9",  "normal_R1", "normal_R2", "normal_R3", "normal_R4", "normal_R5", "normal_R7", "normal_R8", "normal_R9")
## Pheatmap
## Heatmap for all DEGs
pheatmap(logCPM, cluster_cols = FALSE , gaps_col = 9,  fontsize_row = 15, cellwidth = 15, labels_row = rep(" ", nrow(logCPM)) )
pheatmap(DEGs.heatmap, cluster_cols = FALSE , gaps_col = 9,  fontsize_row = 15, cellwidth = 15, labels_row = rep(" ", nrow(logCPM)), scale = "row" )
## Heatmap top 50 most significant DEGs
pheatmap(logCPM[1:50,], cluster_cols = FALSE , gaps_col = 9,  fontsize_row = 15, cellwidth = 15)

## ComplexHeatmap
Heatmap(logCPM, column_km = 2, row_labels = rep(" ", nrow(logCPM)))
Heatmap(logCPM[1:50,], column_km = 3, row_labels = rep(" ", nrow(logCPM[1:50,])))

## ----- Plots of individual gene counts ----- ##
DEGs_plot <- DEGs$X
DEGs.counts <- DEGs.heatmap[DEGs$X %in% DEGs_plot,]
selected_DEGs.counts <- DEGs.counts[sample(nrow(DEGs.counts), 5), ]
selected_DEGs.counts$external_geneid <- rownames(selected_DEGs.counts)
long <- melt(selected_DEGs.counts, id = "external_geneid", variable_name = "Sample")
ggplot(data = long, aes(x = variable, y = value, fill = external_geneid)) + geom_bar(position = "dodge", stat = "identity") +
  theme_bw()

message("+-------------------------------------------------------------------------------+")
message("+       ClusterProfiler: Analysis of differentially expressed genes             +")
message("+-------------------------------------------------------------------------------+")

## ------ GDM vs LEAN ----- ##

## --- Load data (if necessary) --- ##
DEGs <- read.csv2("~/Documents/GitHub/QuantSeq-data-analysis/Results/DEGs_GDMvsLEAN/pls_DEGs_GDMvsLEAN.csv")
all.genes <- read.csv2("~/Documents/GitHub/QuantSeq-data-analysis/Results/DEGs_GDMvsLEAN/pls_allgenes_GDMvsLEAN.csv")

## Change format in order to run KEGG pathway enrichment analysis
search_kegg_organism('hsa', by='kegg_code')
EN.DEGs <- queryMany(DEGs$X, scopes="ensembl.gene", fields="entrezgene", species="human", retrunall = TRUE)
EN.DEGs <- EN.DEGs$entrezgene

## --- KEGG pathways --- ##
kk <- enrichKEGG(gene = EN.DEGs , organism = 'hsa', pvalueCutoff = 0.05, pAdjustMethod = "fdr")
kk <- setReadable(kk, org.Hs.eg.db, keyType="ENTREZID")
head(kk, n = 10)
dotplot(kk)
KEGG <- kk@result[kk@result$p.adjust < 0.05,]

## -- BP: biological process -- ##
BP <- enrichGO(gene = DEGs$external_geneid,keyType ="SYMBOL", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff =0.05)
head(BP, n = 10)
dotplot(BP)
BP_results <- BP@result[BP@result$p.adjust < 0.05,]


## -- MF: Molecular function -- ##
MF <- enrichGO(gene = DEGs$external_geneid,keyType = "SYMBOL", OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
head(MF, n = 10) 
dotplot(MF)
MF_results <- MF@result[MF@result$p.adjust < 0.05,]

## -- CC: Cellular component -- ##
CC <- enrichGO(gene =  DEGs$external_geneid ,keyType = "SYMBOL", OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05)
head(CC, n = 10)
dotplot(CC)
CC_results <- CC@result[CC@result$p.adjust < 0.05,]

## --- Compare clusters of regulation --- ##
Up <- DEGs$X[DEGs$logFC > 1]
Down <- DEGs$X[DEGs$logFC < -1]
Degs_clusters <- list(up = Up, down = Down)
## -- GO enrichment analysis by group -- ##
xx <- compareCluster(Degs_clusters, fun="enrichGO", OrgDb='org.Hs.eg.db', pvalueCutoff=0.05, keyType = 'ENSEMBL')
xx <- setReadable(xx, org.Hs.eg.db, keyType="ENSEMBL")
dotplot(xx)
comparecluster.result <- as.data.frame(xx)

## Checking direction of genes
g <- c("NSG1/LAMP1/LTV1/LITAF/GRN/CHID1/PSAP/IFITM3/ATG9A/HLA-DRA/LAMTOR2/SLC30A2/SQSTM1/CLN3/LAPTM4A/LAMTOR1/HLA-DRB1/VAMP8/SCARB2/RNF13/GIMAP5")
g <- strsplit(g, "/")
g <- c(g[[1]])
reg <- DEGs[DEGs$external_geneid %in% g,]


## Plots
BP_results$Regulation <- "Up"
BP_results$label <- strsplit(BP_results$geneID, "/")
MF_results$Regulation <- "Up"
MF_results$Regulation[c(6,8)] <- "Down"
CC_results$Regulation <- "Up"
CC_results$Regulation[37] <- "Down"
ggplot(data = BP_results, aes(x = -log10(pvalue), y = Description)) + geom_col(fill = "#99CC99") + theme_minimal() + xlab("-log10(P-value)") + ylab("GO terms: biological process")
ggplot(data = MF_results, aes(x = -log10(pvalue), y = Description, fill = Regulation)) + geom_col() + theme_minimal() + scale_fill_brewer(palette="Blues", name = "Regulation") + theme_minimal()+ xlab("-log10(P-value)") + ylab("GO terms: molecular function") + theme(text = element_text(size = 15))

ggplot(data = CC_results, aes(x = -log10(pvalue), y = Description, fill = Regulation)) + geom_col() 


ggplot(data = KEGG, aes(x = -log10(pvalue), y = Description)) + geom_col(fill = "#6699CC") + theme_minimal()+ xlab("-log10(P-value)") + ylab("KEGG pathways") + theme_minimal() + xlab("-log10(P-value)") + ylab("GO terms: cellular component")+ scale_fill_brewer(palette="Reds", name = "Regulation") + theme(text = element_text(size = 15))

## OVERLAP WITH LITERATURE (Figure 10)
qs <- data.frame(DEGs$external_geneid)
colnames(qs) <- "Genes"
## Transcriptome studies
public <- read_excel("Data/DEGs_Studies_full_list.xlsx")
pub1 <- na.omit(data.frame(public[,1]))
colnames(pub1) <- "Genes"
pub2 <- na.omit(data.frame(public[,2]))
colnames(pub2) <- "Genes"
pub3 <- na.omit(data.frame(public[,3]))
colnames(pub3) <- "Genes"
pub4 <- na.omit(data.frame(public[,4]))
colnames(pub4) <- "Genes"

all <- unique(rbind(pub1,pub2,pub3,pub4))
merge1 <- merge(qs, all, by = "Genes")

## Methylation ##
public <- read_excel("Data/DMGs_Studies_full_list.xlsx")
pub1
pub1 <- na.omit(data.frame(public[,1]))
colnames(pub1) <- "Genes"
pub2 <- na.omit(data.frame(public[,2]))
colnames(pub2) <- "Genes"
pub3 <- na.omit(data.frame(public[,3]))
colnames(pub3) <- "Genes"
all2 <- unique(rbind(pub1,pub2,pub3))
merge2 <- merge(qs, all2, by = "Genes")

x <- list(
  QuantSeq = c(qs$Genes), 
  Transcriptome = c(all$Genes), 
  Methylome = c(all2$Genes)
)
ov <- ggVennDiagram(x)
ov$labels
t <- merge(merge1, merge2, by = "Genes")
WriteXLS::WriteXLS(t, "Common_genes_in-house&published_data. xlsx")
