setwd("D:/Sharif University/7th Semester/Bioinfromatics/Project/")

library(GEOquery)
library(limma)
library(umap)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(plyr)
library(ica)
library(ca)
library(Rtsne)
library(uwot)

series <- "GSE48558"
platform <- "GPL6244"

gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group membership for all samples
gsms <- paste0("1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XXXXX",
               "XXXXXXXXXXXXXXXXXX0X0XXX0X0000X0XX00XX00X0X0X0X0X0",
               "XXX0XXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000110111",
               "00000000000000000000")
sml <- strsplit(gsms, split="")[[1]]

ex <- exprs(gset)

##### Quality Control

# 1.boxplot
max(ex)
pdf("Results/boxplot.pdf")
boxplot(ex)
dev.off()
pdf("Results/boxplot_zoomed.pdf",width=85)
boxplot(ex)
dev.off()

# 2.pheatmap
pdf("Results/pheatmap.pdf",width = 15,height = 15)
pheatmap(cor(ex))
dev.off()

## Principle Component Analysis

# 1. Total
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
# 2. pc$x first2
plot(pc$x[,1:2])
dev.off()

# 3. scaled expressions
ex.scale <- t(scale(t(ex),scale=F))
pc <- prcomp(ex.scale)
pdf("Results/PC_scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

# 4. samples pca
gr <- c(rep('AMLP', 13), rep("B_ALL", 4), rep("T_ALL", 2), "B_ALL", "T_ALL","B_ALL", "B_ALL", "T_ALL", "B_ALL", "B_ALL", rep("T_ALL", 2), rep("B_ALL", 5), "T_ALL", "T_ALL", "B_ALL", "T_ALL", "BP", "T_ALL","AMLCL", "Granul_Norm", "BP", "T_ALL", "AMLCL", "Granul_Norm", "BP", "T_ALL", "AMLCL", "BP", rep("AMLCL", 2), "BP",rep("AMLCL", 2), "BP", rep("AMLCL", 2), "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "BP", "B_ALL", "Granul_Norm", "B_ALL", "Granul_Norm", "Mono_Norm", "Mono_Norm", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "B_ALL", rep("T_Cell_Norm", 2), "AMLCL", "B_ALL", rep("T_Cell_Norm", 2), "AMLCL", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "CD34_Norm", "T_ALL", "TP","AMLCL", "CD34_Norm", "T_ALL", "TP","AMLCL", "CD34_Norm", "T_ALL", "TP", "AMLCL", "TP", rep("BP", 9), "TP", rep("BP", 7), rep("TP", 8), rep("Granul_Norm", 7), "AMLP", "AMLP", "T_Cell_Norm", rep("AMLP", 3), rep("B_Cell_Norm", 7), "T_Cell_Norm", rep("Mono_Norm", 4), "Granul_Norm", rep("T_Cell_Norm", 7))
pcr <- data.frame(pc$r[,1:3],Group = gr)
pdf("Results/PCA_samples.pdf")
ggplot(pcr,aes(x=PC1,y=PC2,color=Group))+geom_point(size=3)+theme_bw()
dev.off()


## Correlation Heatmap
pdf("Results/pheatmap2.pdf",width = 15,height = 15)
pheatmap(cor(ex),labels_col = gr,labels_row = gr)
pheatmap(cor(ex.scale),labels_col = gr,labels_row = gr)
dev.off()

## Bonus ##

# ica
ica <- icafast((ex.scale), 2, center = TRUE, maxit = 100,tol = 1e-6)

ic <- data.frame(ICA1 = ica$Y[, 1],ICA2 = ica$Y[, 2])

pdf("Results/ICA.pdf")
ggplot(ic,aes(x=ICA1,y=ICA2))+geom_point()+theme_bw()
dev.off()

gr <- c(rep('AMLP', 13), rep("B_ALL", 4), rep("T_ALL", 2), "B_ALL", "T_ALL","B_ALL", "B_ALL", "T_ALL", "B_ALL", "B_ALL", rep("T_ALL", 2), rep("B_ALL", 5), "T_ALL", "T_ALL", "B_ALL", "T_ALL", "BP", "T_ALL","AMLCL", "Granul_Norm", "BP", "T_ALL", "AMLCL", "Granul_Norm", "BP", "T_ALL", "AMLCL", "BP", rep("AMLCL", 2), "BP",rep("AMLCL", 2), "BP", rep("AMLCL", 2), "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "AMLCL", "BP", "B_ALL", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "BP", "B_ALL", "Granul_Norm", "B_ALL", "Granul_Norm", "Mono_Norm", "Mono_Norm", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "B_ALL", rep("T_Cell_Norm", 2), "AMLCL", "B_ALL", rep("T_Cell_Norm", 2), "AMLCL", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "B_Cell_Norm", "B_ALL", "T_Cell_Norm", "AMLCL", "CD34_Norm", "T_ALL", "TP","AMLCL", "CD34_Norm", "T_ALL", "TP","AMLCL", "CD34_Norm", "T_ALL", "TP", "AMLCL", "TP", rep("BP", 9), "TP", rep("BP", 7), rep("TP", 8), rep("Granul_Norm", 7), "AMLP", "AMLP", "T_Cell_Norm", rep("AMLP", 3), rep("B_Cell_Norm", 7), "T_Cell_Norm", rep("Mono_Norm", 4), "Granul_Norm", rep("T_Cell_Norm", 7))
ica_p <- data.frame(ICA_1 = ica$Q[1,],ICA_2 = ica$Q[2,],Group=gr)
pdf("Results/ICA_samples.pdf")
ggplot(ica_p,aes(x=ICA_1,y=ICA_2,color=Group,label=gr,col=classifiaction))+geom_point()+theme_bw()
dev.off()


# tsne
tsne_cereals_num <- Rtsne(ex.scale,pca = FALSE, perplexity = 10,theta = 0.0)

tsne_num <- data.frame(TSNE1 = tsne_cereals_num$Y[, 1],TSNE2 = tsne_cereals_num$Y[, 2])

ggplot(tsne_cereals_num, aes(
  x = TSNE1, y = TSNE2,
  label = label, col = classification
)) +
  geom_point() 


# umap
umap_num <- umap(ex.scale,n_neighbors = 15,min_dist = 1, spread = 5)
umap_num <- data.frame(
  UMAP1 = umap_num[, 1],
  UMAP2 = umap_num[, 2])

ggplot(umap_num, aes(x = UMAP1, y = UMAP2,)) + geom_point() 


# autoencoder
library(ruta)
x <- as.matrix(ex.scale)
ruta_num <- autoencode(scale(x), 2, type = "robust", activation = "tanh")

yruta_num <- data.frame(RUTA1 = ruta_num[, 1],RUTA2 = ruta_num[, 2])

pdf("Results/RUTA.pdf")
ggplot(ruta_num, aes(x = RUTA1, y = RUTA2)) +geom_point()+theme_bw()
dev.off()

ruta_p <- data.frame(ruta_1 = ruta_num$Q[1,],ruta_2 = ruta_num$Q[2,],Group=gr)
pdf("Results/RUTA_samples.pdf")
ggplot(ruta_p,aes(x=ruta_1,y=ruta_2,color=Group,label=gr,col=classifiaction))+geom_point()+theme_bw()
dev.off()
 


##### Differential Expression Analysis

# assign samples to groups and set up design matrix
gs <- factor(gr)
gset$description <- gs
groups <- levels(gs)
design <- model.matrix(~description + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model


# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[2], groups[6], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("Gene.symbol","ID","adj.P.Val","logFC"))
write.table(tT, "Results/AML_Normal.txt", row.names=F, sep="\t",quote = F)

##Making Gene with increase of expression
aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(aml.up$Gene.symbol)
aml.up.genes <- unique(as.character(strsplit2(aml.up.genes, "///")))
write.table(aml.up.genes, "Results/AMLP_Normal_up.txt",quote = F, row.names = F, col.names = F)

##Making Gene with decrease of expression
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(aml.down$Gene.symbol)
aml.down.genes <- unique(as.character(strsplit2(aml.down.genes, "///")))
write.table(aml.down.genes, "Results/AMLP_Normal_down.txt", row.names = F, col.names = F, quote = F)



