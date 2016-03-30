#! /usr/bin/Rscript

load("../result/GSE13355_Psoriasis.RData")

library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# data preprocess
GSE13355_GEP <- GSE13355.Explist.filtered$GSE13355$x

datExpr <- t(GSE13355_GEP)
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# metwork construct
net = blockwiseModules(datExpr, power = 7,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "GSE13355",
                       verbose = 3)


table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)

gene_cluster_info <- moduleColors
names(gene_cluster_info) <- rownames(GSE13355_GEP)


# Pathway analysis based on clusters
library(cogena)
source("./lib/enrichTable.R")
source("./lib/heatmapTable.R")
source("./lib/raw_clEnrich.R")

annoGMT <- "c2.cp.kegg.v5.0.symbols.gmt.xz"
annofile <- system.file("extdata", annoGMT, package="cogena")
clEnrich_res <- raw_clEnrich(gene_cluster_info, annofile=annofile)


heatmapTable(clEnrich_res, maintitle="GSE13355 WGCNA Pathway Analysis")



# GO analysis
library(org.Hs.eg.db)
eg_Symbol <- toTable(org.Hs.egSYMBOL)

allLLIDs <- match(rownames(GSE13355_GEP), eg_Symbol$symbol)

GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, nBestP = 10)
tab = GOenr$bestPTerms[[4]]$enrichment

save.image("../result/WGCNA/GSE13355_WGCNA.RData")

