### 2021-8-13

# Project TF Experiment Amethyst-2

# Check out data, subset data from mice
# old data was corrupted, new Seurat version available

#-----------------------------------------------------------------------------------------
### load data 
# set directory
setwd("")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"

### load packages
library(dplyr)          # v1.0.6
library(Seurat)         # v4.0.2
library(ggplot2)        # v3.3.3
#library(ggplotify)      # v0.0.5
#library(cowplot)        # v1.1.0
#library(RColorBrewer)   # v1.1-2
#library(gplots)         # v3.1.0
#library(pheatmap)       # v1.0.12
#library(reshape2)       # v1.4.4

# set seed
set.seed(456)

## Get single cell data
# mice
mice <- readRDS(file = "" )

mice <- UpdateSeuratObject(mice)

#-----------------------------------------------------------------------------------------
# View(head(mice@meta.data))
# UMAPPlot(mice)
# 
# # IL1BWTeYFPpos
# experiment <- "IL1BWTeYFPpos"
# 
# coi <- rownames(mice@meta.data[mice@meta.data$orig.ident == "IL1BWTeYFPpos",])
# 
# endomt <- subset(mice, cells = coi)
# 
# FeaturePlot(endomt, "Acta2")
# FeaturePlot(endomt, "Myh11")
# FeaturePlot(endomt, "Vim")
# FeaturePlot(endomt, "Tagln")
# 
# # get rid of "lost"cells
# # cluster 9, 10, 11, 12
# endomt <- subset(endomt, idents = c(1,2,3,4,5,6,8))
# 
# UMAPPlot(endomt)
# FeaturePlot(endomt, "Acta2")
# FeaturePlot(endomt, "Cdh5")

# #------------------------------------------------------------------------------------------
# ### stuff
# #cdh5.pos <- subset(endomt, expression = "Cdh5" > 1)
# 
# max(GetAssayData(endomt)["Cdh5",])
# min(GetAssayData(endomt)["Cdh5",])
# 
# cut1 <- WhichCells(endomt, expression = Cdh5 > 1)
# cut2 <- WhichCells(endomt, expression = Cdh5 > 2)
# cut3 <- WhichCells(endomt, expression = Cdh5 > 3)
# 
# UMAPPlot(endomt)
# FeaturePlot(endomt, features = "Cdh5", cells = cut1)
# FeaturePlot(endomt, features = "Cdh5", cells = cut2)
# FeaturePlot(endomt, features = "Cdh5", cells = cut3)
# 
# FeaturePlot(endomt, features = "Cd34")
# FeaturePlot(endomt, features = "Pecam1")
# 
# FeaturePlot(endomt, features = c("Cd34", "Acta2"), blend = T)
# FeaturePlot(endomt, features = "Zeb1")
# FeaturePlot(endomt, features = "Zeb2")
# FeaturePlot(endomt, features = "Snai1")
# FeaturePlot(endomt, features = "Snai2")
# FeaturePlot(endomt, features = "Runx2")
# FeaturePlot(endomt, features = "Wwtr1")
# 
# FeaturePlot(endomt, features = c("Myh11", "Acta2"), blend = T)
# max(GetAssayData(endomt)["Acta2",])
# min(GetAssayData(endomt)["Acta2",])
# 
# plot(density(GetAssayData(endomt)["Acta2",]))
# pdf(paste(result.folder, "/", Sys.Date(), "_mice_markers.pdf", sep = ""), useDingbats = FALSE)
# VlnPlot(endomt, "Acta2")
# VlnPlot(endomt, "Myh11")
# VlnPlot(endomt, "Pecam1")
# VlnPlot(endomt, "Cd34")
# dev.off()
# 
# DEG <- FindAllMarkers(endomt)
# DEG <- DEG[DEG$p_val_adj <= 0.05,]
#------------------------------------------------------------------------------------------
View(head(mice@meta.data))
UMAPPlot(mice)

unique(mice@meta.data$origin)
# mice[["plot.origin"]] <- as.numeric(factor(mice@meta.data$origin))
# mice[["plot.origin"]] <- (factor(mice@meta.data$origin))
# FeaturePlot(mice, features = "plot.origin")

unique(mice@meta.data$orig.ident)
table(mice@meta.data$orig.ident, mice@meta.data$origin)

# endothelial_eYFPpos
experiment <- "Endothelia_IL1B_WT_eYFP_Positive"

coi <- rownames(mice@meta.data[mice@meta.data$origin == "Endothelia_IL1B_WT_eYFP_Positive",])

endomt <- subset(mice, cells = coi)
UMAPPlot(endomt)
FeaturePlot(endomt, "Acta2")
FeaturePlot(endomt, "Cdh5")

# # get rid of "lost"cells # why are they lost?
# # cluster 9, 10, 11, 12
#endomt <- subset(endomt, idents = c(1,2,3,4,5,6,8))
UMAPPlot(endomt)
dim(endomt)

VlnPlot(endomt, "Acta2")
VlnPlot(endomt, "Myh11")
VlnPlot(endomt, "Pecam1")
VlnPlot(endomt, "Cd34")

# ## group together idents
# # 1 = SMC
# # 8 = EC
# # 2 and 3 = transitional cluster 2
# # 4,5 and 6 = transitional cluster 1
# endomt[["seurat_clusters"]] <- Idents(endomt)
# 
# endomt[["new.groups"]] <- NA
# 
# endomt@meta.data[endomt@meta.data$seurat_clusters == 1, "new.groups"] <- "Smooth Muscle Cells"
# endomt@meta.data[endomt@meta.data$seurat_clusters == 8, "new.groups"] <- "Endothelial Cells"
# endomt@meta.data[endomt@meta.data$seurat_clusters == 2 | endomt@meta.data$seurat_clusters == 3 , "new.groups"] <- "Transitional Cluster 2"
# endomt@meta.data[endomt@meta.data$seurat_clusters == 4 | endomt@meta.data$seurat_clusters == 5 | endomt@meta.data$seurat_clusters == 6 , "new.groups"] <- "Transitional Cluster 1"
# 
# table(endomt@meta.data$new.groups
# )


### recluser
# find variable features
seubset <- endomt

seubset <- FindVariableFeatures(seubset, nfeatures = 500)

top10 <- head(VariableFeatures(seubset), 10)
plot1 <- VariableFeaturePlot(seubset)
plot2 <- LabelPoints(plot = plot1, points = top10, repel =T)
#plot1 + plot2
plot2

# scale all features
seubset <- ScaleData(seubset, features = rownames(seubset))

# run PCA with variable features
seubset <- RunPCA(seubset, features = VariableFeatures(object = seubset), verbose = F)
print(seubset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seubset, dims = 1:2, reduction = "pca")

DimPlot(seubset, reduction = "pca")
View(seubset)

# elbow plot
ElbowPlot(seubset)


# clustering
#spare <- seubset
seubset <- spare
seubset <- FindNeighbors(seubset, dims= 1:16)

# find clusters
seubset <- FindClusters(seubset, resolution = 0.7)
seubset <- FindClusters(seubset, resolution = 1.1)

# compute umap
seubset <- RunUMAP(seubset, dims = 1:15)

UMAPPlot(seubset)
VlnPlot(seubset, "Acta2")
VlnPlot(seubset, "Pecam1")
VlnPlot(seubset, "Tagln2")
VlnPlot(seubset, "Cdh5")


FeaturePlot(seubset, "Acta2")
FeaturePlot(seubset, "Pecam1")
FeaturePlot(seubset, "Cdh5")

table(seubset@meta.data$OriginalClusters, seubset@meta.data$seurat_clusters)

# decide on meta data to keep
keep <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "origin" , "percent.mito","percent.hemo","tissue","orig.cell", 
          "OriginalClusters", "integrated_snn_res.0.7", "seurat_clusters")

seubset@meta.data[,!c(names(seubset@meta.data) %in% keep)] <- NULL
str(seubset@meta.data)

# scale all features in RNA assay
seubset <- ScaleData(seubset, features = rownames(seubset), assay = "RNA")
VlnPlot(seubset, "Acta2", assay = "RNA", slot = "scale.data")
VlnPlot(seubset, "Pecam1", assay = "RNA")
VlnPlot(seubset, "Cdh5", assay = "RNA")

# To keep this simple: You should use the integrated assay when trying to 'align' cell states that are shared across datasets 
# (i.e. for clustering, visualization, learning pseudotime, etc.)
# You should use the RNA assay when exploring the genes that change either across clusters, trajectories, or conditions.
DEG <- FindAllMarkers(seubset, assay = "RNA")
DEG <- subset(DEG, p_val_adj < 0.05)
dim(DEG)

#saveRDS(seubset, file = "20210813_endomt.Mice.subsetted.Owens.LS.RDS")
