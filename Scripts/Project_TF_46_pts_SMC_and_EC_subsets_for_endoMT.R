### subset the EC and SMCs

# R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"

### load packages
library(Seurat)         # V4.0.2
library(dplyr)          # v1.0.6
library(plyr)           # v1.8.6
library(SingleR)        # v1.6.1
library(ggplot2)        # v3.3.3
library(pheatmap)       # v1.0.12

set.seed(456)

seubset <- subset(seuset, idents =c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II", "ACTA2+ Smooth Muscle Cells" ))

table(Idents(seubset))
FeaturePlot(seubset, "ACTA2")
FeaturePlot(seubset, "MYH11")
FeaturePlot(seubset, "CD34")
FeaturePlot(seubset, "PECAM1")

### recluser
# find variable features
seubset <- FindVariableFeatures(seubset)

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
#View(seubset)
# PCA is performed on RNA assay

ElbowPlot(seubset)


# clustering
#spare <- seubset
seubset <- spare
seubset <- FindNeighbors(seubset, dims= 1:15)

seubset <- FindClusters(seubset, resolution = 1.2)
UMAPPlot(seubset)
VlnPlot(seubset, "ACTA2")
VlnPlot(seubset, "PECAM1")
VlnPlot(seubset, "TAGLN")

#dims14res1.1 <- seubset

DEG <- FindAllMarkers(seubset, only.pos = T)
DEG <- subset(DEG, p_val_adj < 0.05)
View(DEG)


writeClipboard(DEG[DEG$cluster == 3, "gene"])
writeClipboard(DEG[DEG$cluster == 4, "gene"])
writeClipboard(DEG[DEG$cluster == 6, "gene"])
writeClipboard(DEG[DEG$cluster == 7, "gene"])



# 3,4,6,7 = SMC
VlnPlot(seubset, "MGP")
VlnPlot(seubset, "PDGFA")
VlnPlot(seubset, "ACTA2")
VlnPlot(seubset, "MYH11")
UMAPPlot(seubset, label = T)
UMAPPlot(seubset, group.by = "new.ident")

seubset <- RenameIdents(seubset,
                        "0" = "EC1",
                        "1" = "EC2",
                        "2" = "EC3",
                        "3" = "SMC1",
                        "4" = "SMC2",
                        "5" = "EC4",
                        "6" = "SMC3",
                        "7" = "SMC4" )

UMAPPlot(seubset, label = T)
Idents(seubset) <- factor(x = Idents(seubset), levels = sort(levels(seubset)))
UMAPPlot(seubset, label = T)

# library(ReactomePA)
# library(org.Hs.eg.db)
# library(AnnotationDbi)
# library(biomaRt)
# set mart
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# reactome.resulst <- lapply(0:7, function(ident){
#   # get DEGs
#   current.DEG <- DEG <- subset(DEG, cluster == ident)
#   current.DEG <- current.DEG$gene
#   
#   # translate to entrez
#   path_genes <- getBM(                                # set source
#                       attributes = c("hgnc_symbol","entrezgene_id"),         # set symbol type
#                       values = current.DEG,                                    # which genes to extract
#                       mart = mart)
#   
#   pa <- enrichPathway(current.DEG)
#   
# })
# 
# symbols <- mapIds(org.Hs.eg.db, keys = current.DEG,
#                   column = c('ENTREZID'), keytype = 'SYMBOL')
# symbols <- symbols[!is.na(symbols)]
# symbols <- symbols[match(current.DEG, names(symbols))]
# rownames(RNA_pimcs_diff_stab_time) <- symbols
# row.names(RNA_pimcs_diff_stab_time)[is.na(rownames(RNA_pimcs_diff_stab_time))]<-ens[is.na(rownames(RNA_pimcs_diff_stab_time))]

saveRDS(seubset, file = "20210812.46.patients.Koen.cleaned.endoMT.subset.RDS")