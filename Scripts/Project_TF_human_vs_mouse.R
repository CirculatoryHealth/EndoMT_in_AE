### 2022-06-13

# Project TF Experiment Hexagonite - 2

# SingleR for human vs mice
# human vs in vitro
# in vitro vs mice

# what are the populations from humans that match mice or vice versa

#-----------------------------------------------------------------------------------------
### Settings

# set directory
setwd("")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"

### load packages
library(Seurat)         # V4.0.3 # specifically updated for singleR
library(plyr)           # v1.8.6
library(dplyr)          # v1.0.6
library(SingleR, quietly = T)        # v1.6.1
library(ggplot2)        # v3.3.3
library(pheatmap)       # v1.6.1
library(SingleCellExperiment) # v1.14.1
library(celldex)        # v1.2.0
library(enrichR)        # v3.0
library(biomaRt)        # v2.48.0

# set seed
set.seed(456)

# load data
seuset <- readRDS(file = "20210812.46.patients.Koen.cleaned.endoMT.subset.RDS")
Idents(seuset) <- factor(x = Idents(seuset), levels = sort(levels(seuset)))
seuset[["current.ident"]] <- Idents(seuset)

mice <- readRDS("20210813_endomt.Mice.subsetted.Owens.LS.RDS")
mice.cols <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702","#BB3E03")

#-----------------------------------------------------------------------------------------
UMAPPlot(mice, label = T)
UMAPPlot(seuset, label = T)

#-----------------------------------------------------------------------------------------

### run normal Single R

# the count matrix is stored in pbmc[["RNA"]]@counts . ... with the raw (non-normalized data)
# can also use seuratobj[["RNA"]]@data
#mice.sce <- as.SingleCellExperiment(mice)
mice.sce <- DietSeurat(mice)
mice.sce <- as.SingleCellExperiment(mice.sce)

# reference dataset
mice.ref <- MouseRNAseqData()

single.mice <- SingleR(test = mice.sce, 
                       ref  = mice.ref,
                       assay.type.test= 1,
                       labels =mice.ref$label.main)
plotScoreHeatmap(single.mice)

# assuming that the order of the cells is the same
mice[["singleR"]] <- single.mice@listData[["labels"]]
table(mice@meta.data$seurat_clusters, mice@meta.data$singleR)


#-----------------------------------------------------------------------------------------
# Human data

# transform into singleCellExperiment
seuset.sce <- DietSeurat(seuset) # otherwise it complains
seuset.sce <- as.SingleCellExperiment(seuset.sce)
#View(seuset.sce)

# load the reference
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

single.human <- SingleR(test = seuset.sce, 
                        ref  = hpca.se,
                        assay.type.test= 1,
                        labels = hpca.se$label.main
)
plotScoreHeatmap(single.human)

# assuming that the order of the cells is the same
seuset[["singleR"]] <- single.human@listData[["labels"]]
seuset[["active.ident"]] <- seuset@active.ident
table(seuset@meta.data$active.ident, seuset@meta.data$singleR)

#-----------------------------------------------------------------------------------------
### human vs mice

single.human.vs.mice <- SingleR(test = seuset.sce, 
                                ref  = mice.sce,
                                assay.type.test= 1,
                                labels = mice.sce@colData@listData[["ident"]]
)
plotScoreHeatmap(single.human.vs.mice)


# assuming that the order of the cells is the same
seuset[["singleR.mice"]] <- single.human.vs.mice@listData[["labels"]]
table(seuset@meta.data$active.ident, seuset@meta.data$singleR.mice)


single.mice.vs.human <- SingleR(test = mice.sce,
                                ref = seuset.sce,
                                assay.type.test = 1,
                                labels = seuset.sce@colData@listData[["ident"]])
plotScoreHeatmap(single.mice.vs.human)


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

single.human2 <- SingleR(test = seuset@assays[["SCT"]]@data, 
                         ref  = hpca.se,
                         assay.type.test= 1,
                         labels = hpca.se$label.main
)
plotScoreHeatmap(single.human2)

# assuming that the order of the cells is the same
seuset[["singleR2"]] <- single.human2@listData[["labels"]]
table(seuset@meta.data$active.ident, seuset@meta.data$singleR2)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
### plot human DEG in mice

# EC1: FN1+ EC
# EC2: activated EC
# EC3: SULF1+ EC
# EC4: ACKR1+ EC
# 
# SMC1: fibroblast
# SMC2: migratory SMC
# SMC3: TBX2+
# SMC4: contractile SMC


DEG <- FindAllMarkers(seuset, only.pos = T)
DEG <- subset(DEG, p_val_adj < 0.05)

# translate genes
convertMouseGeneList <- function(x){
  require("biomaRt")
  # human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") # server error, fixed like this
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  # genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, 
  #                  attributesL = c("hgnc_symbol"), martL = human, uniqueRows = T)
  # 
  # return(genesV2)
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, 
                   attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
  
  return(genesV2)
}

## translate gene list
# make reference list so you only have to access the biomaRt server once
mouse.genes <- convertMouseGeneList(DEG$gene)

# merge
DEG <- merge(DEG, mouse.genes, by.x = c("gene"), by.y = c("HGNC.symbol"))


# calculate module scores
human.populations <- c("EC1","EC2", "EC3", "EC4", "SMC1", "SMC2", "SMC3", "SMC4")
module.scores <- lapply(human.populations, function(p){
  
# subset genes per population
  set <- subset(DEG, cluster == p, select = MGI.symbol)
  set <- set$MGI.symbol
  n.total.genes <- length(set)
  
  # if set length in object > 5 proceed
  set <- set[c(set %in% rownames(mice))]
  n.genes.set <- length(set)
  n.genes <- paste(n.genes.set,n.total.genes, sep = "/")
  
  if(n.genes.set >= 5){
    seubset <- AddModuleScore(object = mice,
                              features = list(set),
                              name = "test",
                              nbin = 10,
                              ctrl = 100)
    
    p1 <- VlnPlot(seubset, "test1", cols = mice.cols)  + 
      labs(title = p, subtitle = n.genes) +      # set title and subtitle
      theme(plot.subtitle = element_text(hjust = 0.5)) +    # adjust placements of subtitle
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) +   # adjust x axis
      theme(legend.position = "none")# +
      #ylim(-0.5, 1)  #adjust y axis
    
    p1
    
    return(p1)
    
  } else {
    return(NULL)
  }
  
  
})
names(module.scores) <- human.populations
#View(module.scores)

patchwork::wrap_plots(module.scores, 3,3)

pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_human_DEG_projected_in_Mice_populations.pdf", sep = ""), height = 11.45 , width = 11.45, useDingbats = FALSE)
patchwork::wrap_plots(module.scores, 3,3)
dev.off()

#-----------------------------------------------------------------------------------------
### pick through mouse DEG

DEG2 <- FindAllMarkers(mice, only.pos = T)
DEG2 <- subset(DEG2, p_val_adj < 0.05)
dim(DEG2)

dbs_mod<-c("GO_Biological_Process_2021")

# test <- enrichr(DEG2$gene[1:50],databases = dbs_mod)
# plotEnrich(test[[1]])

mouse.idents <- c(0,1,2,3,4,5,6)
mouse_EnrichR <- lapply(mouse.idents, function(i){
  
  genelist <- DEG2[DEG2$cluster == i, "gene"]
  
  dff <- enrichr(genelist,databases = dbs_mod)
  p <- plotEnrich(dff[[1]]) +
    ggtitle(label = paste("population", i, sep = " "))
  
  return(p)
})
names(mouse_EnrichR) <- mouse.idents
View(mouse_EnrichR)


# online EnrichR
writeClipboard(DEG2[DEG2$cluster == 0, "gene"])
writeClipboard(DEG2[DEG2$cluster == 1, "gene"])
writeClipboard(DEG2[DEG2$cluster == 2, "gene"])
writeClipboard(DEG2[DEG2$cluster == 3, "gene"])
writeClipboard(DEG2[DEG2$cluster == 4, "gene"])
writeClipboard(DEG2[DEG2$cluster == 5, "gene"])
writeClipboard(DEG2[DEG2$cluster == 6, "gene"])

## other markers (as in opal)
## neovascular markers
VlnPlot(mice, "Ackr1", cols = mice.cols)
VlnPlot(mice, "Sulf1", cols = mice.cols)
VlnPlot(mice, "Pecam1", cols = mice.cols)

# markers from human
VlnPlot(mice, "Fn1", cols = mice.cols)
VlnPlot(mice, "Tbx2", cols = mice.cols)


# pan et al
SEM <- c("Ly6a", "Vcam1", "Ly6c1") #Ly6a, Vcam1, and Ly6c1
VlnPlot(mice, SEM, cols = mice.cols)

# contractile 
con <- c("Des", "Smtn", "Myh11", "Acta2")
VlnPlot(mice, con, cols = mice.cols)

# hu et al
FB <- c("Apod", "Cfd", "Fbln1", "Dcn", "Sfrp2" )
VlnPlot(mice, FB, cols = mice.cols)
MyoFB <- c("S100b", "Mpz", "Gpm6b","Plp1", "Nrxn1")
VlnPlot(mice, MyoFB, cols = mice.cols)

# wirka et al
Fib <- c("Pi16", "Clec3b", "Gpx3", "Serping1", "Cygb", "Dpep1", "Smoc") # Pi16 Clec3b Gpx3 Serping1 Cygb Dpep1 Smoc
DotPlot(mice, features = Fib) + theme_test()

#Fib2 <- c("GSN", "COMP", "CHAD", "WIF1", "THBS1", "FMOD", "ANGPT17", "PRG4") # Gsn Comp Chad Wif1 1500015O10Rik Thbs1 Fmod Angptl7 Prg4 
Fib2 <- c("Gsn", "Comp", " Chad", "Wif1", "1500015O10Rik", "Thbs1", "Fmod", "Angptl7", "Prg4")
DotPlot(mice, features =  Fib2) + theme_test()

#MSMC <- c("PRG4", "SPP1", "IBSP", "FN1", "COL2A1", "LCN1", "LUM", "TIMP1") # Prg4 Spp1 Ibsp Fn1 Col2a1 Lcn2 Lum Timp1 Tnfrsf11b
MSMC <- c("Prg4", "Spp1", "Ibsp", "Fn1", "Col2a1", "Lcn2", "Lum"," Timp1", "Tnfrsf11b")
DotPlot(mice, features =  MSMC) + theme_test()

# Alex et al
pericytes <- c("Pdgfrb", "Adamts1", "Vtn", "Cspg4")
DotPlot(mice, features =  pericytes) + theme_test()
VlnPlot(mice, "Adamts1", cols = mice.cols)
VlnPlot(mice, "Cspg4", cols = mice.cols)
# CD146
VlnPlot(mice, "CD146", cols = mice.cols)


# more other
VlnPlot(mice, "Acta2", cols = mice.cols)
VlnPlot(mice, "Myh11", cols = mice.cols)

# mouse 0 = ECM producing cells
# mouse 1 = Fibroblasts
# mouse 2 = (neo)Vascular ECs 1
# mouse 3 = (neo)Vascular ECs 2
# mouse 4 = Snai1+ cells
# mouse 5 = activated ECs
# mouse 6 = Contractile SMCs