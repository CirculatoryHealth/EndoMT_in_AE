### 2022-04-07

# Project TF Experiment Ruby - 2

# mice pseudotime for paper (all in one script)
# with multiple lineages option

#-----------------------------------------------------------------------------------------
### load data 
# set directory
setwd("")

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 4.0.2 (2020-06-22)

### load packages
library(dplyr)          # v1.0.6
library(Seurat)         # v4.0.2 #v.4.0.3
library(ggplot2)        # v3.3.3
library(slingshot, quietly = T) # V2.0.0
library(enrichR)        # v3.0
library(ggpubr)         # v0.4.0
library(patchwork)      # v1.1.1

# set seed
set.seed(456)

# load data
seuset <- readRDS("20210813_endomt.Mice.subsetted.Owens.LS.RDS")

# k means
kmeans <- read.delim2("clustering_9.txt")
names(kmeans) <- c("symbol", "all", "control", "TGFb", "TNFa", "mix")
#conditions <- c("control", "TGFb", "TNFa", "mix", "all")
conditions <- c("control", "TGFb", "TNFa", "mix")
n.kmeans <- 9


# enrichR database
dbs_mod<-c("GO_Biological_Process_2018","GO_Biological_Process_2015")
#dbs_mod<-c("GO_Biological_Process_2021","GO_Biological_Process_2018")

#-----------------------------------------------------------------------------------------
### try slingshot
# make meta data
meta.data <- data.frame(cells = rownames(seuset@meta.data),
                        seurat_clusters = as.character(seuset@meta.data$seurat_clusters),
                        OriginalClusters = as.character(seuset@meta.data$OriginalClusters)
)

mice.cols <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702","#BB3E03")
meta.data[meta.data$seurat_clusters == 0, "colours"] <- "#005F73" 
meta.data[meta.data$seurat_clusters == 1, "colours"] <- "#0A9396"
meta.data[meta.data$seurat_clusters == 2, "colours"] <- "#94D2BD"
meta.data[meta.data$seurat_clusters == 3, "colours"] <- "#E9D8A6"
meta.data[meta.data$seurat_clusters == 4, "colours"] <- "#EE9B00"
meta.data[meta.data$seurat_clusters == 5, "colours"] <- "#CA6702"
meta.data[meta.data$seurat_clusters == 6, "colours"] <- "#BB3E03"

# mouse 0 = ECM producing cells
# mouse 1 = Fibroblasts
# mouse 2 = (neo)Vascular Endothelial Cells 1
# mouse 3 = (neo)Vascular Endothelial Cells 2
# mouse 4 = Snai1+ cells
# mouse 5 = 
# mouse 6 = Contractile SMCs

# https://bustools.github.io/BUS_notebooks_R/slingshot.html
# we can embed the information from seurat object directly into singecellexperiment from slingshot
# possibility to identify the start cluster

# sds <- slingshot(Embeddings(seuset, "umap"), 
#                  clusterLabels = meta.data$seurat_clusters, 
#                  end.clus = '1', 
#                  stretch = 1)

sds <- slingshot(Embeddings(seuset, "umap"), 
                 clusterLabels = meta.data$seurat_clusters, 
                 end.clus = '6', 
                 stretch = 1)

# plot
plot(sds@elementMetadata@listData[["reducedDim"]], 
     col= meta.data$colours, 
     pch = 16, 
     cex = 0.75)
lines(SlingshotDataSet(sds), 
      lwd=2, 
      type = 'lineages', 
      col = 'black')


lin2 <- getLineages(sds@elementMetadata@listData[["reducedDim"]],
                    clusterLabels = meta.data$seurat_clusters,
                    #start.clus = "4",
                    end.clus = "6"
)

# lin2 <- getLineages(sds@elementMetadata@listData[["reducedDim"]],
#                     clusterLabels = meta.data$seurat_clusters
#                     #start.clus = "4",
#                     #end.clus = "6"
# )

# Constructing smooth curves and ordering cells
crv1 <- getCurves(lin2)
crv1

#pdf(paste(result.folder, "/", Sys.Date(), "_projectTF_mice_UMAP_lineages_pop_4.pdf", sep = ""), useDingbats = FALSE)

plot(sds@elementMetadata@listData[["reducedDim"]], 
     col= meta.data$colours, 
     pch = 16, 
     cex = 0.75)

plot(sds@elementMetadata@listData[["reducedDim"]], 
     col= meta.data$colours, 
     pch = 16, 
     cex = 0.75)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')


#-----------------------------------------------------------------------------------------
### translate to mouse symbols

# translate human HGNC to mice 
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
mouse.genes <- convertMouseGeneList(kmeans$symbol)

# merge
#lookup <- merge(kmeans, mouse.genes, by.x = c("symbol"), by.y = c("HGNC.symbol"))
kmeans <- merge(kmeans, mouse.genes, by.x = c("symbol"), by.y = c("HGNC.symbol"))


#dim(lookup)
dim(mouse.genes)
dim(kmeans)
head(kmeans)
# you loose about 400 genes, quite a chunk 

#write.csv2(lookup, "kmeans_9_mouse_human_translated.txt")


#-----------------------------------------------------------------------------------------
### Pseudotime and module scores

# for each K-means, calculate the module scores
# n genes not really representative
module.scores <- lapply(conditions, function(name){
    ### per condition, plot module score per K means module
  current.set <- kmeans[,c(name, "MGI.symbol")]
  
  # per k
  res <- lapply(1:n.kmeans, function(k){
    # select all genes that are in this group
    set <- current.set[current.set[,1] == k, "MGI.symbol"]
    n.total.genes <- length(set)
    
    # see if they are in the data
    set <- set[c(set %in% rownames(seuset))]
    n.genes.set <- length(set)
    n.genes <- paste(n.genes.set,n.total.genes, sep = "/")
    
    seuset <- AddModuleScore(object = seuset,
                             features = list(set),
                             name = "test",
                             nbin = 10,
                             ctrl = 100)
    
    scores <- data.frame(geneset = seuset@meta.data$test1,
                         row.names = rownames(seuset@meta.data))
    names(scores) <- name
    
    
    return.list <- list(scores, n.genes)
    names(return.list) <- c("scores", "n genes")
    return(return.list)
    
    
  })
  
  names(res) <- c(1:n.kmeans)
  return(res)
  
})
names(module.scores) <- conditions
#View(module.scores)
#-----------------------------------------------------------------------------------------
### the lineages

crv1@metadata[["lineages"]]
# interesting lineages:
# lin2 and lin3

### function that per lineage
# extracts the cells from the given lineage 
# plots the module scores per k means against pseudotime
# calculates EnrichR and plots

# weights
weight <- data.frame(crv1@assays@data@listData[["weights"]])

# population labels
populations <- c("0", "1", "2","3", "4", "5", "6", "7")

# plot limits
y.max = 0.8
y.min = -0.4

plotlist <- function(lin){
  # subset the cells based on weight
  # if a cell from that lineage has a weight > 0.5, include it.
  weight1 <- weight[,lin, drop = F]
  weight1 <- weight1[weight1[,1] > 0.5, , drop = F]
  # this will introduce some NA values later, because there are also weights for cells not fully belonging to the lineage
  
  #dim(weight1)
  # cells of interest
  coi <- rownames(weight1)
  
  # subset seuset 
  seubset <- subset(seuset, cells = coi)
  
  
  ## plotting rules
  # determine levels for plotting
  chosen.levels <- c(crv1@metadata[["lineages"]][[lin]])
  # determine colours for plotting
  chosen.colours <- mice.cols[c(as.numeric(chosen.levels) + 1)]
  # determine population labels
  population.labels <- populations[c(as.numeric(chosen.levels) + 1)]
  
  # take only the cells that belong to the lineage clusters
  coi <- rownames(seubset@meta.data[seubset@meta.data$seurat_clusters %in% chosen.levels,])
  seubset <- subset(seuset, cells = coi)
  
  
  # scores 
  scores <- lapply(conditions, function(c){
    # get all corresponding k means
    input <- lapply(1:n.kmeans, function(k){
      input.ggplot <- module.scores[[c]][[k]][["scores"]]
      names(input.ggplot)[1] <- "scores"
      
      # subset input based on selected cells
      input.ggplot <- input.ggplot[rownames(seubset@meta.data),,drop = F]
      
      # plot per condition 
      # add other information
      # population
      input.ggplot$current.ident <- seubset@meta.data$current.ident
      input.ggplot$population <- seubset@meta.data$seurat_clusters
      input.ggplot$population.fctr <- factor(input.ggplot$population, levels = chosen.levels)
      # pseudotime
      pseudotime <- data.frame(crv1@assays@data@listData[["pseudotime"]])
      input.ggplot$pseudotime <- pseudotime[rownames(input.ggplot),lin]
      # cells
      input.ggplot$cells <- rownames(input.ggplot)
      # kmeans
      input.ggplot$k.means <- k
      # condition
      input.ggplot$conditon <- c
      
      return(input.ggplot)
    })
    input.ggplot <- do.call(rbind.data.frame, input)
    names(input.ggplot)[1] <- "scores"
    # remove the NA values here
    input.ggplot <- na.omit(input.ggplot)
    
    
    #plot
    p1 <- ggplot(input.ggplot, aes(x = pseudotime, y = scores)) +
      geom_point(aes(colour = population.fctr), size = 0.8) + # put colour aes here, so the smooth function works on all data
      ylim(y.min, y.max) +
      theme_bw() +
      theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
      ylab(label = element_text("score")) +
      #xlab(label = element_blank()) +
      scale_colour_manual(values = chosen.colours, labels = population.labels) +
      #geom_smooth(method = "lm", se = FALSE, formula= y~x, colour = "black") +
      #stat_cor(method = "pearson", label.y = y.max - 0.1) +
      labs(colour = "Population") 
    
    p1.1 <- p1 + facet_grid(rows = vars(k.means)) + ggtitle(label = c)
    
    p2 <- p1 + 
      stat_cor(method = "pearson", label.y = y.max - 0.1) +
      geom_smooth(method = "lm", se = FALSE, formula= y~x, colour = "black")
    
    p2.2 <- p2 + facet_wrap(facets = vars(k.means)) + ggtitle(label = c)
    
    p2.2
    
    ### enrichR
    ## per condition, plot module score per K means module
    current.set <- kmeans[,c(c, "MGI.symbol")]
    
    # per k
    res <- lapply(1:n.kmeans, function(k){
      # select all genes that are in this group
      set <- current.set[current.set[,1] == k, "MGI.symbol"]
      n.total.genes <- length(set)
      
      # see if they are in the data
      set <- set[c(set %in% rownames(seubset))]
      n.genes.set <- length(set)
      n.genes <- paste(n.genes.set,n.total.genes, sep = "/")
      
      # enrichR
      require(enrichR)
      dff <- enrichr(set,databases = dbs_mod)
      dff <- dff[["GO_Biological_Process_2018"]]
      #dff <- dff[[1]]
      dff$condition <- c
      dff$k.means <- k
      
      return.list <- list(n.genes, dff)
      names(return.list) <- c("n genes", "enrichR")
      return(return.list)
    })
    names(res) <- 1:n.kmeans
    
    ## get only first 5 pathways
    input2 <- lapply(1:n.kmeans, function(k){
      data <- res[[k]][["enrichR"]]
      # get only first 5 pathways
      data <- data[1:5,]
    })
    input.ggplot2 <- do.call(rbind.data.frame, input2)
    
    # calculate proportion of overlap (makes it easier to fix the axis as well)
    split <- lapply(input.ggplot2$Overlap, function(x){
      to.split <- strsplit(x, "[/]")
      to.split <- as.numeric(unlist(to.split))
      to.split <- to.split[1]/to.split[2]
    })
    input.ggplot2$ratio <- unlist(split)
    
    # shorten term for plotting
    # remove GO term
    short_terms <- sub("\\(.*)", "", input.ggplot2$Term)
    # pick first 40 characters
    short_terms <- strtrim(short_terms, 40)
    input.ggplot2$short_term <- short_terms
    
    # plot
    p3 <- ggplot(input.ggplot2, aes(x = ratio, y = tidytext::reorder_within(short_term, Adjusted.P.value, k.means, sep = "..."), fill = Adjusted.P.value)) + # https://www.geeksforgeeks.org/how-to-reorder-barplots-with-facetting-with-ggplot2-in-r/
      # use reoder_within to sort the bar graph on Term according to p value for the facet k means
      geom_col() + 
      #coord_flip() + 
      theme_bw() +
      scale_fill_continuous(low = "#FFC300", high = "#F06543") +
      #theme(axis.text.y=element_blank()) +
      xlim(0, 0.55) +
      scale_y_discrete(limits=rev) + # reverse the order of y axis plotting
      ylab(label = element_blank()) +
      labs(fill = "Adj. p value") +
      guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
    
    #p3
    
    p3.1 <- p3 + facet_grid(rows = vars(k.means), scales = "free_y") + # free y to plot only terms belonging to that K
      ggtitle(label = c)
    
    #p3.1 
    return.list <- list(p2.2, p3.1, res, coi, input.ggplot)
    names(return.list) <- c("pseudotime", "enrichR", "enrichR info", "cells of lineage", "input for pseudotime")
    return(return.list)
  })
  names(scores) <- conditions
  return(scores)
  #print("done")
}

### lin 1
lin1 <- plotlist("Lineage1")
#View(lin1)

# wide pdf
pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_mice_pseudotime_vs_", n.kmeans,"_all_conditions_lineage1_0.5.wide.pdf", sep = ""), height =8.25 , width = 11.45, useDingbats = FALSE)
lin1[["control"]][["pseudotime"]] + lin1[["TGFb"]][["pseudotime"]] + lin1[["TNFa"]][["pseudotime"]] + lin1[["mix"]][["pseudotime"]] +
  plot_layout(guides='collect') + plot_annotation(title = "Lineage 1") &
  theme(legend.position='bottom') 
dev.off()

# enrichR
pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_mice_EnrichR_", n.kmeans,"_all_conditions_lineage1.tall.pdf", sep = ""), height = 11.45 , width = 8.25, useDingbats = FALSE)
lin1[["control"]][["enrichR"]] + lin1[["TGFb"]][["enrichR"]] + lin1[["TNFa"]][["enrichR"]] + lin1[["mix"]][["enrichR"]] +
  plot_layout(guides='collect') + plot_annotation(title = "Lineage 1") &
  theme(legend.position='bottom')
dev.off()


### lin 2
lin2 <- plotlist("Lineage2")
#View(lin2)

# wide pdf
pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_mice_pseudotime_vs_", n.kmeans,"_all_conditions_lineage2_0.5.wide.pdf", sep = ""), height =8.25 , width = 11.45, useDingbats = FALSE)
lin2[["control"]][["pseudotime"]] + lin2[["TGFb"]][["pseudotime"]] + lin2[["TNFa"]][["pseudotime"]] + lin2[["mix"]][["pseudotime"]] +
  plot_layout(guides='collect') + plot_annotation(title = "Lineage 2") &
  theme(legend.position='bottom') 
dev.off()

# enrichR
pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_mice_EnrichR_", n.kmeans,"_all_conditions_lineage2.tall.pdf", sep = ""), height = 11.45 , width = 8.25, useDingbats = FALSE)
lin2[["control"]][["enrichR"]] + lin2[["TGFb"]][["enrichR"]] + lin2[["TNFa"]][["enrichR"]] + lin2[["mix"]][["enrichR"]] +
  plot_layout(guides='collect') + plot_annotation(title = "Lineage 2") &
  theme(legend.position='bottom')
dev.off()

### lin 3
lin3 <- plotlist("Lineage3")
#view(lin3)

# wide pdf
pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_mice_pseudotime_vs_", n.kmeans,"_all_conditions_lineage3_0.5.wide.pdf", sep = ""), height =8.25 , width = 11.45, useDingbats = FALSE)
lin3[["control"]][["pseudotime"]] + lin3[["TGFb"]][["pseudotime"]] + lin3[["TNFa"]][["pseudotime"]] + lin3[["mix"]][["pseudotime"]] +
  plot_layout(guides='collect') + plot_annotation(title = "Lineage 3") &
  theme(legend.position='bottom') 
dev.off()

# enrichR
pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_mice_EnrichR_", n.kmeans,"_all_conditions_lineage3.tall.pdf", sep = ""), height = 11.45 , width = 8.25, useDingbats = FALSE)
lin3[["control"]][["enrichR"]] + lin3[["TGFb"]][["enrichR"]] + lin3[["TNFa"]][["enrichR"]] + lin3[["mix"]][["enrichR"]] +
  plot_layout(guides='collect') + plot_annotation(title = "Lineage 3") &
  theme(legend.position='bottom')
dev.off()

#-----------------------------------------------------------------------------------------

# do a happy dance