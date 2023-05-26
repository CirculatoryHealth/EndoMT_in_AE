### 2022-02-14

# Project TF Experiment Quartz

# In vitro data from Arjan
# Pathway analysis from K means + expression data

#

#-----------------------------------------------------------------------------------------
### load data 
# set directory
setwd()

# create folder
result.folder <- paste(Sys.Date(), "results", sep = " ")
dir.create(result.folder, showWarnings = FALSE)

# R version 4.0.2 (2020-06-22)

# R version 4.0.2 (2020-06-22)

### load packages
library(dplyr)          # v1.0.6
library(tidyr)          # v1.1.3
#library(Seurat)        # v4.0.3
library(ggplot2)        # v3.3.3
library(enrichR)        # v3.0.0
library(patchwork)      # v1.1.1
library(matrixStats)    # v0.59.0
library(readr)          # v1.4.0
library(tibble)         # v3.1.2
library(enrichR)        # v3.0
library(pcaExplorer)    # v2.18.0

# set seed
#set.seed(456)

## data
# k means
## k means
kmeans <- read.delim2("clustering_9.txt")
names(kmeans) <- c("symbol", "all", "control", "TGFb", "TNFa", "mix")
#conditions <- c("control", "TGFb", "TNFa", "mix", "all")
conditions <- c("control", "TGFb", "TNFa", "mix")
n.kmeans <- 9

# expression data
expression_data <- read_tsv("gene_expression_df.txt")
metadata <- read_tsv("metadata.txt")
deg_genes <- read_tsv("deg_genes.txt")
n_clusters <- n.kmeans

dim(expression_data)
dim(metadata)

# enrichR database
dbs_mod<-c("GO_Biological_Process_2018","GO_Biological_Process_2015")

#-----------------------------------------------------------------------------------------
# ### run Arjan's code
# # 3.2 K-means on all samples
# # # center with 'colMedians()'
# center_colmedians <- function(x) {
#   tx <- t(x)
#   xcenter = colMedians(tx)
#   txc <- tx - rep(xcenter, rep.int(nrow(tx), ncol(tx)))
#   t(txc)
# }
# 
# # # apply it
# expression_data2 <- expression_data
# expression_data[-48] <- center_colmedians(as.matrix(expression_data[-48]) )
# 
# # expression_data <- expression_data2
# # expression_data[-48] <- scale(expression_data[-48])
# 
# kmeans_expr_data <- expression_data %>% 
#   filter(symbol %in% deg_genes$symbol) %>% 
#   distinct(symbol, .keep_all = TRUE) %>% 
#   remove_rownames() %>%
#   column_to_rownames("symbol")
# # %>% 
# # pivot_longer(cols = -symbol, names_to = "samples", values_to = "expression") %>% 
# # pivot_wider(names_from = symbol, values_from = expression) %>% 
# # column_to_rownames("samples")
# 
# set.seed(30921)
# c1 <- kmeans(kmeans_expr_data, centers = n_clusters)
# 
# # pull out cluster definitions
# all_stim <- c1$cluster
# all_stim <- all_stim %>% 
#   as.data.frame %>% 
#   setNames("all")
# 
# kclusterdf <- c1$cluster %>% 
#   as.data.frame() %>% 
#   rename(kcluster = 1) %>% 
#   rownames_to_column("symbol")
# 
# metadata <- metadata %>% 
#   mutate_at(vars(duplicate, time_point), list(as.factor))
# # metadata$time_point <- as.factor(metadata$time_point)
# 
# # 3.3 K-means plot - k-means on all samples
# kmeans_plot_data <- kmeans_expr_data %>% 
#   rownames_to_column("symbol") %>% 
#   pivot_longer(cols = -symbol, names_to = "sample_name", values_to = "expression_value") %>% 
#   left_join(kclusterdf, by = "symbol") %>% 
#   left_join(metadata 
#             # %>% 
#             #   rownames_to_column("sample_name_rows")
#             ,
#             by = c("sample_name" = "sample_name_rows") )
# 
# kmeans_plot_data %>% 
#   ggplot(aes(x = time_point, y = expression_value)) +
#   # geom_boxplot(aes(colour = duplicate) ) +
#   geom_boxplot() +
#   facet_grid(kcluster~stimulus, scales = "free_y")
# 
# 
# # 3.4 K-means on control samples
# 
# kmeans_expr_data_c <- expression_data %>% 
#   filter(symbol %in% deg_genes$symbol) %>% 
#   distinct(symbol, .keep_all = TRUE) %>% 
#   remove_rownames() %>%
#   column_to_rownames("symbol") %>% 
#   select(ends_with("C"))
# # %>% 
# # pivot_longer(cols = -symbol, names_to = "samples", values_to = "expression") %>% 
# # pivot_wider(names_from = symbol, values_from = expression) %>% 
# # column_to_rownames("samples")
# 
# set.seed(30921)
# c1 <- kmeans(kmeans_expr_data_c, centers = n_clusters)
# 
# # pull out cluster definitions
# control_stim <- c1$cluster
# control_stim <- control_stim %>% 
#   as.data.frame %>% 
#   setNames("control")
# 
# kmeans_plot_data %>% 
#   ggplot(aes(x = time_point, y = expression_value) ) +
#   # geom_boxplot(aes(colour = duplicate) ) +
#   geom_boxplot() +
#   facet_grid(kcluster~stimulus, scales = "free_y")
# 
# # I get it somewhat

#-----------------------------------------------------------------------------------------
### plot the k means
# steal some from Arjan;s code, 
# # center with 'colMedians()'
center_colmedians <- function(x) {
  tx <- t(x)
  xcenter = colMedians(tx)
  txc <- tx - rep(xcenter, rep.int(nrow(tx), ncol(tx)))
  t(txc)
}

# # apply it
expression_data2 <- expression_data
expression_data[-48] <- center_colmedians(as.matrix(expression_data[-48]) )

# expression_data <- expression_data2
# expression_data[-48] <- scale(expression_data[-48])

# this is the expression data per sample?
kmeans_expr_data <- expression_data %>% 
  filter(symbol %in% deg_genes$symbol) %>% 
  distinct(symbol, .keep_all = TRUE) %>% 
  remove_rownames() %>%
  column_to_rownames("symbol")
# %>% 
# pivot_longer(cols = -symbol, names_to = "samples", values_to = "expression") %>% 
# pivot_wider(names_from = symbol, values_from = expression) %>% 
# column_to_rownames("samples")

set.seed(30921)
c1 <- kmeans(kmeans_expr_data, centers = n_clusters)

kclusterdf <- c1$cluster %>% 
  as.data.frame() %>% 
  rename(kcluster = 1) %>% 
  rownames_to_column("symbol")

metadata <- metadata %>% 
  mutate_at(vars(duplicate, time_point), list(as.factor))
# metadata$time_point <- as.factor(metadata$time_point)

kmeans_plot_data <- kmeans_expr_data %>% 
  rownames_to_column("symbol") %>% 
  pivot_longer(cols = -symbol, names_to = "sample_name", values_to = "expression_value") %>% 
  left_join(kclusterdf, by = "symbol") %>% 
  left_join(metadata 
            # %>% 
            #   rownames_to_column("sample_name_rows")
            ,
            by = c("sample_name" = "sample_name_rows") )



# pff okay just use the established k means and the centered expression to plot per condition.
# because then I at least can plot something
#-----------------------------------------------------------------------------------------
# how many k means?

# i just stole this, no idea https://2-bitbio.com/2017/10/clustering-rnaseq-data-using-k-means.html

scaledata <- kmeans_expr_data

wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata,
                                     centers=i)$withinss)

plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

## make ggplot elbow plot
input.ggplot <- data.frame(ss = wss,
                           clusters = c(1:20))

p <- ggplot(input.ggplot, aes(x = clusters, y = ss)) +
  geom_line(linetype = 2, color = "#310052", size = 0.5) +
  geom_point(size = 4, color = "#25003D") +
  theme_bw() +
  xlab(label = 'Number of clusters') +
  ylab(label = "Within groups sum of squares") 

p

#pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_elbowplot.pdf", sep = ""), useDingbats = FALSE)
#p
#dev.off()

#-----------------------------------------------------------------------------------------
### lets do this

# centered expression data:
#head(expression_data2)
# where are the genes
head(kmeans_expr_data)
dim(kmeans_expr_data)

# the meta data to couple
head(metadata)
dim(metadata)

# kmeans output
head(kmeans)
dim(kmeans)

# make the same
unique(metadata$stimulus)
metadata[metadata$stimulus == "TGFb_TNFa","stimulus"] <- "mix"

#-----------------------------------------------------------------------------------------
### enrichR

# lets see what we got anyways
enrichlist <- lapply(conditions, function(c){
  # select the right condition
  current.set <- kmeans[,c(c, "symbol")]
  
  # perform enrichR per k means
  res <- lapply(1:n.kmeans, function(k){
    # get genes belonging to current k
    set <- current.set[current.set[,1] == k, "symbol"]
    
    # perform enrichR
    dff <- enrichr(set,databases = dbs_mod)
    
    # add stuff for later
    # use the 2018 one anyways
    dff <- dff[["GO_Biological_Process_2018"]]
    #dff <- dff[["GO_Biological_Process_2015"]]
    dff$condition <- c
    dff$k.means <- k
    
    return(dff)
    
  })
  names(res) <- 1:n.kmeans
  return(res)
 
})
names(enrichlist) <- conditions
#View(enrichlist)


## plot per condition
# plot the expression per condition and add the pathway analysis results
plotlist <- lapply(conditions, function(c){
  # select the right condition
  current.set <- kmeans[,c(c, "symbol")] 
  # change name for plotting purposes
  names(current.set) <- c("kcluster", "symbol")
  
  # steal code from arjan
  metadata.plot <- metadata %>% 
    mutate_at(vars(duplicate, time_point), list(as.factor))
  
  means_plot_data <- kmeans_expr_data %>% 
    rownames_to_column("symbol") %>% 
    pivot_longer(cols = -symbol, names_to = "sample_name", values_to = "expression_value") %>% 
    left_join(current.set, by = "symbol") %>% 
    left_join(metadata 
              # %>% 
              #   rownames_to_column("sample_name_rows")
              ,
              by = c("sample_name" = "sample_name_rows") )
  
  means_plot_data <- subset(means_plot_data, stimulus == c) 
  
  p1 <- ggplot(means_plot_data, aes(x = time_point, y = expression_value)) +
    geom_boxplot(outlier.size= 0.5) +
    theme_bw() +
    xlab(label = "Time point") +
    ylab(label = "Expression value")
  
  #p1
  
  p1.1 <- p1 + facet_grid(vars(kcluster), scales = "free") #+ ggtitle(label = paste("condition", c, sep = " "))
  p1.1
  
   ## get all the corresponding enrichR information
  input2 <- lapply(1:n.kmeans, function(k){
    data <- enrichlist[[c]][[k]]
    # select only first 5 colums
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
  p2 <- ggplot(input.ggplot2, aes(x = ratio, y = tidytext::reorder_within(short_term, Adjusted.P.value, k.means), fill = Adjusted.P.value)) + # https://www.geeksforgeeks.org/how-to-reorder-barplots-with-facetting-with-ggplot2-in-r/
    # use reoder_within to sort the bar graph on Term according to p value for the facet k means
    geom_col() + 
    #coord_flip() + 
    theme_bw() +
    scale_fill_continuous(low = "#FFC300", high = "#F06543") +
    #theme(axis.text.y=element_blank()) +
    xlim(0, 0.55) +
    scale_y_discrete(limits=rev) + # reverse the order of y axis plotting
    ylab(label = element_blank()) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
  
  # p2
  
  p2.1 <- p2 + facet_grid(rows = vars(k.means), scales = "free_y") # free y to plot only terms belonging to that K
  #p2.1
  
  #combine plots
  wrapped <-  wrap_plots(p1.1, p2.1) + guide_area() + plot_layout(widths = c(2,1), guides = 'collect') + plot_annotation(title = c)
  #wrapped <-  wrap_plots(p1.1, p2.1) + guide_area() + plot_layout(widths = c(2,1)) + plot_annotation(title = paste(c, "condition", sep = " "))
  
  return(wrapped)
  
})
names(plotlist) <- conditions

#plotlist[["control"]]
# 
# pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_in_vitro_k_means_", n.kmeans, "_and_enrichR.pdf", sep = ""), useDingbats = FALSE, width=11, height=8.5)
# print(plotlist)
# dev.off()
# 
# 
# pdf(paste(result.folder, "/", Sys.Date(), " test2.pdf", sep = ""), useDingbats = FALSE, width=11, height=8.5)
# wrapped
# dev.off()

plotlist2 <- lapply(conditions, function(c){
  # select the right condition
  current.set <- kmeans[,c(c, "symbol")] 
  # change name for plotting purposes
  names(current.set) <- c("kcluster", "symbol")
  
  # steal code from arjan
  metadata.plot <- metadata %>% 
    mutate_at(vars(duplicate, time_point), list(as.factor))
  
  means_plot_data <- kmeans_expr_data %>% 
    rownames_to_column("symbol") %>% 
    pivot_longer(cols = -symbol, names_to = "sample_name", values_to = "expression_value") %>% 
    left_join(current.set, by = "symbol") %>% 
    left_join(metadata 
              # %>% 
              #   rownames_to_column("sample_name_rows")
              ,
              by = c("sample_name" = "sample_name_rows") )
  
  means_plot_data <- subset(means_plot_data, stimulus == c) 
  
  p1 <- ggplot(means_plot_data, aes(x = time_point, y = expression_value)) +
    geom_boxplot(outlier.size= 0.5) +
    theme_bw() +
    xlab(label = "Time point") +
    ylab(label = "Expression value")
  
  #p1
  
  #p1.1 <- p1 + facet_wrap(vars(kcluster), scales = "free") + ggtitle(label = paste("condition", c, sep = " "))
  p1.1 <- p1 + facet_wrap(vars(kcluster), scales = "free") + ggtitle(label = c)
  p1.1
  
  ## get all the corresponding enrichR information
  input2 <- lapply(1:n.kmeans, function(k){
    data <- enrichlist[[c]][[k]]
    # select only first 5 colums
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
  p2 <- ggplot(input.ggplot2, aes(x = ratio, y = tidytext::reorder_within(short_term, Adjusted.P.value, k.means), fill = Adjusted.P.value)) + # https://www.geeksforgeeks.org/how-to-reorder-barplots-with-facetting-with-ggplot2-in-r/
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
  # p2
  
  p2.1 <- p2 + facet_grid(rows = vars(k.means), scales = "free_y") + # free y to plot only terms belonging to that K
    ggtitle(label = c)
    #p2.1
  
  # #combine plots
  # wrapped <-  wrap_plots(p1.1, p2.1) + guide_area() + plot_layout(widths = c(2,1), guides = 'collect') + plot_annotation(title = paste(c, "condition", sep = " "))
  # #wrapped <-  wrap_plots(p1.1, p2.1) + guide_area() + plot_layout(widths = c(2,1)) + plot_annotation(title = paste(c, "condition", sep = " "))
  
  return.list <- list(p1.1, p2.1)
  names(return.list) <- c("k means", "EnrichR")
  
  return(return.list)
  
})
names(plotlist2) <- conditions

## in vitro
# 8-1/4 x 11-3/4 in
pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_in_vitro_k_means_", n.kmeans, "_and_enrichR_seperate.pdf", sep = ""), height =11.45 , width =8.25, useDingbats = FALSE)
print(plotlist2)
dev.off()

# tall pdf
# pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_in_vitro_k_means_", n.kmeans,"_all_conditions_tall.pdf", sep = ""), height =11.45 , width =8.25, useDingbats = FALSE)
# plotlist2[["control"]][["k means"]] + plotlist2[["TGFb"]][["k means"]] +
#   plotlist2[["TNFa"]][["k means"]] + plotlist2[["mix"]][["k means"]]
# dev.off()

# wide pdf
pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_in_vitro_k_means_", n.kmeans,"_all_conditions_wide.pdf", sep = ""), height =8.25 , width =11.45, useDingbats = FALSE)
plotlist2[["control"]][["k means"]] + plotlist2[["TGFb"]][["k means"]] +
  plotlist2[["TNFa"]][["k means"]] + plotlist2[["mix"]][["k means"]]
dev.off()


## enrichR
pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_in_vitro_k_means_", n.kmeans,"_all_conditions_EnrichR_tall.pdf", sep = ""), height =11.45 , width =8.25, useDingbats = FALSE)
plotlist2[["control"]][["EnrichR"]] + plotlist2[["TGFb"]][["EnrichR"]] +
  plotlist2[["TNFa"]][["EnrichR"]] + plotlist2[["mix"]][["EnrichR"]] +
  plot_layout(guides='collect') &
  theme(legend.position='bottom', legend.key.size = unit(0.25, 'cm'), legend.text = element_text(size=6), legend.direction = "vertical")
dev.off()

# pdf(paste(result.folder, "/", Sys.Date(), "_project_TF_in_vitro_k_means_", n.kmeans,"_all_conditions_EnrichR_wide.pdf", sep = ""), height =8.25 , width =11.45, useDingbats = FALSE)
# plotlist2[["control"]][["EnrichR"]] + plotlist2[["TGFb"]][["EnrichR"]] +
#   plotlist2[["TNFa"]][["EnrichR"]] + plotlist2[["mix"]][["EnrichR"]] + plot_layout(guides='collect') &
#   theme(legend.position='bottom')
# dev.off()


save(n.kmeans, conditions, enrichlist, plotlist, plotlist2, file = "In_vitro_9_k_means_enrichR.RData")


## EnrichR table for K means
gettable <- lapply(conditions, function(c){
  
  input2 <- lapply(1:n.kmeans, function(k){
    data <- enrichlist[[c]][[k]][["enrichR"]]
    # get only first 10 pathways
    data <- data[1:10,]
  })
  input.ggplot2 <- do.call(rbind.data.frame, input2)
  
})

Enrich.table <- do.call(rbind.data.frame, gettable)

write.csv2(Enrich.table, "Go_terms_top_10_in_vitro_9_k_means.txt")


# fin



###
# housekeeping genes


input.pheatmap <- as.data.frame(expression_data[expression_data$symbol %in% c("GAPDH", "B2M", "ACTB", "YWHAZ", "HPP1"),])
input.pheatmap <- input.pheatmap[-5,]
rownames(input.pheatmap) <- input.pheatmap$symbol
input.pheatmap <- input.pheatmap[,-length(input.pheatmap)]

pheatmap::pheatmap(input.pheatmap)
pheatmap::pheatmap(input.pheatmap, cluster_cols = F)

#-----------------------------------------------------------------------------------------
### PCA plot

expression_data3 <- as.data.frame(expression_data)
dim(expression_data3)
expression_data3 <- expression_data3[!is.na(expression_data3$symbol),]
dim(expression_data3)
rownames(expression_data3) <- expression_data[,"symbol"]



