################################################################################
#                                PACKAGES TO LOAD                              #
################################################################################

cat("\n* General packages...\n")
install.packages.auto("credentials")
library("credentials")
credentials::set_github_pat()

install.packages.auto("readr")
install.packages.auto("optparse")
install.packages.auto("tools")
install.packages.auto("dplyr")
install.packages.auto("tidyr")
install.packages.auto("tidylog")
library("tidylog", warn.conflicts = FALSE)
install.packages.auto("naniar")

# To get 'data.table' with 'fwrite' to be able to directly write gzipped-files
# Ref: https://stackoverflow.com/questions/42788401/is-possible-to-use-fwrite-from-data-table-with-gzfile
# install.packages("data.table", repos = "https://Rdatatable.gitlab.io/data.table")
library(data.table)

install.packages.auto("tidyverse")
install.packages.auto("knitr")
install.packages.auto("DT")

# for plotting
install.packages.auto("qqman")
install.packages.auto("forestplot")
install.packages.auto("pheatmap")

# for meta-analysis
install.packages.auto("meta")
install.packages.auto("bacon")

install.packages.auto("reshape2")

install.packages.auto("ggpubr")
install.packages.auto("patchwork")
install.packages.auto("corrr")

# Installation of ggcorrplot()
# --------------------------------
if(!require(devtools))
  install.packages.auto("devtools")
devtools::install_github("kassambara/ggcorrplot")

library(ggcorrplot)
install.packages.auto("PerformanceAnalytics")
install.packages.auto("GGally")
library(GGally)

install.packages.auto("haven")
install.packages.auto("tableone")

install.packages.auto("survival")
install.packages.auto("survminer")
install.packages.auto("Hmisc")

# Install the devtools package from Hadley Wickham
install.packages.auto('devtools')

cat("\n* Genomic packages...\n")
install.packages.auto("GenomicFeatures")
install.packages.auto("GenomicRanges")
install.packages.auto("SummarizedExperiment")
install.packages.auto("DESeq2")
install.packages.auto("org.Hs.eg.db")
install.packages.auto("mygene")
install.packages.auto("TxDb.Hsapiens.UCSC.hg19.knownGene")
install.packages.auto("org.Hs.eg.db")
install.packages.auto("AnnotationDbi")
install.packages.auto("EnsDb.Hsapiens.v86")
install.packages.auto("EnhancedVolcano")

# Install the annotation tables
library("devtools")
devtools::install_github("stephenturner/annotables")
library(dplyr)
library(annotables)

# alternative chart of a correlation matrix
# --------------------------------
# Alternative solution https://www.r-graph-gallery.com/199-correlation-matrix-with-ggally.html
install.packages.auto("GGally")

# Quick display of two cabapilities of GGally, to assess the distribution and correlation of variables 
library(GGally)



