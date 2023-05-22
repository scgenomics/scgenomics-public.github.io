## this file is sourced by all the R files so that knitr can be done per lesson

## used often so we want call without prefixes: 
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggraph)
library(patchwork)
library(clusterProfiler) # enrichment analysis
library(org.Hs.eg.db) # gene annotation databases
library(SCutils) # our R package with various convenience data and functions.
library(SingleR) # automatic celltype identification
library(infercnv) # finding CNV's in scRNAseq data
#### used rarely, don't load but use prefixes explicitly
## library(ggExtra)
## library(Polychrome)
## library(viridis)
## library(clustree)
## 
#### maybe later:
## library(SeuratWrappers) # e.g. for velocity
## library(ggthemes)
## library(ggpubr)
## library(clusterProfiler)
## library(zeallot)
## library(zoo)
## library(uuutils)
## library(SingleR)
## library(celldex)
## library(DBI)
## library(topGO)
## library(org.Hs.eg.db)
