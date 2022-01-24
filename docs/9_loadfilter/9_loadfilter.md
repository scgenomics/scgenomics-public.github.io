---
layout: default
title: Loading and filtering
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson ; not the first time-->

# Loading and filtering data
## Preliminary stuff


```r
# libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clusterProfiler) # enrichment analysis
library(org.Hs.eg.db) # gene annotation databases
library(SingleR) # automatic celltype identification
library(infercnv) # finding CNV's in scRNAseq data
library(clustree) # for visualizing how clusters split when increasing resolution
library(ggraph) # needed for clustree
```


```r
# libraries
CFG <- list() # global

CFG$data_dir <- "/opt/course/day2" # where to ready data from
CFG$output_dir <- "."              # where to write data to

CFG$random_seed <- 2033
```

## Loading data

This time we provide you with a ready made Seurat object. Empty wells,
and control features have already been removed, but otherwise the
data is raw. 


```r
## Load it:
file <- paste0(CFG$data_dir,"/srat-raw.rds")
srat <- readRDS(file=paste0(CFG$data_dir,"/srat-raw.rds"))
```

## Filter cells

The next step, a QC-step really, is to remove low quality cells.

## Minimum transcripts


```r
## show the distribution of transcipts, by plate:
v <- VlnPlot(srat, "nCount_RNA", group.by = 'plate_id')

v

## same, but using a logarithmic scale:
v + scale_y_continuous(trans='log2')

## following might be reasonable:

CFG$min_txpts = 1000
CFG$max_txpts = Inf  # meaning Infinity
```

## Max transcripts

Some of the plates have cells that might seem outliers, but they are
be different per plate. You could get rid of any cell with e.g. >=
50,000 transcripts, or you throw out the cells in the top 0.5% of
transcript counts.  They turn out not to matter much so we will not do
that this time.

<!-- for code throwing out 0.5% top percentile per plate, see prior to 21-Jan-2022 07:19:48 -->

## Mitochondrial percentage

<!-- @@ say something about myotic cells naturally expressing high mito -->


```r
v <- VlnPlot(srat, "percent_mito", group.by = 'plate_id')

v

CFG$max_pctmito = 60
```

## Final selection


```r
## show our choices:
FeatureScatter(srat, feature1='nCount_RNA', feature2='percent_mito', pt.size=2)  + 
  geom_hline(yintercept = CFG$max_pctmito, linetype=2) +
  geom_vline(xintercept = c(CFG$min_txpts, CFG$max_txpts),linetype=2) + 
  scale_x_continuous(trans='log2')

## how much do we have currently:
dim(srat)

## what will we loose:
dim(subset(srat, nCount_RNA < CFG$min_txpts))
dim(subset(srat, percent_mito > CFG$max_pctmito))

## and do the subsetting
srat <- subset(srat,
               subset =
                 nCount_RNA >= CFG$min_txpts
               & percent_mito <= CFG$max_pctmito ) 

## how much do we have now:
dim(srat)

# maybe save
file  <- paste0(CFG$output_dir,'/srat-cellsfiltered.rds')
saveRDS(file=file, srat)
```



<!-- lastly, save the complete sesssion for the next time -->

