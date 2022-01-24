---
layout: default
title: Normalization and Dimensional Reduction
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


# Normalization and Dimensional Reduction
## Normalization


```r
##@ first LogNormalize, is needed anyway

srat <- NormalizeData(srat, normalization.method = "LogNormalize")

srat <- ScaleData(srat, features = rownames(srat), verbose = FALSE)

srat <- FindVariableFeatures(srat)

##@ then SCTransform:
srat <- suppressWarnings(SCTransform(srat, vars.to.regres = NULL,
  verbose = FALSE))
```

## Dimensional Reduction

As before, we will reveal (this time partial) type  information of these cells to you, 
to make it a bit more interesting. 


```r
## read the externally provided type information and add it to the object.
file <- paste0(CFG$data_dir, "/rms-types.rds")
types <- readRDS(file)
srat <- AddMetaData(srat, col.name='type', metadata=types)

## run the PCA
srat <- RunPCA(srat, verbose = TRUE, npcs = 50)

  DimPlot(srat,
    reduction = "pca", pt.size = 1, group.by = "plate_id")  |
  DimPlot(srat,
    reduction = "pca", pt.size = 1, group.by = "patient")  |
  DimPlot(srat,
    reduction = "pca", pt.size = 1, group.by = "type")
```

## Choosing the dimensions


```r
ElbowPlot(srat, reduction="pca", ndims=50)
##=> around 30?

DimHeatmap(srat, dims = (1:24), cells=100)

## and the next 24:
DimHeatmap(srat, dims = 24 + (1:24), cells=100)
##=> also around 30

CFG$ndims=30
```

## UMAP

Let's do the UMAP and check cell cycle, percentage mitochondria, hemo's and stress 
if they are might be confouding.


```r
data(cc.genes)

srat<- CellCycleScoring(object = srat,
        s.features =   cc.genes$s.genes,
        g2m.features = cc.genes$g2m.genes)

stress_genes <- readRDS(paste0(CFG$data_dir, '/stress_genes.rds'))
srat <- AddModuleScore(srat, features=list(stress=stress_genes), name='stress')

set.seed(CFG$random_seed)

srat <- RunUMAP(srat, dims=1:CFG$ndims )

## show the results:
p_plate <- DimPlot(srat, reduction = "umap", pt.size=0.5,
  group.by='plate_id')

p_patient <- DimPlot(srat, reduction = "umap", pt.size=0.5,
  group.by='patient')

p_type <- DimPlot(srat, reduction = "umap",pt.size=0.5,
  group.by='type')

p_phase <- DimPlot(srat, reduction = "umap", pt.size=0.5,
  group.by='Phase')

( p_plate | p_patient ) / ( p_type | p_phase ) 

## and also the continue variables:
##@ first need colors
CFG$gradient.colors = viridis::viridis(101, direction= -1)

FeaturePlot(srat, pt.size=0.5, 
  feature=c('nCount_RNA', 'percent_mito',
  'pct_hemo', 'stress1'), order=TRUE,
  cols=CFG$gradient.colors)
```

**Challenge**: See if ribosomal protein coding genes are a
confounder, see file `ribo_genes.rds` in `CFG$data_dir`.


<!-- lastly, save the complete sesssion for the next time -->

