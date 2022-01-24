---
layout: default
title: Removing confounders
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->



# Removing confounders

You can remove the usual confounders if you like (see the intestinal
organoid lessons), but this time that is not enough.  Within in the main
clusters, there are cell cycle-induced sub-clusters.

The reason is that although the canonical cell cycle genes are gone,
there are quite a few additional ones, and they
drive to UMAP differences.  We can find them by their correlation with
the cell cycle ModuleScores.  We give you the function
`derivedCellcycleGenes` to find them. This function is in file
`correlators.R` in your data directory which you have to `source()`.

It returns a `list()` with two vectors of gene names. Both sets 
have to be removed from the variable genes.



```r
data(cc.genes)

file <- paste0(CFG$data_dir, "/stress_genes.rds")
stress_genes <- readRDS(file)

file <- paste0(CFG$data_dir, "/hemo_genes.rds")
hemo_genes <- readRDS(file)

## get the external code:
source(paste0(CFG$data_dir, "/correlators.R"))

## this function plots an overview but also returns a list
additional <- derivedCellcycleGenes(srat)

## see what's in the list:
show(additional)

## combine these sets (also including stress_genes and hemoglobin genes)
remove <- unique(c(cc.genes$s.genes, cc.genes$g2m.genes,
  additional$S.derived, additional$G2M.derived, 
  stress_genes, hemo_genes))

## check how many:
length(intersect(VariableFeatures(srat) , remove))
```

## Redo the Dimensional Reduction

We can now throw out the confounders. 
After this, we have to do the complete dimensional reducion over again!


```r
## do the removal:
VariableFeatures(srat) <- setdiff( VariableFeatures(srat), remove)

## how much now:
length(VariableFeatures(srat))

## IMPORTANT: the VariableFeature are kept per assay,
## but we only removed them for the DefaultAssay, which
## is 'SCT'. We also must remove them for the RNA assay
DefaultAssay(srat) <- 'RNA'
VariableFeatures(srat) <- setdiff( VariableFeatures(srat), remove)
DefaultAssay(srat) <- 'SCT'

srat <- RunPCA(srat, verbose = TRUE, npcs = 50)

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

## and also the continuous variables:
FeaturePlot(srat, pt.size=0.5, 
  feature=c('nCount_RNA', 'percent_mito',
  'pct_hemo', 'stress1'), order=TRUE,
  cols=CFG$gradient.colors)
```

Compare with the previous umaps to see that we indeed got rid of confounding.

**Challenge**: Make side-by-side UMAP plots before and after deconfounding.

**Challenge**: Is sex likely to be a confounder for these samples? 
You can use the `XIST` and `TSIX`  genes to find out.

<!-- lastly, save the complete sesssion for the next time -->

