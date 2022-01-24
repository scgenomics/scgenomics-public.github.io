---
layout: default
title: Clustering
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->






# Clustering

The next step is clustering. Try a few different resolutions.


```r
srat <- FindNeighbors(srat, dims = 1:CFG$ndims)

srat <- FindClusters(srat, resolution = 0.4, algorithm = 1)

srat <- FindClusters(srat, resolution = 0.8, algorithm = 1)

srat <- FindClusters(srat, resolution = 1.2, algorithm = 1)

columns <- grep("SCT_snn_res", names(srat@meta.data), value = TRUE)

##@ show numbers of cells per cluster (top row shows the
##@ cluster ids)
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(srat@meta.data[[col]]))
}
```

<!-- note: you can also pass an array!  -->


## Choosing the clustering resolution


```r
## let's plot the differnt clusterings side by side
## to see the differeces:
p_clusbyreso <- list()

for (reso in c("0.4", "0.8", "1.2")) {

  clus.column <- paste0("SCT_snn_res.", reso)
  clusters <- sort(unique(srat@meta.data[[clus.column]]))

  ## switching color scheme to make it look different:
  clus.colors <- Polychrome::alphabet.colors(length(clusters))
  ## clus.colors <-
  ## Polychrome::kelly.colors(length(clusters))
  names(clus.colors) <- as.character(clusters)

  p <- DimPlot(srat, pt.size = 1, group.by = clus.column, cols = clus.colors) +
    labs(title = reso)
  p_clusbyreso[[reso]] <- p
}

p_clusbyreso[[1]] |  p_clusbyreso[[2]] | p_clusbyreso[[3]] 
##=> 0.8 is fine?

## set the identies to the clusters coming from the 0.8 resolution:

  clus.column <- paste0("SCT_snn_res.", 0.8)
  Idents(srat) <- clus.column
  
  clusters <- sort(unique(srat@meta.data[[clus.column]]))
  ## switching color scheme to make it look different:
  cluster.colors <- Polychrome::alphabet.colors(length(clusters))
  names(cluster.colors) <- as.character(clusters)

##@ and show result with the other umaps to get your bearings:
  p_clus <- DimPlot(srat, reduction = "umap", pt.size=0.5, 
     group.by=clus.column, cols=cluster.colors)

( p_plate | p_patient ) / ( p_type | p_clus )

rm(p_clusbyreso)
gc()
```

## Using clustree

Although the best clustering has always to be informed by the biology,
we often need a first guess.  As with PCA, there are several 'rules'
to choose the 'best' resolution for the clustering.

The `clustree` library may be able to help decide this.  It depicts the
'splitting' of clusters as the resolution increases, with
interconnections where cells jump between the main clusters.  A
reasonable resolution is one were the clusters stop splitting up and
there is little 'jumping around' for a few resolutions steps.

For this it is useful to create a larger range of resolutions.
We can then plot them using the `clustree` function. 


```r
resolutions <- seq(0.4, 2.0, 0.2)

for (reso in resolutions  )
    srat <- FindClusters(srat, resolution = reso, algorithm = 1)

clustree::clustree(srat, prefix = "SCT_snn_res.", node_size=2, edge_width=0.5)

clus.column <- paste0("SCT_snn_res.", 0.8)
Idents(srat) <- clus.column
```

**Challenge**: Given these splits in the clustree diagram, see if you can find meaningful gene expression differences that explain the splits.

<!-- lastly, save the complete sesssion for the next time -->

