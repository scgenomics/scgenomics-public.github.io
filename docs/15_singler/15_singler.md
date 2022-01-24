---
layout: default
title: Cell typing with SingleR
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


<!--  make it look like we had the p_plate etc. by
redoing it here but not showing. Including in rda file is too big
-->





# Cell typing using SingleR

We already have some idea of what each cluster is, we will now subject it
to [SingleR](https://doi.org/10.1038/s41590-018-0276-y)'s automatic cell
type identification procedure.


```r
## load the reference data:
## hpca <- HumanPrimaryCellAtlasData()
file <- paste0(CFG$data_dir, '/hpca.rds') # cached verssion for speed
hpca <- readRDS(file=file)

## check how many (and/or which) types:
length( unique(hpca$label.main) )
length( unique(hpca$label.fine) )

## following takes 2 min or so:
singler_hpca <- SingleR(test = GetAssayData(srat, assay = "RNA", slot = "data"),
   ref = hpca, labels = hpca$label.main)

SingleR::plotScoreHeatmap(singler_hpca,
                          show.labels = TRUE, max.labels=100,
                          show.pruned = FALSE,
                          order.by="clusters",
                          clusters = srat@meta.data[[clus.column]], 
                          annotation_colors = list(Clusters=cluster.colors)
                          )

## show previous plot to get our bearings:
( p_plate | p_patient ) / ( p_type | p_clus )

## do the  asssignment:
  type2 <- rep("unknown", nrow(srat@meta.data))
  names(type2) <- rownames(srat@meta.data)
  
  type2[ WhichCells(srat, idents = 0  )  ] <- 'ESC-0'
  type2[ WhichCells(srat, idents = 1  )  ] <- 'T'
  type2[ WhichCells(srat, idents = 2  )  ] <- 'ESC-2'
  type2[ WhichCells(srat, idents = 3  )  ] <- 'ESC-3'
  type2[ WhichCells(srat, idents = 4  )  ] <- 'NK'
  type2[ WhichCells(srat, idents = 5  )  ] <- 'Mono'
  type2[ WhichCells(srat, idents = 6  )  ] <- 'ESC-6'
  type2[ WhichCells(srat, idents = 7  )  ] <- 'ESC-7'
  type2[ WhichCells(srat, idents = 8  )  ] <- 'B'
  type2[ WhichCells(srat, idents = 9  )  ] <- 'B'
  type2[ WhichCells(srat, idents = 10  )  ] <- 'Mesen'

srat <- AddMetaData(srat, col.name='type2', metadata=type2)

typenames = unique(type2)
### type2.colors <- c(T="red", ESC="black", B="blue", NK="DarkGreen", Mesen="brown", Mono="orange")
type2.colors <- Polychrome::alphabet.colors(length(typenames))
names(type2.colors) <- typenames

p_type2 <- DimPlot(srat, reduction = "umap", pt.size=0.5, 
  group.by='type2', cols=type2.colors)

## show the final result:
( p_patient | p_clus ) / ( p_type | p_type2 )
```

**Challenge**: do the exact same celltyping, now using the
`hpca$label.fine` instead of `hpca$label.main`

<!-- lastly, save the complete sesssion for the next time -->

