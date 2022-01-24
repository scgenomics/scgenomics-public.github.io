---
layout: default
title: Differential Gene Expression
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


# Differential gene expression

A good way to establish or confirm the celltypes of a cluster is to have
Seurat find, for each cluster we have, the genes that are differentially
expressed in that cluster as compared to all the cells not in that
cluster.

The function `FindAllMarkers` does exactly that. There is one important
caveat: `FindAllMarkers` and similar functions (currently) work better
by using the `LogNormalized` data. We 
therefore explicitly pass
`assay='RNA'` as an argument to `FindAllMarkers`. 

**Warning**: make sure the `RNA` assay is present and log-normalized and
scaled. If it is not, `FindAllMarkers` will use the raw data. If your
data was not yet log-normalized, explicitly set the `DefaultAssay` to
`RNA` first, otherwise it will overwrite the `SCT` assay!

`FindAllMarkers` returns a `data.frame` with p-value, multiple-testing
corrected p-value (`p_val_adj`), average log(fold change), and the
fraction (not percentage) of cells where the feature is detected in the
first and second group.

Often this analysis is restricted to look only for up-regulation.  To
skip marginal results and to speed things up, genes not detected in
fewer than `min.pct` fraction of the cells in either population or
differing less than `logfc.threshold` in fold-change on average are
never tested.

The function takes fairly long so it is a good idea to save the results 
straight away.

We will then plot the top 5 of most significant genes per cluster and
a subsample of cells of that cluster, as a heatmap. 

An alternative representation is the DotPlot, which we also show.


```r
column <- 'SCT_snn_res.1.2'
Idents(srat) <-  column

clusters <-sort(unique(srat@meta.data[[column]]))
cluster.colors <- Polychrome::alphabet.colors(length(clusters))
names(cluster.colors)=as.character(clusters)

de_genes <- FindAllMarkers(srat, assay='RNA', 
                           only.pos = TRUE,
                           min.pct = 0.1,
                           logfc.threshold = log(1.5) )

## save just in case
file <- paste0(CFG$output_dir, "/intestorgs_de_genes.txt")
write.table(de_genes, file=file, sep = "\t", 
    quote = FALSE, row.names = T, col.names = NA, na = "")

de_genes %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5

DoHeatmap(subset(srat, downsample = 50),
          features = top5$gene,
          group.colors=cluster.colors,
          assay='RNA',
          slot='scale.data')  + theme(axis.text.y=element_text(size=6))

## Another way of showing this:
DotPlot(srat, features=unique(top5$gene), cols=c('blue', 'red'), 
   col.min=0, col.max=1, dot.scale=3) + 
   coord_flip()  # genes are rows, columns are clusters
## If the labels are unreadable, 'add' the following to the plot function invocation
## + theme(axis.text.y=element_text(size=6))

## See the marker gene exercise for 'markerplot' object with the overview of the markers
```

<!-- markerplot, see 6-markergenes/6-markergenes.Rmd, was lost , will not recreate it here -->

## Cell type assignment

As you can see in the heatmap, some clusters (such as the
entero-endocrine cells) are clear-cut, whereas other clusters share
their top-5 genes with other clusters.

We will now assign the 'definitive' cell types, based on the expression
of the marker genes and the automatic clustering. 

Since you may have made different choices regarding the cutoffs,
confounders and random seeds it is not guaranteed that you will see the
exact same clusters, colors and layout that our example code has
produced. 

We therefore first plot the way things look like with our choices, and
proceed from there.



```r
p_bytype <- DimPlot(srat,
  pt.size = 1, group.by = "type", cols = type.colors) + labs(title = "types")

p_byphase <-  DimPlot(srat, pt.size = 1, group.by = "Phase")

column <- "SCT_snn_res.1.2"
clusters <-sort(unique(srat@meta.data[[column]]))
cluster.colors <- Polychrome::alphabet.colors(length(clusters))
names(cluster.colors) <- as.character(clusters)

p_bycluster <-  DimPlot(srat, pt.size=1, group.by=column, cols=cluster.colors ) +  labs(title=column)

##@ show side by side

p_bytype | p_byphase  | p_bycluster 
```

It is perhaps logical to start with the stem cells. 

### cluster 1

The SMOC2 stem cell marker is most clearly expressed in cluster 1. Let's assign
all of cluster 1 to *intestinal stem cells* (ISC).

### clusters 4,6,8,9

The stem cell cluster branches off into two 'lobes'. The right lobe,
consisting of clusters 4,6,8 and 9 can be seen to be cycling in the
Phase plot. The upper half of this combined cluster has high expression
of the MKI67 proliferation marker.

Let's lump them together into *transit amplifying cells* (TA)

### cluster 0,3,5

Cluster  and  have clear expression of early enterocyte marker FABP1, 
and cluster 3 too but less so. Cluster 3, the 'tip' of this lobe 
has high expression of mature enterocyte marker KRT20.  Lets's therefor call
cluster 0 and 5 *early enterocytes* (EEC) and cluster 3 *mature enterocyte* (EC).

### cluster 7

HES6, a known marker of secretory progenitors coincides with cluster 7.

However, based on  the expression of Tuft cell marker AVIL, a few of those
cells are probably Tuft cells. More on this below.

### cluster 10

CHGA, a marker of enteroendocrine cells is very speciffically expressed only 
in cluster 10, making it an obvious assignment.

### cluster 2

Cluster 2 has clear expression of Goblet marker MUC2, in our case
represented by the muc2-mneongreen construct. However, using our Paneth1
cell module score, it is also clear that the tip of that part of the
UMAP is probably Paneth instead of Goblet. 

## Assigning the types

Our final judgement is:

| cluster 1: stem cell (ISC)
| cluster 4,6,8,9: transit amplifying (TA)
| cluster 0,5: early enterocyte (early EC)
| cluster 3: enterocyte (EC)
| cluster 7: secretory progenitor, high AVIL expression: Tuft
| cluster 10: enteroendocrine (EEC)
| cluster 2: Goblet, partly Paneth, high Paneth module score: Paneth

Let's create a new meta.data column 'type2' to hold these assignments.
We can use the `WhichCells` function to select the cells.  For the Tuft
and Paneth cells we will select on both the cluster number (the `idents`
argument) as well as the normalized expression (`expression` argument).


```r
## set up empty vector with the names, then fill them.
types2 <- rep("unknown", nrow(srat@meta.data))
names(types2) <- rownames(srat@meta.data)

## do the assignment
types2 [ WhichCells(srat, idents=1) ] <- 'ISC'
types2 [ WhichCells(srat, idents=c(4,6,8,9) ) ] <- 'TA'
types2 [ WhichCells(srat, idents=c(0,5) ) ] <- 'early EC'
types2 [ WhichCells(srat, idents=3) ] <- 'EC'
types2 [ WhichCells(srat, idents=7) ] <- 'secretory prog.'
## override those with high AVIL expression:
types2 [ WhichCells(srat, idents=7, expression = { AVIL > 2 } ) ] <- 'Tuft'
types2 [ WhichCells(srat, idents=10) ] <- 'EEC'
types2 [ WhichCells(srat, idents=2) ] <- 'Goblet'
## override those with high Paneth expression
types2 [ WhichCells(srat, idents=2, expression =  { Paneth1 > 1 }) ] <- 'Paneth'

srat <- AddMetaData(srat, col.name='type2', metadata=types2)

## show it:

type2.colors <- Polychrome::dark.colors( length(unique(types2)))
names(type2.colors)<-unique(types2)

p_bytype2 <- DimPlot(srat,
  pt.size = 1, group.by = "type2", cols = type2.colors) + labs(title = "type2")

( p_bytype | p_byphase ) / (  p_bycluster  | p_bytype2 )
```

<!-- lastly, save the complete sesssion for the next time -->


