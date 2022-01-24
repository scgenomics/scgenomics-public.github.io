---
layout: default
title: Dimensional Reduction
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


# Dimensional Reduction

We now have a set of cells and a set of variable genes to work with,
normalized in two different ways.  However, we cannot make much sense of
cells if each one is described by the 3000 variable genes. It is in fact
a 3000-dimensional mathematical space whose dimensions must be reduced.

## Principal Component Analysis

The first step is to apply Principal Component Analysis. 
The normalized expression counts ( *m* genes x *n* cells ) are replaced by
a lower number of *k* principal components x  *n* cells). 

The principal components (or axes) are new fictitious dimensions
composed of weighted sums of the original gene coordinates.  These new
dimensions are sorted by importance. The first few dimensions
(e.g. 1.. 30) represent most of the relevant actual information. The
remainder is not informative and in fact represents noise.  Throwing out
the higher dimensions reduces noise, but it also makes life much easier
for subsequent steps such as UMAP and clustering.

Let's then run the PCA with a high number of dimensions (we can always
ignore the higher ones).

We have run two different normalizations, each with different sets of
varying features, so the principal component analyses will also be
different.  They are currently two different assays withing the same
object, but there can be only one PCA. Let's therefore copy the
Seurat object and switch it back to use `RNA` as the default assay,
rather than `SCT`. 


```r
## copy and switch it back to the RNA assay:
  srat_rna <- srat
  DefaultAssay(srat_rna) <- 'RNA'
  
  srat_rna <- RunPCA(srat_rna, verbose=TRUE, npcs=50)

## the srat object had the SCT normalization, so we can continue
  srat <- RunPCA(srat, npcs=50)
  
## Show the cells along the first two principal axes using different normalizations,
## again plotting them side by side using patchworks' |-operator:

lognorm <- DimPlot(srat_rna, reduction = "pca",  cols=type.colors, 
      pt.size=1, group.by='type') + labs(title="LogNormalize")

sct <-DimPlot(srat, reduction = "pca", cols=type.colors,
    pt.size=1, group.by='type') + labs(title="SCTransform")
    
lognorm | sct
```

Based on the above you can see that:
1. There is some structure in the data. If there would be just some ellipsoidal blob, there is something wrong
1. The scales are 'flipped'. This is because the signs of axes is arbitrary
1. The scales are different. This due to the different normalization.
1. In the `SCTransform` case, the separation between different celltypes is a bit more orthogonal. 

## Deciding the number of dimensions

How many dimensions to use is again a judgement call, but it's good to
play around. There are many tools you can use to base your choice of
components on.

### Scree plot

The simplest one is the scree plot, which shows you the fraction of the
total variance per principal component. Seurat has something similar
called `ElbowPlot` which shows the standard deviation (not fraction of
variance) per component. The 'knee' in this plot, or alternatively the
component where the slope becomes roughly constant is often a good first
guess.


```r
## The ggplot's patchwork has a nifty way to compose multipanel graphs
## by using '|' for putting things next to each other, and '/' to 
## stack things vertically (think of an arithmetic division line)

lognorm <-  ElbowPlot(srat_rna, reduction="pca", ndims=50) + labs(title="LogNorm") 
sct <- ElbowPlot(srat, ndims=50, reduction="pca") + labs(title="SCTransform")

lognorm | sct

## alternatively : lognorm / sct
```

Looking at these plots it would seem that around 25 principal components 
would seem to be OK, but please play around.

<!-- 
These elbow plots don't show it, but if you look at the fraction of
variance you will see that for this data,
`SCTransform` will put 50% of the total variance in the first 50 PC's,
whereas `LogTransform` only gets to 5% in them. This massive improvement
is a sign that SCTransform is the superior method.
-->

### PC scores heatmap 

The *loadings* of PCA are the weights with which a gene contributes to a
principal component. Let's take, per principal component, the 30 or so
genes having the most positive and most negative loadings and the 100 or
so cells with the strongest scores in these principal axes. If you plot
their scores as a heatmap it shows, per principal component, the origin
and strength of the primary discriminants of the differences between the
cells.  This is done using the `DimHeatmap` function.

Expect strong differences between relevant genes in the first few
components, with the picture getting blurier for the higher components.
Ignore the 'sign' of the scores as they are arbitrary. 


```r
##
DimHeatmap(srat_rna, dims = 1:12, cells= 100)

DimHeatmap(srat, dims = 1:12, cells= 100)
```

**Challenge:** Show the behaviour in the higher Principal Compenents.


<!-- Discuss:

@@ some of the genes, esp. muc2-mneongreen

-->

Again, around 25 dimensions would seem to be enough. 

## UMAP: Uniform Manifold Approximation and Projection

Based on the above, we could settle for e.g. 25 dimensions but this is
still way too much to visualize. We need to compress these 25
dimensions further to just 2 dimensions in order to plot them. This has
to be done in such a way that, by and large, the similarities and
dissimilarities of cells (i.e. the distances between the cells in the
plot) are preserved after the compression.  UMAP (and previously t-SNE)
is one of the most frequently used methods to achieve this.

### Random seeds

The UMAP method itself is fairly involved and it is also slightly
non-deterministic, in that the computation has to start with a random
seed, which by default is different every time you run it.  You can fix
this random seed by a call to `set.seed(some_fixed_integer)` before the
call to `RunUMAP`. Try out a few different values for the seed and see
the difference in the resulting map.

### Neighbours

As a first step in the dimensionality reduction, UMAP creates local
'maps' of the original data. They are based on the nearest neighbours of
each point. By specifying the number of neighours using the
`n.neighhors` parameters you can change how many neighhors UMAP will
take into account. This will determined how locally or globally UMAP will
preserve the overall structure of the data. Using a low `n.neighbors`
(say, 10) it will only preserve very local structure, whereas using a
high `n.neighbors` (e.g., 50) it will give more emphasis to preserving
the global structure.  Play around with this parameter.

Be aware that in addition to the random seed and other UMAP parameters,
the dimensional reduction is of course also dependent on the data
itself. UMAP is very sensitive to this: even small changes in the data
such as adding or removing a few genes from the list of variable genes
can completely alter the layout. This is perfectly normal (also: there
is no 'correct layout') but during the analysis of a project it is often
better to fixate the layout by sticking to a fixed seed, saving the
seurat object or the UMAP coordinates to disk and reading them back in
manually.

## Running UMAP

Remember that we have two different normalizations, so we have to run
two separate UMAP's too. 


```r
## CFG$ndims <- 25
## CFG$random_seed <- 2033

set.seed(CFG$random_seed)

srat_rna <- RunUMAP(srat_rna, dims=1:CFG$ndims )

srat <- RunUMAP(srat, dims=1:CFG$ndims )

lognorm <- DimPlot(srat_rna, reduction = "umap", pt.size=1, 
  group.by='type', cols=type.colors) +
  labs(title='lognorm') 
  
sct <-  DimPlot(srat, reduction = "umap", pt.size=1, 
  group.by='type', cols=type.colors) +
  labs(title='SCT')

## show next to each other:
lognorm  | sct

## after this we will leave out the 'reduction = "umap"' for brevity,
## DimPlot first searches for umap, then tsne, then pca
```

If you look closely you can see that `LogTransform` splits both the
Goblet cells and the Secretory Progenitor cells into two separate
clusters. Based on this we could conclude, again, that `SCTransform` is
superior. 

<!-- make this an exercise -->

Let's play with the  `n.neighbors` parameter to see how things change.


```r
### play around with the n.neighbors parameter:
srat_ngbs5 <- RunUMAP(srat, dims=1:CFG$ndims, n.neighbors=5)

ngbs30 <- DimPlot(srat, pt.size=1, group.by='type', cols=type.colors )+
  labs(title='SCT, n.ngbs=30') 

ngbs5 <- DimPlot(srat_ngbs5, pt.size=1, group.by='type', cols=type.colors) +
  labs(title='SCT, n.ngbs=5')

ngbs30 | ngbs5

## we don't need srat_ngbs5 after this so clean it up:
rm(srat_ngbs5)
gc()
```

The default number of 30 neighbours works well for many cases so we
will keep using that here. For smaller datasets it can help to reduce
the number of neighbours to ~ square root of the number of cells 
([Kobak 2019 &amp; Berens](https://pubmed.ncbi.nlm.nih.gov/31780648/)).


## Other methods

There are other methods of dimensionality reduction avaiable, such as
t-SNE.  Depending on your data they may bring out different aspects of
your data more prominently.  You can also load your own dimensionality
reduction from outside Seurat.  We will not go in to that here.

**Challenge**: try out tSNE for these data.

**Challenge**: if you're adventurous you can try embedding the data in 3 dimensions
  rather than two.

<!-- lastly, save the complete sesssion for the next time -->

