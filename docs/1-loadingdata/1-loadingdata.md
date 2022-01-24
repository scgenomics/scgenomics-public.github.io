---
layout: default
title: Loading data
---

<!-- auto TOC from just-the-docs theme -->

<!-- stuff to make Rmarkdown do what we want:  -->



<!-- load complete state from previous lesson: none -->

# Loading data

## Preliminary stuff

Many stages have to make use of common things like \*
libraries-to-be-loaded and files to be `source`-ed \* directories where
to find data, annotations and functions \* parameters such as cutoff
values \* gene, sample and cell annotation and definitions of gene lists
\* color schemes for plots \* etc. etc. You may be tempted to just rely
on `R`'s workspace and saving sessions to keep track of this, but this
easily gets messy and irreproduceable.

For this reason it is useful to have all the global variables in one
place, typically at the top of the script. Here we will use the `CFG`
(an R `list`) for that. It also helps finding back their definition and
actual usage easily in scripts and in Rstudio.

We will define a few variables we will be needing up-front, along with
the libraries we will be needing

### Libraries and functions


```r
library(Seurat) # our main 'environment' 
library(Matrix) # needed for working with sparse matrices
library(dplyr)  # general data manipulation
library(ggplot2) # plotting functions
library(patchwork) # composing multi-panel plots
```

### Config variables


```r
CFG <- list() # global

CFG$data_dir <- "~/data/course/day1/"   # where to read data from
CFG$output_dir <- "~/tmp/day1"              # where to write data to

# explained later
CFG$ndims <- 25                         
CFG$random_seed <- 2033
```

### Annotation


```r
CFG$mito_pattern <- "^MT-" # needed to recognize human mitochondrial genes
```

### Parameters such as cutoffs


```r
# Cutoffs used during quality control for filtering genes (explained later)
CFG$min_txpts <- 5000 
CFG$max_txpts <- 60000
CFG$max_pctmito <- 25

# Graphical parameters

CFG$gradient.colors = viridis::viridis(101, direction= -1)

## NOTE: many colors change if the number of clusters change, so
## they are not fixed here.
```

Loading the data depends on the format of the data. For `10X` data it's
easiest to use the `Seurat::Read10X`. It requires the full name of the
directory `outs/filtered_feature_bc_matrix` which contains the files
`matrix.mtx.gz features.tsv.gz barcodes.tsv.gz`


```r
lib <- "LX010" # our library
dir <- paste0(CFG$data_dir, "/", lib, "/filtered_feature_bc_matrix")
counts <- Read10X(data.dir =dir) # this is of type 'dgCMatrix', a sparse matrix
```

`counts` now contains a so-called *sparse matrix* with all the data.
<!-- @callout: sparse Matrix, also requiring Matrix::colSums() -->

How many cells were loaded? <!-- @challenge: dim(counts) -->

How many txpts are there per cell, on average? Let's also plot the
distribution of the counts to get a 'feel' for the data.


```r
cnts_per_cell <- Matrix::colSums(counts) # note use of Matrix::colSums

summary(cnts_per_cell)

df <- data.frame(counts=cnts_per_cell)

ggplot(df, aes(x=counts)) +
  geom_density() +
  geom_rug() +
  scale_x_continuous(trans='log2')
```

## Changing cell labels

When we are merging data from different sources, identical cell names
from different sources may clash, so we often have to prepend the
sample name to the cell names. Let's do that here too, even though we
are not merging data from different sources.


```r
  cellnames <- colnames(counts)
  colnames(counts) <- paste0(lib, "_", cellnames)
```

## Metadata

Metadata is a fancy term for 'data about the data'. In Seurat it means
'additional data per cell'. The metadata is used very often during an
analysis (e.g. for plotting). Metadata such as the total number of
transcripts per cell (see above) is calculated and added to the big
`SeuratObject` we are about to create. However some useful metadata is
not included automatically (or even known) by Seurat, so we have to do
it ourselves.

For QC purposes, we routinely need: 1. percentage mitochondrial reads 1. log2(1 + transcript counts) 1. log2(1 + gene
counts) 1. library of origin -- identical for all cells in one library

We set up a `data.frame` to hold this information and will use it later.
To get a very quick impression, let's plot the relationship between
number of genes and number of transcripts.

<!-- @callout mitochondrials --->


```r
mitos <- grep(CFG$mito_pattern, rownames(counts), value=TRUE) # search if there are any mitochondrial genes

## show them:
mitos

percent_mito <- 100 * Matrix::colSums(counts[mitos, ]) / cnts_per_cell # note the use of Matrix::colSums again

## get the non-mito's for easy selection of 'the rest of the genes'
nuclear <- setdiff(rownames(counts), mitos)

log2_transcript_counts <- log2(1+Matrix::colSums(counts[nuclear , ]))

log2_feature_counts <- log2(1+Matrix::colSums(counts[nuclear, ] > 0))

## set up a data.frame:
meta <- data.frame(percent_mito = percent_mito,
  log2_counts = log2_transcript_counts,
  log2_features = log2_feature_counts,
  lib=rep(lib, ncol(counts)))

## let's show how features ( = genes ) relate to number of transcripts
ggplot(meta, aes(x=log2_counts, y=log2_features)) + geom_point(color='blue')
```

## Merging data

We currently have only one library, so merging is not needed.

## Creating the Seurat object

Before creating the Seurat object it is good to ignore genes that we are
not interested in, for example the mitochondrial genes.


```r
counts <- counts[ nuclear , ]

srat <- Seurat::CreateSeuratObject(counts=counts,
    project = 'intestinal_organoids',
    meta = meta)

## Later this we won't use the counts object anymore, so 
## we should get rid of it and cleanup
rm(counts)
gc()
```

The Seurat object contains everything needed for the analysis. Although
you often interact with the Seurat object using Seurat's functions (e.g.
`GetAssayData` to get the data) it is useful to have look at the
internal structure of the object.

Typing its name gives a quick summary, using the `str` function on it
exposes the full structure.


```r
srat # just an overview

str(srat)
```

<!-- .callout: slots and list members -->

The toplevel are 13 predefined slots, with `assays` and `meta.data` the
most important ones. `assays` is a `list` with currently just one
member: `RNA`, which in turn has several predefined slots to keep the
actual data (both raw and normalized) and auxiliary data.

For more information on the structure of they Seurat object, see
<https://github.com/satijalab/seurat/wiki/Seurat>.

The other important and often-used slot is `meta.data` which contains
the data.frame we created and supplied to the `CreateSeuratObject`
function. Have a quick look.


```r
head(srat@meta.data) ## short-hand: head(srat)
```

When creating the object, Seurat has added additional columns
`orig.ident`, `nCount_RNA` (the number of transcripts) and
`nFeature_RNA` (the number of genes expressed).

Pretty much all parts of the `srat` object can be changed or extended.
Often this is done by calling a Seurat function on the `srat` object
which returns a new object with changed extended internal parts, e.g.

| `srat <- SCTransform(obj = srat)`

will do the normalization, putting the normalization results into
`srat@assays$SCT@data` and

| `srat <- FindClusters(obj = srat, ...)`

adds a column to `srat@meta.data` containing the cluster number each
cell was assigned to.

<!-- lastly, save the complete session for the next time -->



