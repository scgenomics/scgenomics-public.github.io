---
layout: default
title: Normalization
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


# Normalization

Normalization is the act of making the data look like it was measured
accurately, correcting as much as possible the measurement errors
introduced by differences in:

1. reagent concentrations in the droplets
1. efficiency of the amplifications inside the droplet and during library prep
1. temperature during handling
1. etc. etc.

We will discuss two methods: `LogNormalize` and `SCTransform`. Results
of both normalization methods are stored separately inside the Seurat
object, so we can easily compare them.

## Variable Features: informative variability

In connection with normalization we also have to introduce the notion of
'variable features'.  Features (in our case: genes) whose expression is
more or less constant acrosss all cells are not useful to establish
differences between cells. They are best ignored because they increase
the computational load while mostly contributing noise.  We therefore
have to decide which genes are 'most varying'. But more highly expressed
genes automatically have a higher standard deviation, so the measure of
variability needs to be scaled somehow (hence the connection with
normalization). Having determined this 'scaled variability' and ranked
the genes accordingly, the top few thousand genes are marked as being
variable for use use in subsequent steps.  `LogNormalize` does this in a
seperate step; `SCTransform` has it built-in.

## DefaultAssay

A word here about the so-called `DefaultAssay`.  This is a mechanism in
Seurat to simplify dealing with various data sources.  Up till now we
have been dealing with the `RNA` assay, but this is about to change.

We mention the `DefaultAssay` here because this mechanism is used to
keep the the `LogNormalize` and `SCTransform` results separate.  To
[quote one of the
developers](https://github.com/satijalab/seurat/issues/3043#issuecomment-632754495)
"[The `defaultAssay` is] the currently active assay in the Seurat
object. The Seurat object can store multiple independent assays, with
the requirement that the same cells are present across assays. This can
be used for different experimental methods that measure multiple data
modalities per cell (eg, CITE-seq measures RNA and protein, SNARE-seq
measures RNA and ATAC), or storing different processed forms of the data
(eg, integrated data, log-normalized data, SCTransform-normalized data)
within a single object."  For more information see
[https://github.com/satijalab/seurat/wiki/Seurat](https://github.com/satijalab/seurat/wiki/Seurat).


Function `DefaultAssay` can be used to get and set the currently
active assay.  Some methods, including `SCTransform`, set the
`DefaultAssay` automatically. 

## LogTransform

The easiest way to normalize is to pretend that each cell originally had
exactly the same total number of transcripts (by default, this number is
10,000 transcripts, but the number is unimportant).

The scaling the transcripts of gene *g* is simply

&nbsp;&nbsp;&nbsp;&nbsp;norm(g) = 10,000 * g / ( sum of all g )

This number is subsequently logged (hence *`Log`*`Normalize`) since the
log scale works better for expression data.

An equal number of transcripts is not realistic for samples with wildly
different celltypes. Although the method is now superceded by
`SCTransform` (see below), it works reasonably well. 
We recommend to always run `LogTransform` regardless because it is 
required when determininig differential gene expression between groups
of cells. 

<!--  also because srat@assays$RNA@data is otherwise just raw! 
and also because doing later without setting DefaultAssay = RNA
will overwrite the SCTransform assay  $%^&*()
-->

### Additional steps: scaling and finding variable genes

`LogNormalize` needs to be followed by `ScaleData` and `FindVariableFeatures`.

`ScaleData` Z-transforms the normalized data. That is, per gene, the
number of standard deviations away from the mean is calculated. This is
needed for the next steps, dimensional reduction.  The reason to do
normalization and scaling in separate steps is that the latter can also
be used to regress out (i.e. correct for) the influence of e.g. cell
cycle or percentage mitochondrial reads. We advise against this
regressing-out as it often makes the data worse.

For now we simply rely on the defaults for the `LogNormalize` normalization.
Because we want to compare the behaviour `LogNormalize` and `SCTransform`
we will use different Seurat objects for both.

### Cell types

Before we normalize we are going to give you the cell types
of the current set. The cells are:

| Paneth cells - secret antimicrobial peptides
| goblet cells - secret mucus
| Tuft cells - chemosensing/immunomodulatory cytokine-secreting 
| EC:  enterocytes -  the intestinal absorptive cells
| EEC: enteroendocrine cells -  hormone secreting cells
| ISC: intestinal stem cells
| TA: transient-amplifying 
| NA  not asssigned

You will later infer that yourself, but for now it makes it easier to
judge the normalizations as well as the dimensionality reductions.
The celltypes as determined by Jeff are available from the `celltypes.rds`
file. 

We will also switch to colors that are more easily discernable than
Seurat's default colors.

<!-- NOTE: colors could be put in the CFG list but they are
determined dynamically and change too often, so I don't. -->


```r
## load the cell types from disk:

file <- paste0(CFG$data_dir, '/celltypes.rds')
celltypes <- readRDS(file)
srat <- AddMetaData(srat, col.name='type', metadata=celltypes)

types <- sort(unique(srat@meta.data$type))
type.colors <- Polychrome::dark.colors( length(types ))
names(type.colors) <- types
## use as DimPlot( ... , cols=type.colors, )
```

Next, do the normalization:


```r
srat <- NormalizeData(srat,normalization.method = "LogNormalize")

srat <- ScaleData(srat, features = rownames(srat), verbose=FALSE)

srat <- FindVariableFeatures(srat)
```

The results of the `LogNormalize` normalization go the slot `data` and of the scaling to
`scale.data`, in this case of the `RNA` assay. If needed you can
retrieve them with using the `GetAssayData` function.

<!-- 

There are a number of 'slots':

 * `counts`: the raw data
 * `data`: log(1+normalized(counts)), often used for visualisation color
 * `scale.data`: Z-scores of `data` with any covariates regressed out, used for

The `SCT` assay uses similar slots with the same names

-->

## SCTransform

As explained in the lecture, the normalization depends on the number of
transcripts per cell in a non-linear way.  The `SCTransform`-method
takes this into account, and takes care of all the details.  By
specifying a `vars.to.regresss` argument it can also 'regress out'
(i.e. correct for) other variables such as mitochondrial content. As
already mentioned, we don't use this option as it tends to make 
the data worse.

The `SCTransform`-normalized data gives rise to a new assay: `SCT`.
That is, the normalization results are stored in a different part of the
Seurat object, and the `DefaultAssay` is automatically set to `SCT` (see
above). The previous assay, `RNA` is still present in the object, and
can be 'reactivated' if needed by doing `DefaultAssay <- 'RNA'`.
You can often pass the assay to be used explicitly as an argument to
many functions calls. 

`SCTransform` also scales the data and sets the variable features, so
unlike `LogNormalize`, you should not call `ScaleData` and
`FindVariableFeatures` afterwards.


```r
## SCTransfrom is very chatty, so we wrap it in a supressWarnings() call:
srat <- suppressWarnings(SCTransform(srat,
               vars.to.regres=NULL, 
               verbose=FALSE,
               variable.features.n=3000))

## note that the assay has now changed:
DefaultAssay(srat)

## Since SCTransform normalization takes a bit longer, it may be useful to store the
## resulting object:

## file=paste0(CFG$output_dir, "/srat-normalized.rds")
## saveRDS(file=file,srat)
## To read the object back in, do
## srat <- readRDS(file=file.rds,srat)
```

### Resulting variable genes

The function `VariableFeatures` gives an impression of which 
genes are selected as being variable, i.e. informative. 
You can label the most highly variable ones using 
the `LabelPlot` function. 


```r
## select the most highly variable genes
## `VariableFeatures` returns them sorted by decreasing variance.)

top10 <- head(VariableFeatures(srat), 10)

plot1 <- VariableFeaturePlot(srat)

plot2 <- LabelPoints(plot = plot1, points = top10)

plot2
```

<!-- 
`RNA` assay: 
* `counts`: the raw counts
* `data` : log(1 + normalized(`counts`))
* `scale.date`: as `data` but scaled and centered

`SCT` assay: 
* `counts`: back-transformed data (see lecture)
* `data` : log(1 + `counts`)
* `scale.date`: Pearson residuals

<!-- lastly, save the complete sesssion for the next time -->

