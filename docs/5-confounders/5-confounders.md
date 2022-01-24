---
layout: default
title: Confounding
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


# Confounders

In the previous steps, the variable genes were determined as part of
the normalization so as to avoid the noise and computational burden of
little-varying genes. 

We also may need to avoid the influence of factors that shape the data
without being informative. E.g., genes involved in the cell cycle could
lead to unrelated cells looking similar just because they are in
S-phase, leading to clustering by cell cycle phase rather than by
celltype. 

## Identifying confounding factors

The cell cycle is often a confounder, and the `LogNormalize`
splitting of the Goblet cells and the Secretory Progenitor cells shown
above is in fact due to that, as we will show now.

Seurat's `CellCycleScoring` method calculates an 'average S-phase score'
and also an 'average G2M-phase score' per cell. Based on these scores,
the most likely cell cycle phase is determined as being G1, G2M or S. The
phase-specific scores and the final phases are added to the `meta.data`
so they can be plotted.

The `CellCycleScoring` method needs sets of S-phase and G2M-phase
specific genes as input.  [Tirosh et al. Science
2016](https://dx.doi.org/10.1126/science.aad0501) established such sets,
and they are available from within the Seurat package. Doing
`data(cc.genes)` attaches a `list` with members `s.genes` and
`g2m.genes`. 


```r
data(cc.genes)

## add cell cycle phase estimates to both objects
srat_rna <- CellCycleScoring(object = srat_rna,
        s.features = cc.genes$s.genes,
        g2m.features = cc.genes$g2m.genes)

srat<- CellCycleScoring(object = srat,
        s.features =   cc.genes$s.genes,
        g2m.features = cc.genes$g2m.genes)

## show the UMAP with type and phase information for LogNormalized data:
title <- labs(title='LogNormalize')
p_type <-  DimPlot(srat_rna, pt.size=1, group.by='type', cols=type.colors)
p_phase <- DimPlot(srat_rna, pt.size=1, group.by='Phase') + title 

p_type | p_phase

## ... and same for SCTransformed data:
title <- labs(title='SCTransform')
p_type <- DimPlot(srat, pt.size=1, group.by='type', cols=type.colors )
p_phase <-  DimPlot(srat, pt.size=1, group.by='Phase') + title

p_type | p_phase
```

<!--  code to rename the missing genes see prior to 19-Jan-2022 21:01:20 -->

The `LogNormalize`'d data is more sensitive to confounding by cell-cycle phase, but
the SCTransformed data is not immune either: the S-phase specific genes of
*Goblet* and *Secretory Progenitor* are still clustering together.

<!-- maybe *not* remove cell cycle?  --->

For now, let's get rid of the LogNormalized 'rna' Seurat object and continue
with the other one, caling the garbarge collector (`gc`) to free up
space.


```r
rm(srat_rna)
gc()
```

### Module scores: 'average expression' of a group of genes

Let's look at stress as another source of confouding. For this we need
some kind of 'average expression of stress genes'. This can be
obtained using the so-called `ModuleScore` method. We indirectly
already used it when calling `CellCycleScoring` method.

We need a set of genes that is indicative of stress. The genes annotated
with the GO term "response to unfolded protein" will suffice, and are
available from the file `stress_genes.rds` . `AddModuleScore` adds a new
column to the `meta.data` which we can plot.

<!-- for literal genes see version prior to 13-Jan-2022 00:12:58  --->

```r
stress_genes <- readRDS(paste0(CFG$data_dir, '/stress_genes.rds'))

## show them:
stress_genes

srat <- AddModuleScore(srat, features=list(stress=stress_genes), name='stress')

## check the addition to the meta.data
head(srat@meta.data)

## as with the type colors, we will use a different color scheme that
## may work better. 
## An important other argument that we are passing is  order=TRUE.
## It makes sure that high values are printed 'in front' so they
## are not obscured

title <- labs(title='SCTransform')
p_type <-  DimPlot(srat, pt.size=1, group.by='type', cols=type.colors)
p_stress <- FeaturePlot(srat, pt.size=1, feature='stress1', order=TRUE,
    cols=CFG$gradient.colors) + title

p_type | p_stress
```

You can see a tight cluster in the early EC cells, so as with the cell-cycle effect, 
we may need to address this. 

## Dealing with confounding

There are a number of ways to deal with such confounding factors.  The
most drastic is to throw out all affected cells, which is probably too
drastic here. Seurat allows you to 'regress out' the effect, which means
subtracting the estimated confounder effect from all the expression
data. In our experience this does not work well, so we instead opt for a
simpler approach: ignoring genes that are known to be involved in the
confounding process.

Note that we do not _delete_ the confounders from the dataset. We
only exclude them from the list of variable genes. This way
these genes are still available for analyses and plotting. They are however
not taken into account anymore for the dimensionality reduction and
clustering.

After removing the offending genes from the list of variable features,
we have have to redo the PCA and the UMAP. We subsequently plot the UMAP
and judge by eye if there appears to be any influence of the confounders
on the structure of the data.


```r
remove <- unique( c(cc.genes$s.genes, cc.genes$g2m.genes, stress_genes) ) 

length(intersect(VariableFeatures(srat), remove)) # how many?

srat_clean <- srat # copy so we can compare

VariableFeatures(srat_clean) <- setdiff( VariableFeatures(srat_clean), remove)

## re-PCA and re-UMAP and compare

srat_clean <- RunPCA(srat_clean)

srat_clean <- RunUMAP(srat_clean, dims=1:CFG$ndims)

## Let's plot old and new, types and cell cycle phase in a four-way plot

old_type = DimPlot(srat, pt.size=1, group.by='type', cols=type.colors) + labs(title='old')
old_phase=  DimPlot(srat, pt.size=1, group.by='Phase')
new_type = DimPlot(srat_clean, pt.size=1, group.by='type', cols=type.colors) + labs(title='clean')
new_phase = DimPlot(srat_clean, pt.size=1, group.by='Phase')

## do the compositing. Note the use of parentheses to group things
( old_type | old_phase )  /  (new_type | new_phase )

## same but now with stress (note the reuse of previous plots)

old_stress = FeaturePlot(srat, pt.size=1, feature='stress1', cols=CFG$gradient.colors, order=TRUE)
new_stress = FeaturePlot(srat_clean, pt.size=1, feature='stress1', cols=CFG$gradient.colors, order=TRUE)

( old_type | old_stress )  /  (new_type | new_stress )
```

It looks like the cell cycle influence is somewhat resolved. In
particular, one of the big cycles (between TA and Goblet) is now broken,
which makes biological sense if the route from intestinal stem cells to
Goblet cells is always via secretory progenitors, never via
transit-amplifying cells.

## General 'deconfounding'

Confounding effects can be various, and depend on the sample at
hand. Some factors should probably be checked for by default:

* cellcyle
* stress
* metabolic activity, using ribosomal protein coding genes as a proxy
* sex

A general approach if you suspect a certain process may be confounding is to:

1. lookup the genes involved in that process using annotation databases such as Gene Ontology
1. create a ModuleScore of the genes involved
1. FeaturePlot the resulting ModuleScore
1. Check if there are clusters that are especially high for that score
1. If so, remove the genes from the VariableFeatures 
1. Redo the dimensional reduction
1. redo the FeaturePlot of the ModuleScore
1. If the cells with a high ModuleScore don't cluster anymore (or less), it worked
1. If not, undo the gene removal as these gene

In the above case we have removed some of the cell cycle influence, and while
stress is there, the observed clustering may not be caused only by stress.

There are more elaborate ways to get rid of confounding. A good method
is to not just look at the ModuleScore of the list of genes *known* to
be involved in a certain process, but also to check any genes
*correlating* with this ModuleScore. We will not go into this here.

Let's keep the 'deconfounded' object. Remember that the cleaning was done
on the `SCT` assay, but it also needs to be done on the `RNA`
(i.e. LogNormalized) assay which we need to use later on.

<!-- this is actually unneccesary -->


```r
### overwrite the old, un-deconfounded object with the clean one:
srat <- srat_clean 

## clean it up, it's big and unneeded
rm(srat_clean)
gc()

## IMPORTANT: do the same removal for the RNA assay
DefaultAssay(srat) <- 'RNA'
VariableFeatures(srat) <- setdiff( VariableFeatures(srat), remove)
DefaultAssay(srat) <- 'SCT'
```

<!-- lastly, save the complete sesssion for the next time -->

