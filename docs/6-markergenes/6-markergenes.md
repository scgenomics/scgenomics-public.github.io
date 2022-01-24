---
layout: default
title: Marker genes
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


# Plotting known marker genes.

Probably the most well-known marker gene for the development of the
intestine is Lgr5.  Let's plot it next to our UMAP. Let's do the same
for MUC2 which in our case is represented by the `muc2_mneongreen`
construct (the other one, `defa5_ires2_dsred`, has very low mRNA
expression although the protein product was found).


```r
p_type <- DimPlot(srat, pt.size=1, group.by='type', cols=type.colors)

## show LGR5 next to type plot
gene="LGR5"
p_lgr5 <- FeaturePlot(srat, pt.size=1, features=gene, 
    cols=CFG$gradient.colors, order=TRUE) + labs(title=gene)
p_type  | p_lgr5

gene="muc2-mneongreen"
p_muc2 <- FeaturePlot(srat, pt.size=1, features=gene, 
   cols=CFG$gradient.colors, order=TRUE) + labs(title=gene)
p_type  | p_muc2
```

You can clearly see expression of the MUC2 construct in 
the Goblet+Paneth cell cluster in the top-left. 

You can use `FeaturePlot` to show the expression of
several markers side by side. Canonical markers for Paneth cells are 
DEFA5, DEFA6, PLA2G2A, PRSS2, REG3A, and ITLN2. Let's see if
we can find them back.


```r
## Paneth:

paneth_markers <-  c('defa5-ires2-dsred',
  'DEFA6', 
  'PLA2G2A', 
  'PRSS2', 
  'REG3A', 
  'ITLN2')
# note that DEFA5 was removed from the genome so we use the chimeric gene

FeaturePlot(srat, pt.size=1, features=paneth_markers, order=TRUE,
  cols=CFG$gradient.colors)
```

The canonical markers are expressed in one 'corner' of the UMAP, but they
are also a bit confused.  We may be able to clean this by creating a
ModuleScore out of these markers, just like we did when looking at
cell cycle and stress.


```r
## create the modulescore:
srat <- AddModuleScore(srat, features=list(paneth=paneth_markers), name='Paneth')

## we can treat the modulescore as if it were a gene. Let's replace the
## defa5_ires2_dsred marker with it, and see if things improve now:

paneth_markers[1] <- 'Paneth1'          # replaces 'defa5-ires2-dsred'

FeaturePlot(srat, pt.size=1, cols=CFG$gradient.colors, order=TRUE,
  features=paneth_markers)
```

Well, maybe. It's at least a bit clearer than DEFA6, the canonical
Paneth cell marker.  Some of the markers for the other cell types are:

<!-- from ~/git/PanethAnalysis/R/plotting.R  (search for 'markers') -->

<!--  See e.g. Smith et al. J. Physiol. 2016
https://physoc-onlinelibrary-wiley-com.proxy.library.uu.nl/doi/full/10.1113/JP271651

more markers here: "Human Intestinal Organoids Maintain Self-Renewal Capacity 
and Cellular Diversity in Niche-Inspired Culture Condition"
Fujii et al. Cell Stem Cell, Volume 23

https://pubmed.ncbi.nlm.nih.gov/30526881/ 

-->

| Stem cells:                     SMOC2 
| Proliferating cells:            MKI67 
| Early enterocytes:              FABP1 
| Enterocytes:                    KRT20 
| Secretory progenitors:          HES6
| Goblet cells:                   MUC2
| Paneth cells:                   DEFA6  # let's replace that with our module score
| Enteroendocrine cells:          CHGA
| Tuft cells:                     AVIL

Let's see what their expression is.


```r
markers = c(
  SMOC2=               'SMOC2: Stem cells' ,
  MKI67=               'MKI67: Proliferating' ,
  FABP1=               'FABP1: Early enterocytes' ,
  KRT20=               'KRT20: Enterocytes' ,
  HES6=                'HES6: Secretory progen.' ,
  `muc2-mneongreen`=   'MUC2: Goblet' ,
  Paneth1=             'Paneth modulescore' ,
  CHGA=                'CHGA: Enteroendocrine' ,
  AVIL =               'AVIL: Tuft'
)

## create plot showing the expression of all our markers. Assign it to
# an object for later reference
markerplot <- FeaturePlot(srat, pt.size=1, order=TRUE,
   features=names(markers), 
   cols=CFG$gradient.colors)

## override the panel titles for readability
for(i in 1:length(markers))
    markerplot[[i]] <- ( markerplot[[i]] + labs(title = markers[i]) )

## show it:
markerplot
```

You can clearly see the UMAP 'making sense': each part represents a
different cell type. We can nearly assign  them 'by hand',
but that is infeasible. We need to cluster them first.

<!-- lastly, save the complete sesssion for the next time -->

