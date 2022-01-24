--- 
layout: default
title: Hemoglobin removal
--- 


<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


# Hemoglobin removal 

When plotting the number of genes versus the number of transcripts per plate,
something unusual is going on. We have to check this.


```r
plates <- unique(srat@meta.data$plate_id)

## for convenience we will plot everything in a loop:
for (plate in plates) { 

df <- subset( srat@meta.data, plate_id==plate)

p =  ggplot(df,
       aes(x=nCount_RNA, y=nFeature_RNA, color=plate_id)) + geom_point() + 
       labs(title=plate)

## for some reason, ggplot2 plot objects won't plot automaticlly inside a loop.
## We therefore have explicitly call it with show():
show(p)

## and/or logarithmic:

p + scale_x_continuous(trans='log2') + 
   scale_y_continuous(trans='log2')

, eval=FALSE}
```

These 'off-diagonal' cells express the same number of transcripts with
fewer genes, i.e. they are more specialized. A typical example of such
cells are erythrocytes containing predominantly hemoglobin mRNA's.
The cells we see here maybe be proper cells that are contaminated 
with ambient hemoglobin mRNA from erythrocytes.

Let's show this and filter them out, using `hemo_genes.rds`, a
list of heme genes provided in your data directory.



```r
## read the file in:
file <- paste0(CFG$data_dir, "/hemo_genes.rds")
hb_genes <- readRDS(file)
hb_genes <- intersect(hb_genes, rownames(srat)) # in case not all are expressed

hb_counts <- Matrix::colSums(srat@assays$RNA@counts[ hb_genes, ])

srat <- AddMetaData(srat,
                    col.name = "log2hb_genes",
                    metadata=log2(1+hb_counts))

srat <- AddMetaData(srat,
                    col.name = "pct_hemo",
                    metadata=100*hb_counts/srat@meta.data$nCount_RNA )

## now show the previous plot along side a plot of the hemo content:
for (plate in plates) {

  df <- srat@meta.data[  srat@meta.data$plate_id == plate, ]
  
  p_ngenes <- ggplot(df, aes(x=log2_counts, y=log2_features, color= pct_hemo)) + 
    geom_point(size=1, alpha=0.5) +
    scale_color_gradient(low="blue", high="red", limits=c(0,5), oob=scales::squish) +
    labs(title=plate)

  p_pcthemo <- ggplot(df, aes(x=log2_counts, y=pct_hemo, color= log2_features)) + 
     geom_point(size=1, alpha=0.5) +
     scale_color_gradient(low="blue", high="red", limits=c(2,16), oob=scales::squish)

## again using show() because we're inside a loop
show(p_ngenes | p_pcthemo)

, eval=FALSE}

## set a maximum to get rid of cells, and do the selection
CFG$max_pct_hemo <- 5                   #

dim(srat)
dim(subset(srat, pct_hemo > CFG$max_pct_hemo))

srat <- subset(srat, pct_hemo <= CFG$max_pct_hemo)

## maybe save
file <- paste0(CFG$output_dir, "/srat-cellsfiltered.rds")
saveRDS(file=file, srat)
```

Post scriptum: a way to find explore peculiarties of groups of cells is
to use the `CellSelector ` to pick cells from a plot. It works roughly
as follows:


```r
plot <- DimPlot(srat, reduction = "umap") 

selected <- CellSelector(plot = plot)  # this launches the cell picker

## afterwards,  selected contains the names of the picked cells.
## If your override their Identity you can use FindMarkers  to find differentially
## expressed genes:
Idents(srat, cells = selected ) <- "maybe_erythro"

markers <- FindMarkers(srat, ident.1 = "maybe_erythro", 
   only.pos = TRUE, etc.) 

head(markers) 
```


<!-- lastly, save the complete sesssion for the next time -->

