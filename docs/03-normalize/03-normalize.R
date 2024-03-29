### autogenerated from R -e 'knitr::purl(03-normalize.Rmd)'
### autogenerated from R -e 'knitr::purl(03-normalize.Rmd)'
### autogenerated from R -e 'knitr::purl(03-normalize.Rmd)'
### autogenerated from R -e 'knitr::purl(03-normalize.Rmd)'
## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
suppressMessages(library(rmarkdown))
knitr::opts_chunk$set(echo = TRUE,
  fig.width=9,
  fig.height=6,
  cache=FALSE, # only use TRUE if you're really sure. Also check xfun::cache_rds
  collapse=TRUE,  # put misc. outputs from one chunk in one block
  tidy=TRUE, # tidies up the R code
  tidy.opts=list(arrow=TRUE, 
                 indent=2,
                 width.cutoff = 60))


## ----loadsession, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source('../libs.R')
load("../02-filtercells/session.rda")
gc()



## ----load-celltypes-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## load the cell types from disk:

file <- paste0(CFG$data_dir, '/celltypes.rds')
celltypes <- readRDS(file)
srat <- AddMetaData(srat, col.name='type', metadata=celltypes)

types <- sort(unique(srat@meta.data$type))
type.colors <- Polychrome::dark.colors( length(types ))
names(type.colors) <- types
## use as DimPlot( ... , cols=type.colors, )



## ----do-log-norm----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

srat <- NormalizeData(srat,normalization.method = "LogNormalize")

srat <- ScaleData(srat, features = rownames(srat), verbose=FALSE)

srat <- FindVariableFeatures(srat)



## ----do-sct, output.lines=(1:20)------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## we wrap the call to SCTransform in suppressMessages to avoid 
## the verbosity. Feel free to leave out.
srat <- suppressMessages(SCTransform(srat, 
                         vars.to.regres=NULL,
                         vst.flavor='v2', ## version that solved some problems
                         verbose=FALSE,
                         variable.features.n=3000))

## note that the assay has now changed:
DefaultAssay(srat)

## Since SCTransform normalization takes a bit longer, it may be useful 
## to store the resulting object. That way you can easily go back if things
## go wrong.

## file=paste0(CFG$output_dir, "/srat-normalized.rds")
## saveRDS(file=file,srat)
## To read the object back in, do
## srat <- readRDS(file=file,srat)



## ----show-var-genes-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## select the most highly variable genes
## `VariableFeatures` returns them sorted by decreasing variance.)

top10 <- head(VariableFeatures(srat), 10)

p.varfeat <- VariableFeaturePlot(srat)

## show plot without labels:
p.varfeat

## add labels:
p.varfeat <- LabelPoints(plot = p.varfeat, points = top10, repel=TRUE)

p.varfeat



## ----savesession, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## save complete state for next lesson
## save.image("session.rda")
if(Sys.getenv("nodump")!="TRUE")
    save(file="session.rda", list=c("srat", ls(pat="col|CFG")))


