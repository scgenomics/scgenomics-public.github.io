
### autogenerated from R -e 'knitr::purl(14_enrich.Rmd)'
## ----setup, include=FALSE, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## suppressMessages(library(rmarkdown))
## knitr::opts_chunk$set(echo = TRUE,
##   fig.width=9,
##   fig.height=6,
##   cache=FALSE, # only use TRUE if you're really sure. Also check xfun::cache_rds
##   collapse=TRUE,  # put misc. outputs from one chunk in one block
##   tidy=TRUE, # tidies up the R code
##   tidy.opts=list(arrow=TRUE,
##                  indent=2,
##                  width.cutoff = 60))


## ----loadsession, include=FALSE, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## source('../libs.R')
## load("../13_clust/session.rda")


## ----de-genes, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## 
## ## first step: find the differentially expressed genes
## ## this takes ~ 2min
## de_genes <- FindAllMarkers(srat, assay='RNA',
##                            only.pos = TRUE,
##                            min.pct = 0.1,
##                            logfc.threshold = log(1.5))
## 
## ## save them:
## file <- paste0(CFG$output_dir, "/rms_de_genes.txt")
## write.table(de_genes, file=file, sep = "\t",
##     quote = FALSE, row.names = T, col.names = NA, na = "")
## 
## ## The resulting data.frame is too unwieldy to work with,
## ## we the top 10 genes per cluster are enough:
## de_genes %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
## ## (after this incantatation, top10 contains the slimmed-down list)
## 
## ## show the contents, maybe you already notice interesting genes:
## head(top10)
## 
## ## show it visually:
## DoHeatmap(subset(srat, downsample = 50),
##           features = top10$gene,
##           group.colors=cluster.colors,
##           assay='RNA',
##           slot='scale.data'
##           ) + theme(axis.text.y=element_text(size=6))
## 


## ----enrichment, eval=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## 
## ## compareCluster crashes on gene symbols, so we have tot translate
## ## to EnsEMBL ids:
## anno <- readRDS(file=paste0(CFG$data_dir, "/annotation.rds"))
## lookup_ensg_id = function(id){anno[id,'gene_id'] }
## 
## ## compareClusters needs a list of clusters, not a
## ## data.frame. Convert our top10:
## df  = as.data.frame(top10)
## clusters = split(df$gene, df$cluster)
## 
## ## show contents:
## show(clusters)
## 
## ## make EnsEMBL id's out of them:
## clusters = lapply(clusters, lookup_ensg_id)
## 
## ## show it:
## show(clusters)
## 
## ## lastly run it. This takes ~ 2 min
## enrichment <- compareCluster(clusters,
##       fun= "enrichGO",
##       OrgDb = "org.Hs.eg.db",
##       keyType = "ENSEMBL",
##       ont = "BP", minGSSize = 3,
##       pAdjustMethod = "BH",
##       pvalueCutoff = 0.01,
##       qvalueCutoff = 0.05,
##       readable = TRUE)
## 
## ## so save it:
## file=paste0(CFG$output_dir, '/enrichment.rds')
## saveRDS(file=file, enrichment)
## 
## ## show results:
## clusterProfiler::dotplot(enrichment,
##    showCategory = 4,
##    title = "GO biological process",
##    font = 5)
## ## output can be simplified with simplify and dropGO, see the clusterProfiler vignette.
## 


## ----savesession, include=FALSE, eval=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## ## save complete state for next lesson
## save(file="session.rda", list=c("srat", ls(pat="col|CFG")))

