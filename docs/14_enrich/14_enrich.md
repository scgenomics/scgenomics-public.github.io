---
layout: default
title: Enrichment analysis
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


# Enrichment analysis

We will now do an enrichment analysis to find out what might be characteristic for our 
clusters. This hopefully will give insight into the role or activity of the cells in the
clusters. This analysis will look for over-representation  of 'annotation terms' in the
genes that significantly up in a cluster, relative to the rest of the cells.

First we have to find those characterisic genes per cluster. This is
done based on Differential Expression analysis. In the next step these
genes are analyzed.  There are many ways of doing this, and many
catalogues of 'annotation terms' to do it.  For simplicity, we will only
use the hypergeometric test and [Gene
Ontology](http://geneontology.org/) *Biological Process* terms here.

## Find differentially expressed genes


```r
## first step: find the differentially expressed genes
## this takes ~ 2min
de_genes <- FindAllMarkers(srat, assay='RNA', 
                           only.pos = TRUE,
                           min.pct = 0.1,
                           logfc.threshold = log(1.5))

## save them:
file <- paste0(CFG$output_dir, "/rms_de_genes.txt")
write.table(de_genes, file=file, sep = "\t", 
    quote = FALSE, row.names = T, col.names = NA, na = "")

## The resulting data.frame is too unwieldy to work with, 
## we the top 10 genes per cluster are enough:
de_genes %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
## (after this incantatation, top10 contains the slimmed-down list)

## show the contents, maybe you already notice interesting genes:
head(top10)

## show it visually:
DoHeatmap(subset(srat, downsample = 50),
          features = top10$gene,
          group.colors=cluster.colors,
          assay='RNA',
          slot='scale.data'
          ) + theme(axis.text.y=element_text(size=6))
```

## Perform enrichment analysis 


```r
## compareCluster crashes on gene symbols, so we have tot translate
## to EnsEMBL ids:
anno <- readRDS(file=paste0(CFG$data_dir, "/annotation.rds"))
lookup_ensg_id = function(id){anno[id,'gene_id'] } 

## compareClusters needs a list of clusters, not a 
## data.frame. Convert our top10:
df  = as.data.frame(top10)
clusters = split(df$gene, df$cluster)

## show contents:
show(clusters)

## make EnsEMBL id's out of them:
clusters = lapply(clusters, lookup_ensg_id)

## show it:
show(clusters)

## lastly run it. This takes ~ 2 min
enrichment <- compareCluster(clusters,
      fun= "enrichGO",
      OrgDb = "org.Hs.eg.db",
      keyType = "ENSEMBL", 
      ont = "BP", minGSSize = 3,
      pAdjustMethod = "BH", 
      pvalueCutoff = 0.01, 
      qvalueCutoff = 0.05,
      readable = TRUE)

## so save it:
file=paste0(CFG$output_dir, '/enrichment.rds')
saveRDS(file=file, enrichment)

## show results:
clusterProfiler::dotplot(enrichment,
   showCategory = 4, 
   title = "GO biological process",
   font = 5)
## output can be simplified with simplify and dropGO, see the clusterProfiler vignette.
```

<!-- question: based on this, what would be cell types?  -->


<!-- lastly, save the complete sesssion for the next time -->

