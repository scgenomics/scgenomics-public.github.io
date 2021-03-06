---
layout: default
title: Calling tumor cells with inferCNV
---

<!-- stuff to make Rmarkdown do what we want:  -->


<!-- load complete state from previous lesson -->


# InferCNV

<!--  code for by patient assignment: see prior to 21-Jan-2022 00:17:11 -->

We have suspicions regarding which clusters are likely to be tumor cells, but
we're not certian. In particular the mesenchymal cluster is unclear, they might 
well be healthy cells.

Some tumor cells harbour chromosomal amplifications or deletions that
can be detected by looking for stretches of adjacent genes that are all,
on average, being up or down relative to a
reference. [infercnv](https://github.com/broadinstitute/inferCNV) is a
tool to do this.  It needs healthy cells as a reference. We can use the
cells with a clear, definite cell type that we just found using SingleR for
this.


```r
## overview of numbers per celltype:
table(srat@meta.data$type2)

## which cells we think we can trust as being healthy:
ref_types <- c("T", "NK", "B", "Mono")

## overview of this:
maybe_tumor <- !(srat@meta.data$type2 %in% ref_types)
table( maybe_tumor )

## infercnv needs them as an external data file in a specific format,
## create that:
df = data.frame(cell=colnames(srat), type=srat@meta.data$type2)


celltypes_file =paste0(CFG$output_dir,"/celltypes.txt")
write.table(df, file=celltypes_file, sep="\t", 
  quote=FALSE, na="", row.names=FALSE, col.names=FALSE)

## infercnv also needs to gene coordinates, can find them
## in the data directory:
geneorder_file = paste0(CFG$data_dir, "/geneorder.txt")

## set up the infercnv object with all necessary information:
infercnv_obj = infercnv::CreateInfercnvObject(delim="\t",
      raw_counts_matrix=GetAssayData(srat, assay='RNA', slot='counts'),
      annotations_file=celltypes_file,
      gene_order_file=geneorder_file,
      ref_group_name = ref_types )

## free up space before we start:
gc()

## infercnv writes its output to disk. It is a fair amount,
## so let's create a separate directory for that and write
## it there:
outdir=paste0(CFG$output_dir, "/infercnv")
dir.create(outdir)

## running takes can take up to 10 min, best start this before the lecture:
infercnv_obj <- infercnv::run(infercnv_obj,
                             cutoff=0.3, 
                             out_dir=outdir,
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             analysis_mode = "samples",
                             num_threads = 3,
                             output_format = "png",
                             window_length = 101,
                             save_final_rds = TRUE,
                             plot_steps=FALSE, 
                             sd_amplifier = 2,
                             HMM=FALSE)
## plot_steps=TRUE may lead to weird crash, see github
## cutoff: min average read counts per gene for reference cells
## sd_amplifier: noise threshold. Higher: less noise but also less signal.

## This took a while so let's save this object:
file <- paste0(CFG$output_dir, "/infercnv_obj.rds")
saveRDS(file=file,infercnv_obj)
```


<!-- {r include-cnv-plot, echo = TRUE, out.width="30%" include=FALSE } -->

<!-- following is needed to get the plot into the knitted product as it is external -->



## Results

The script does not produce output to the *Plots* tab of RStudio, but instead
writes directly to disk, in directory `outdir` that we just created. 
You can access this folder in the *Files* tab. There you can see, 
amongst others, the file `infercnv.png`, which contains the final
inferCNV plot.  You should be able to open it by clicking.

It should look like

![infercnv](figure/infercnv.png)

There is no easy way to extract a CNV-score from the `infercnv_obj`, so currently
we have to do the assignment by eye. It seems safe to say that all the 
non-reference cells are indeed tumor. 

<!-- @@@ note: discuss the HLA cluster ( chr6, p21.3, halfway p-arm) 
@@@ mesen has no CNV's, nice.

ESC-2+3 are patient 410;  close together

ESC-6+7 are patient GRN (note the split, prolly artefact)

ESC-0 is patient 'CPU'
-->


What can you say about the mesenchymal cluster?


<!-- lastly, save the complete sesssion for the next time -->

