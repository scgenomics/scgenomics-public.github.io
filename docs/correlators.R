### functions for finding additional confounders based on corrrelations.
### Written by Philip Lijnzaad <plijnzaad@gmail.com> based on idea by
### Thanasis Margaritis. (Original is in misc-functions.R dated 20-Jan-2022)

metadataCorrelations <- function(srat,var,
                             assay.type='SCT',
                             method='pearson') {
  restore.type <- DefaultAssay(srat)
  
  v <- srat@meta.data[[var]]
  if(is.null(v) || !is.numeric(v))
    stop("Variable ", var, " not found in cell meta data or not numeric")
  var <- v
  data <-srat@assays[[assay.type]]@data  
  
  cor <- cor(var, as.matrix(t(data)), method=method)
  stopifnot(nrow(cor)==1)
  dim(cor) <- NULL 
  names(cor) <- rownames(srat)
  cor[is.na(cor)] <- 0
  
  restore.type -> DefaultAssay(srat)
  return(cor)
} # metadataCorrelations

derivedCellcycleGenes <- function(srat,
                                  assay.type='SCT',
                                  quantile=0.25,
                                  mindiff=0.05,
                                  plot=TRUE) {
  
  data(cc.genes) # built into Seurat
  
  if ( sum(c('S.Score','G2M.Score') %in% names(srat@meta.data)) !=2 )
    stop("Has CellCycleScoring been run on this object??")

  Scor <- metadataCorrelations(srat,
                               var='S.Score',
                               assay.type=assay.type,
                               method='pearson')
  G2Mcor <- metadataCorrelations(srat,
                                 var='G2M.Score',
                                 assay.type=assay.type,
                                 method='pearson')
  
  stopifnot(is.character(names(Scor)) && is.character(names(G2Mcor)))
  stopifnot(identical(names(Scor), names(G2Mcor)))
    
  Sgenes <- intersect(cc.genes$s.genes, names(Scor))
  G2Mgenes <- intersect(cc.genes$g2m, names(Scor))
    
  stopifnot(is.character(Sgenes) && is.character(G2Mgenes))

  cors <- list(Sall = Scor,                       # all 
               S_Sgenes = Scor[Sgenes],             # just the Sgenes
               S_G2Mgenes = Scor[G2Mgenes],
               
               G2Mall = G2Mcor,
               G2M_G2Mgenes = G2Mcor[G2Mgenes],
               G2M_S = G2Mcor[Sgenes]
               )

  ## find cutoffs based on quantiles. Make sure they don't overlap (too much)
  par(pty='m', mfrow=c(1,3))
  Scors <- cors[ grep("^S",names(cors)) ]
  names(Scors) <- c('all genes', 'S_genes', 'G2M_genes')
  S.cutoff <- quantile(Scors$S_genes, quantile)
  
  G2Mcors <- cors[ grep("^G",names(cors)) ]
  names(G2Mcors) <- c('all genes', 'G2M_genes', 'S_genes')
  G2M.cutoff <- quantile(G2Mcors$G2M_genes, quantile)

  if(plot) { 
    boxplot(Scors, ylim=c(0,1), main='correlations to S score', ylab='correlation')
    abline(h=S.cutoff, col = "red", lty = 2)
    boxplot(G2Mcors, ylim=c(0,1), main='correlations to G2M score')
    abline(h=G2M.cutoff, col = "blue", lty = 2)
  }

  CC.cors <- data.frame(cors[c('Sall', 'G2Mall')])
  names(CC.cors) <- c('S', 'G2M')
  
  S.derived <- rownames(subset(CC.cors,
                                   S > S.cutoff
                                   & S - G2M >= mindiff))
  S.derived <- setdiff(S.derived, Sgenes)

  G2M.derived <- rownames(subset(CC.cors,
                                       G2M > G2M.cutoff
                                       & G2M - S >= mindiff))

  G2M.derived <- setdiff(G2M.derived, G2Mgenes)

  if(plot) {
#    par(mfrow=c(1,1), pty='s')
    xylim <- c(-0.2, 0.8)

    cex <- 0.6
    with(CC.cors,
         plot(S ~ G2M, pch = 19, cex = cex-0.1, xlim=xylim, ylim=xylim,
              xlab='G2M cor', ylab='S cor'))
    abline(h=0, v=0, col='grey')
    abline(h = S.cutoff, lty = 2, col = "red", lwd=2)
    abline(v = G2M.cutoff, lty = 2, col = "blue", lwd=2)
    abline(a= mindiff, b=1, col='grey', lwd=2, lty=2)
    abline(a= -mindiff, b=1, col='grey', lwd=2, lty=2)
    ## overplot the genes we found
    with(CC.cors[S.derived, ], 
         points(S ~ G2M, col = "red", pch = 1, cex=cex*2))
    with(CC.cors[G2M.derived, ], 
         points(S ~ G2M, col = "blue", pch = 1, cex=cex*2))
    with(CC.cors[Sgenes, ], 
         points(points(S ~ G2M,col = "red", pch = 19, cex=cex*1.2)))
    with(CC.cors[G2Mgenes, ], 
         points(S ~ G2M,col = "blue", pch = 19, cex=cex*1.2))
    legend(x='topleft', lty=NULL,
           legend=c('S genes', 'additional S',
                    'G2M genes', 'additional G2M'),
           col=c('red', 'red', 'blue', 'blue'),
           pch=c(19, 1, 19, 1))
  }

  return(list(S.derived=S.derived, G2M.derived=G2M.derived))
} # derivedCellcycleGenes
