setGeneric("smutans.title", function(object) standardGeneric("smutans.title"))
setMethod("smutans.title",
          "Smutans",
  function( object ) {
    object@title
  }
)

setGeneric("smutans.title<-", function(object,value) standardGeneric("smutans.title<-"))
setReplaceMethod("smutans.title", 
                 "Smutans",
  function( object, value ) {
    object@title <- value
    object
  }
)

setGeneric("smutans.plotDispersionEstimates", 
  function(object, file="default", 
           xlab="Base mean expression", ylab="Dispersion values", qval=0.05,...) 
           standardGeneric("smutans.plotDispersionEstimates")
)
setMethod("smutans.plotDispersionEstimates", 
          "Smutans", 
  function(object, file="default", 
           xlab="Base mean expression", ylab="Dispersion values", qval=0.05 ) 
  {
    if (file != "default")
    {
      postscript(file, width=6, height=6, horizontal = FALSE, onefile = FALSE, paper = "special")
      oldpar <- par (mar=c(4.5, 4.5, 0.5, 0.5))
    }
    plot( rowMeans( counts( object@cds, normalized=TRUE ) ),
      fitInfo(object@cds)$perGeneDispEsts,
      xlab=xlab, ylab=ylab,
      pch = '.', log="xy" )
    xg <- 10^seq( -.5, 5, length.out=300 )
    lines( xg, fitInfo(object@cds)$dispFun( xg ), col="red" )
    if (file != "default")
    {
      par(oldpar)
      dev.off()
    }
  }
)

setGeneric("smutans.plotDiffExp", 
  function( object, file="default", 
            xlab="Base mean expression", ylab="Log base 2 fold changes", qval=0.05 ) 
    standardGeneric("smutans.plotDiffExp")
)
setMethod("smutans.plotDiffExp", 
          "Smutans",
  function(object, file="default", 
           xlab="Base mean expression", ylab="Log base 2 fold changes", qval=0.05 ) 
  {
    if (file != "default")
    {
      postscript(file, width=6, height=6, horizontal = FALSE, onefile = FALSE, paper = "special")
      oldpar <- par (mar=c(4.5, 4.5, 0.5, 0.5))
    }
    plot(
      object@res$baseMean,
      object@res$log2FoldChange,
      xlab=xlab, ylab=ylab,
      log="x", pch=20, cex=.3,
      col = ifelse( object@res$padj < qval, "red", "black" ) )
    if (file != "default")
    {
      par(oldpar)
      dev.off()
    }
  }
)

setGeneric("smutans.de2type", function(object,...) standardGeneric("smutans.de2type"))
setMethod("smutans.de2type", 
          "Smutans", 
  function(object, type="", condition="", condA, condB) {
    ua159Samples <- pData(object@countDataSet)$type == type
    countsTable <- counts(object@countDataSet)[ , ua159Samples ]
    conds <- pData(object@countDataSet)$condition[ ua159Samples ]
    cds <- newCountDataSet( countsTable, conds )
    cds <- estimateSizeFactors( cds )
    object@cds <- estimateDispersions( cds )
    object@res <- nbinomTest( object@cds, condA, condB )
    object@pval <- object@res$pval
    object@padj <- object@res$padj
  
    object
  }
)

setGeneric("smutans.de2", 
  function(object, compareCondition="yes", 
           type="", condition="", condA="", condB="",
           method="pooled") standardGeneric("smutans.de2")
)
setMethod("smutans.de2", 
          "Smutans", 
  function(object, compareCondition="yes", 
           type="", condition="", condA="", condB="",
           method="pooled" ) {
    if (type=="" && condition=="") {
      design <- pData( object@countDataSet )[,c("type","condition")]
      fullCountsTable <- counts( object@countDataSet )
      cdsFull <- newCountDataSet( fullCountsTable, design )
      cdsFull <- estimateSizeFactors( cdsFull )
      object@cds <- estimateDispersions( cdsFull, method=method )
      if (compareCondition == "yes") {
        fit1 <- fitNbinomGLMs( object@cds, count ~ type + condition )
        fit0 <- fitNbinomGLMs( object@cds, count ~ type  )
      } else {
        fit1 <- fitNbinomGLMs( object@cds, count ~ type + condition )
        fit0 <- fitNbinomGLMs( object@cds, count ~ condition )
      }
      object@pval <- nbinomGLMTest( fit1, fit0 )
      object@padj <- p.adjust( object@pval, method="BH" )
      object@res <- data.frame(fit1,pval=object@pval,padj=object@padj)
    } else if (type!="" && condition=="") {
      ua159Samples <- pData(object@countDataSet)$type == type
      countsTable <- counts(object@countDataSet)[ , ua159Samples ]
      conds <- pData(object@countDataSet)$condition[ ua159Samples ]

      cds <- newCountDataSet( countsTable, conds )
      cds <- estimateSizeFactors( cds )
      object@cds <- estimateDispersions( cds )
      object@res <- nbinomTest( object@cds, condA, condB )
      object@pval <- object@res$pval
      object@padj <- object@res$padj
    } else if (type=="" && condition!="") {
      ua159Samples <- pData(object@countDataSet)$condition == condition
      countsTable <- counts(object@countDataSet)[ , ua159Samples ]
      conds <- pData(object@countDataSet)$type[ ua159Samples ]

      cds <- newCountDataSet( countsTable, conds )
      cds <- estimateSizeFactors( cds )
      object@cds <- estimateDispersions( cds )
      object@res <- nbinomTest( object@cds, condA, condB )
      object@pval <- object@res$pval
      object@padj <- object@res$padj
    }
    object
  }
)

setGeneric("smutans.padj", function(object) standardGeneric("smutans.padj"))
setMethod("smutans.padj",
          "Smutans",
  function( object ) {
    object@padj
  }
)

setGeneric("smutans.pval", function(object) standardGeneric("smutans.pval"))
setMethod("smutans.pval",
          "Smutans",
  function( object ) {
    object@pval
  }
)

setGeneric("smutans.de2WithoutType", 
  function(object,...) standardGeneric("smutans.de2WithoutType")
)
setMethod("smutans.de2WithoutType", 
          "Smutans", 
  function(object, condA, condB) {
    design <- pData( object@countDataSet )[,c("type","condition")]
    fullCountsTable <- counts( object@countDataSet )
    cdsFullB <- newCountDataSet( fullCountsTable, design$condition )
    cdsFullB <- estimateSizeFactors( cdsFullB )
    object@cds <- estimateDispersions( cdsFullB )
    object@res <- nbinomTest( object@cds, condA, condB )
    object@pval <- object@res$pval
    object@padj <- object@res$padj
    object
  }
)

setGeneric("smutans.de2TypeLfc", 
  function(object,...) standardGeneric("smutans.de2TypeLfc")
)
setMethod("smutans.de2TypeLfc", 
          "Smutans", 
  function(object, condA, condB, file="default") {
    cdsBlind <- estimateDispersions( object@cds, method="blind" )
    vsd <- getVarianceStabilizedData( cdsBlind )
    mod_lfc <- (rowMeans( vsd[, conditions(object@cds)==condB, drop=FALSE] ) -
                rowMeans( vsd[, conditions(object@cds)==condA, drop=FALSE] ))

    lfc <- object@res$log2FoldChange
    finite <- is.finite(lfc)
    # table(as.character(lfc[!finite]), useNA="always")

    largeNumber <- 10
    lfc <- ifelse(finite, lfc, sign(lfc) * largeNumber)

    if (file != "default")
    {
      postscript(file, width=6, height=6, horizontal = FALSE, onefile = FALSE, paper = "special")
      oldpar <- par (mar=c(4.5, 4.5, 0.5, 0.5))
    }
    # plot( lfc, mod_lfc, pch=20, cex=.3, col = ifelse( finite, "#80808040", "red" ) )
    plot( lfc, mod_lfc, 
          xlab="Log fold changes", ylab="Moderated log fold changes", 
          pch=20, cex=.3, col = ifelse( finite, "black", "red" ) )
    abline( a=0, b=1, col="#40404040" )
    if (file != "default")
    {
      par(oldpar)
      dev.off()
    }
  }
)

setGeneric("smutans.de2TypeHeatmap", 
  function(object,file="default", nselect=40) standardGeneric("smutans.de2TypeHeatmap")
)
setMethod("smutans.de2TypeHeatmap", 
          "Smutans", 
  function(object,file="default", nselect=40) {
    cdsBlind <- estimateDispersions( object@cds, method="blind" )
    vsd <- getVarianceStabilizedData( cdsBlind )
    select <- order(object@res$pval)[1:nselect]
    colors <- colorRampPalette(c("white","black"))(100)
    # colors <- colorRampPalette(brewer.pal(9,"Blues"))(100)

    if (file != "default")
    {
      postscript(file, width=6, height=6, horizontal = FALSE, onefile = FALSE, paper = "special")
      oldpar <- par (mar=c(4.5, 4.5, 0.5, 0.5))
    }
    # heatmap( vsd[select,], col = colors, scale = "none")
    heatmap.2( vsd[select,], col = colors, scale = "none", trace="none",
               density.info="none" ) 
    if (file != "default")
    {
      par(oldpar)
      dev.off()
    }
  }
)

setGeneric("smutans.de2Clust", 
  function(object,...) standardGeneric("smutans.de2Clust")
)
setMethod("smutans.de2Clust", 
          "Smutans", 
  function(object, file="default", margins=c(10,10)) {
    if (file != "default")
    {
      postscript(file, width=6, height=6, horizontal = FALSE, onefile = FALSE, paper = "special")
      oldpar <- par (mar=c(4.5, 4.5, 0.5, 0.5))
    }
    design <- pData( object@countDataSet )[,c("type","condition")]
    fullCountsTable <- counts( object@countDataSet )
    cdsFull <- newCountDataSet( fullCountsTable, design )
    cdsFull <- estimateSizeFactors( cdsFull )
    cdsFullBlind <- estimateDispersions( cdsFull, method = "blind" )
    vsdFull <- getVarianceStabilizedData( cdsFullBlind )
    dists <- dist( t( vsdFull ) )
    heatmap( as.matrix( dists ),
       symm=TRUE, scale="none", margins=margins,
       col = colorRampPalette(c("black","white"))(100),
       labCol = paste( pData(cdsFullBlind)$condition, pData(cdsFullBlind)$type),
       labRow = paste( pData(cdsFullBlind)$condition, pData(cdsFullBlind)$type)
       )
#    heatmap( as.matrix( dists ), Colv=TRUE,
#       symm=TRUE, scale="none", margins=margins,
#       col = colorRampPalette(c("black","white"))(100),
#       trace="none", density.info="none",
#       labRow = paste( pData(cdsFullBlind)$condition, pData(cdsFullBlind)$type),
#       labCol = paste( pData(cdsFullBlind)$condition, pData(cdsFullBlind)$type),
#       dendrogram = "both", keysize = 1.2
#       )
    if (file != "default")
    {
      par(oldpar)
      dev.off()
    }
  }
)

setGeneric("smutans.de2List", 
  function(object,file="default") standardGeneric("smutans.de2List")
)
setMethod("smutans.de2List", 
          "Smutans", 
  function(object, file="default") {
    if (file == "default") {
      return( object@res )
    } else {
      write.csv(object@res, file=file)
    }
  }
)

setGeneric("smutans.de2Genes", 
  function(object,qval=0.05) standardGeneric("smutans.de2Genes")
)
setMethod("smutans.de2Genes", 
          "Smutans", 
  function(object, qval=0.05) {
    # resSig <- object@res[object@padj < qval,]
    # resSig$id
    rownames(counts(object@cds))[object@padj < qval]
  }
)

setGeneric("smutans.de2Goseq", 
  function(object,qval=0.05, feature.genes, go.genes, cat.desc) 
    standardGeneric("smutans.de2Goseq")
)
setMethod("smutans.de2Goseq", 
          "Smutans", 
  function(object, qval=0.05, feature.genes, go.genes, cat.desc) {
    length.genes = feature.genes$V3 - feature.genes$V2
    assayed.genes = feature.genes$V4

    de.genes <- smutans.de2Genes ( object, qval=qval )
    gene.vector = as.integer(assayed.genes %in% de.genes)
    names(gene.vector) = assayed.genes

    rm(feature.genes)
    pwf = nullp( gene.vector, bias.data=length.genes, plot.fit=FALSE )

    # go.hypergeometric = goseq(pwf,gene2cat=go.genes,method="Hypergeometric") # No length bias correction
    go = goseq(pwf,gene2cat=go.genes,method="Wallenius") # Length bias correction - Approximation
    # go.sample = goseq(pwf,gene2cat=go.genes,method="Sampling",repcnt=10000) # Length bias correction - Sampling

    #go.fdr = go[p.adjust(go$over_represented_pvalue,method="BH")<.05,]
    go.fdr <- data.frame(go, qval=p.adjust(go$over_represented_pvalue,method="BH"))
    go.fdr <- go.fdr[go.fdr$qval < .05,]

    for (i in go.fdr$category)
    {
      cat(go.fdr$qval[go.fdr$category==i], " & ", as.character(cat.desc$V3[cat.desc$V1==i]))
      cat("\\\\\n")
    }
  }
)



