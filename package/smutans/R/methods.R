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

setGeneric("smutans.counts", function(object) standardGeneric("smutans.counts"))
setMethod("smutans.counts",
          "Smutans",
  function( object ) {
    counts( object@countDataSet )
  }
)

setGeneric("smutans.counts<-", function(object,value) standardGeneric("smutans.counts<-"))
setReplaceMethod("smutans.counts", 
                 "Smutans",
  function( object, value ) {
    object@counts <- value
    object
  }
)

setGeneric("smutans.plotDispersionEstimates", function(object,...) standardGeneric("smutans.plotDispersionEstimates"))
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

setGeneric("smutans.plotDiffExp", function(object,...) standardGeneric("smutans.plotDiffExp"))
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

setGeneric("smutans.de2", function(object,...) standardGeneric("smutans.de2"))
setMethod("smutans.de2", 
          "Smutans", 
  function(object, compareCondition="yes", type="", condition="", condA="", condB="") {
    if (type=="" && condition=="") {
      design <- pData( object@countDataSet )[,c("type","condition")]
      fullCountsTable <- counts( object@countDataSet )
      cdsFull <- newCountDataSet( fullCountsTable, design )
      cdsFull <- estimateSizeFactors( cdsFull )
      object@cds <- estimateDispersions( cdsFull )
      if (compareCondition == "yes") {
        fit1 <- fitNbinomGLMs( object@cds, count ~ type + condition )
        fit0 <- fitNbinomGLMs( object@cds, count ~ type  )
      } else {
        fit1 <- fitNbinomGLMs( object@cds, count ~ type + condition )
        fit0 <- fitNbinomGLMs( object@cds, count ~ condition )
      }
      object@pval <- nbinomGLMTest( fit1, fit0 )
      object@padj <- p.adjust( object@pval, method="BH" )
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
  function(object,...) standardGeneric("smutans.de2TypeHeatmap")
)
setMethod("smutans.de2TypeHeatmap", 
          "Smutans", 
  function(object, file="default") {
    cdsBlind <- estimateDispersions( object@cds, method="blind" )
    vsd <- getVarianceStabilizedData( cdsBlind )
    select <- order(object@res$pval)[1:40]
    colors <- colorRampPalette(c("white","darkblue"))(100)

    if (file != "default")
    {
      postscript(file, width=6, height=6, horizontal = FALSE, onefile = FALSE, paper = "special")
      oldpar <- par (mar=c(4.5, 4.5, 0.5, 0.5))
    }
    heatmap( vsd[select,], col = colors, scale = "none")
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
  function(object, file="default") {
    if (file != "default")
    {
      postscript(file, width=6, height=6, horizontal = FALSE, onefile = FALSE, paper = "special")
      oldpar <- par (mar=c(4.5, 4.5, 0.5, 0.5))
    }
    cdsFullBlind <- estimateDispersions( object@cds, method = "blind" )
    vsdFull <- getVarianceStabilizedData( cdsFullBlind )
    dists <- dist( t( vsdFull ) )
    heatmap( as.matrix( dists ),
       symm=TRUE, scale="none", margins=c(10,10),
       col = colorRampPalette(c("darkblue","white"))(100),
       labRow = paste( pData(cdsFullBlind)$condition, pData(cdsFullBlind)$type ) )
    if (file != "default")
    {
      par(oldpar)
      dev.off()
    }
  }
)
