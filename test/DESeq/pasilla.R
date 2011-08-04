#library( "DESeq" )

# Could not find this package.
# library( "pasilla" )
# data( "pasillaGenes" )

loadfly <- function ()
{
  countsTableFly <- read.delim( "fly_RNA_counts.tsv" )

  # Q: Which one is right?
  condsFly <- c( "A", "A", "B", "B" )
  # condsFly <- factor(c( "A", "A", "B", "B" ))

  # Add dummy names to avoid confusion later
  rownames( countsTableFly ) <- paste( "Gene", 1:nrow(countsTableFly), sep="_" )
}

deseq1_fly <- function ()
{
  cdsFly <- newCountDataSet( countsTableFly, condsFly )
  cdsFly <- estimateSizeFactors( cdsFly )
  cdsFly <- estimateVarianceFunctions( cdsFly )
}

deseq2_fly <- function ()
{
  resFly <- nbinomTest( cdsFly, "A", "B" )
}

# Q: Not working
# ls(cdsFly@fitInfo)
plotDispEsts <- function( cds, cond )
{
   plot( 
      rowMeans( counts( cds, normalized=TRUE ) ),
      cds@fitInfo[[cond]]$perGeneDispEsts,
      pch = '.', log="xy" )
   xg <- 10^seq( -.5, 5, length.out=300 )
   lines( xg, cds@fitInfo[[cond]]$dispFun( xg ), col="red" )   
}

# plotDispEsts ( cdsFly, "A")
# print(head(fData(cdsFly)))


plotDE <- function( res )
   plot( 
      res$baseMean, 
      res$log2FoldChange, 
      log="x", pch=20, cex=.3, 
      col = ifelse( res$padj < .1, "red", "black" ) )

# plotDE( resFly )

res <- resFly
cds <- cdsFly
# Histogram of p values
#hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")

# Filter for significant genes at FDR 10%
resSig <- res[ res$padj < 0.1, ]
#print(head( resSig[ order(resSig$pval), ] ))

# Downregulated
#print(head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] ))
# Upregulated
#print(head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] ))

# Write the result.
#write.table( res, file="results.txt" )

# Compare within a treatment.
# Q: normalized is not working.
# ncu <- counts( cdsFly, normalized=TRUE )[ , conditions(cds)=="A" ]
# -------------------------------------------------------------------
#ncu <- counts( cdsFly )[ , conditions(cdsFly)=="A" ]
#plot( rowMeans(ncu), log2( ncu[,2] / ncu[,1] ), pch=".", log="x" )

partialreplicate <- function ()
{
cdsTTU <- cds[ , 1:3] 
pData( cdsTTU )
cdsTTU <- estimateSizeFactors( cdsTTU )
cdsTTU <- estimateVarianceFunctions( cdsTTU )
resTTU <- nbinomTest( cdsTTU, "B", "A" )
plotDE( resTTU )
}

noreplicate <- function ()
{
cds2 <- cds[ ,c("A1", "B1") ]
# Q: estimateDispersions no such function
# cds2 <- estimateDispersions( cds2, method="blind" )
cds2 <- estimateSizeFactors( cds2 )
cds2 <- estimateVarianceFunctions( cds2, method="blind" )
res2 <- nbinomTest( cds2, "A", "B" )
plotDE( res2 )
print(addmargins( table( res_sig = res$padj < .1, res2_sig = res2$padj < .1 ) ))
}

pseudocount <- function ()
{
cdsBlind <- estimateVarianceFunctions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cdsBlind )
mod_lfc <- (rowMeans( vsd[, conditions(cds)=="A", drop=FALSE] ) - 
            rowMeans( vsd[, conditions(cds)=="B", drop=FALSE] ))
lfc <- res$log2FoldChange
finite <- is.finite(lfc)
print(table(as.character(lfc[!finite]), useNA="always"))
largeNumber <- 10
lfc <- ifelse(finite, lfc, sign(lfc) * largeNumber)
plot( lfc, mod_lfc, pch=20, cex=.3,
      col = ifelse( finite, "#80808040", "red" ) )
abline( a=0, b=1, col="#40404040" )
}

clustering <- function ()
{
select <- order(res$pval)[1:40]
colors <- colorRampPalette(c("white","darkblue"))(100)
heatmap( vsd[select,], 
         col = colors, scale = "none")
heatmap( counts(cds)[select,], 
         col = colors, scale = "none")
}

clustering2 <- function ()
{
# cdsFullBlind <- estimateDispersions( cdsFull, method = "blind" )
# cdsFullBlind <- estimateDispersions( cdsFull, method = "blind" )
# vsdFull <- getVarianceStabilizedData( cdsFullBlind )
dists <- dist( t( vsd ) )
heatmap( as.matrix( dists ), 
   symm=TRUE, scale="none", margins=c(10,10),
   col = colorRampPalette(c("darkblue","white"))(100),
   labRow = paste( pData(cdsBlind)$condition, pData(cdsBlind)$type ) )
}
