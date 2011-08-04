loadlibrary <- function ()
{
  library( DESeq )
  library( edgeR ) 
  library( hexbin ) 
  library( latticeExtra ) 
  library( gplots )
  library( geneplotter )
  sessionInfo()
}

loadfly <- function ()
{
  countsTableFly <- read.delim( "fly_RNA_counts.tsv" )
  condsFly <- c( "A", "A", "B", "B" )
  # Add dummy names to avoid confusion later
  rownames( countsTableFly ) <- paste( "Gene", 1:nrow(countsTableFly), sep="_" )

  cdsFly <- newCountDataSet( countsTableFly, condsFly )
  cdsFly <- estimateSizeFactors( cdsFly )
  cdsFly <- estimateVarianceFunctions( cdsFly )
  resFly <- nbinomTest( cdsFly, "A", "B" )
}

loadneural <- function ()
{
countsTableNS <- read.delim( "NeuralStemCellData.tab", row.names=1 )
condsNS <- c( "GNS_dupe", "GNS", "GNS", "GNSL", "NS", "NS" )
}

loadyeast <- function ()
{
countsYeast <- read.delim( "counts_Nagalakshmi_et_al.tab", row.names=1 ) 
}

runedgeR <- function ()
{
  # using DESeq's size factors:

  dglFlyD <- DGEList( counts=countsTableFly, group=condsFly, lib.size=1e7*sizeFactors(cdsFly) )
  dglFlyD <- estimateCommonDisp( dglFlyD )
  dglFlyD <- estimateTagwiseDisp( dglFlyD )

  # test with tag-wise dispersion:

  edgerResFlyDT <- exactTest( dglFlyD, common.disp=FALSE )
  edgerResFlyDTPadj <- p.adjust( edgerResFlyDT$table$p.value, method="BH" )

  # common dispersion, DESeq's size factors:

  edgerResFlyDC <- exactTest( dglFlyD, common.disp=TRUE )
  edgerResFlyDCPadj <- p.adjust( edgerResFlyDC$table$p.value, method="BH" )

  # with the original library sizes:

  # the total number of sequenced counts were as follows:
  libsizesFly <- c( A1=13613446, A2=12106517, B1=12920058, B2=11599459 )

  dglFlyN <- DGEList( counts=countsTableFly, group=condsFly, lib.size=libsizesFly )
  dglFlyN <- estimateCommonDisp( dglFlyN )
  dglFlyN <- estimateTagwiseDisp( dglFlyN )

  # test with tag-wise dispersion:

  edgerResFlyNT <- exactTest( dglFlyN, common.disp=FALSE )
  edgerResFlyNTPadj <- p.adjust( edgerResFlyNT$table$p.value, method="BH" )

  # common dispersion, DESeq's size factors:

  edgerResFlyNC <- exactTest( dglFlyN, common.disp=TRUE )
  edgerResFlyNCPadj <- p.adjust( edgerResFlyNC$table$p.value, method="BH" )
}

alpha <- dglFlyD$common.dispersion

baseVarFunc <- function( cds, cond ) {
   rvf <- rawVarFunc( cds, cond )
   sf <- sizeFactors(cds)[ conditions(cds) == cond ]
   xim <- sum(1/sf) / length(sf)
   function( q ) rvf( q ) + xim * q
}   

diagForA <- varianceFitDiagnostics( cdsFly, "A" )

xscale.components.log <- function( ... ) {
   res = xscale.components.default( ... )
   res$bottom$labels$labels = do.call( expression, 
      lapply( res$bottom$labels$at, function(a)
         substitute( 10^b, list(b=a) ) ) )
   res
}

yscale.components.log = function( ... ) {
   res = yscale.components.default( ... )    
   res$left$labels$labels = do.call( expression, 
      lapply( res$left$labels$at, function(a)
         substitute( 10^b, list(b=a) ) ) )
   res
}

xg <- seq( log10(1/3), 
   log10(max( diagForA$baseMean )), length.out=1000 )

var_mean_fly <- function ()
{
print( xyplot(    
   baseVar ~ baseMean, 
   diagForA[ diagForA$baseMean>0,], 
   panel = function( ... ) { 
      panel.hexbinplot( ... )
      
      # The solid orange line with the fitted base variance
      # (Note that we have to manually switch between natural
      # and logarithmic axis scaling, because hexbin cannot deal
      # with automatic logarithmic scaling):
      
      llines( xg, log10( baseVarFunc(cdsFly,"A")( 10^xg ) ),
         col="orange", lwd=2 )
	 
      # The dashed orange line with edgeR's variance fit,
      # i.e., variance = mean + alpha * mean^2:
      
      llines( xg, log10( 10^xg + alpha * (10^xg)^2 ),
         col="orange", lwd=2, lty="dashed" )
	 
      # For each sample with condition "A", plot a purple line
      # with the shot noise, variance = mean, where the mean
      # needs to be scaled by the size factor:
	 
      for( sf in sizeFactors( cdsFly )[ conditions(cdsFly) == "A" ] )
         llines( xg, log10( 10^xg*sf ), col="purple", lwd=1.5 ) },
	 
   trans=function(x) x^(1/6), 
   inv=function(x) x^6,
   xbins=80,
   scales = list( 
      x = list( log=TRUE, axs="i" ),
      y = list( log=TRUE, axs="i", limits=c( 3e-5,3e8 ), tick.number=8 ) ),
   xlab="mean", ylab="variance",
   xscale.components = xscale.components.log,
   yscale.components = yscale.components.log
) )  
}

# pdf( "cv2_mean_fly.pdf", width=4, height=4 )
cv2_mean_fly <- function ()
{
print(xyplot(    
   I( baseVar/baseMean^2 ) ~ baseMean, 
   diagForA[ diagForA$baseMean>0,], 
   panel = function( ... ) { 
      panel.hexbinplot( ... ) 
      
      # As we convert from a variance to an SCV, we should adjust for
      # bias, as described in Supplementary Note C. The function
      # adjustScvForBias takes care of this. It requires, as second
      # argument, the information how many samples were used for the
      # variance estimate. The 'attr(...)' expression supplies this,
      # returning in this case simply the value 2:
      
      llines( xg, 
         adjustScvForBias( 
	    baseVarFunc(cdsFly,"A")( 10^xg ) / (10^xg)^2, 
	    attr( rawVarFunc(cdsFly,"A")(NA), "size" ) ),
         col="orange", lwd=2 )
	 
      # The other lines are as before, only without the log of
      # the y value but with a division by mean^2.	 
	 
      llines( xg, ( 10^xg + alpha * (10^xg)^2 ) / (10^xg)^2,
         col="orange", lwd=2, lty="dashed" )
	 
      for( sf in sizeFactors( cdsFly )[ conditions(cdsFly) == "A" ] )
         llines( xg, 1 / (10^xg * sf), col="purple", lwd=1.5 ) },
	 
   trans = function(x) x^(1/4), 
   inv = function(x) x^4,
   xbins = 80,
   scales = list( 
      x = list( log=TRUE, axs="i" ),
      y = list( log=FALSE, axs="i", tick.number=8, limits=c(0,.2) ) ),
   xlab = "mean", ylab="squared coefficient of variation",
   xscale.components = xscale.components.log,
))
}

# pdf( "cv2_mean_full_fly.pdf", width=4, height=4 )
cv2_mean_full_fly <- function ()
{
print(xyplot(    
   I( baseVar/baseMean^2 ) ~ baseMean, 
   diagForA[ diagForA$baseMean>0,], 
   panel = function( ... ) { 
      panel.hexbinplot( ... ) 
      llines( xg, 
         adjustScvForBias( 
	    baseVarFunc(cdsFly,"A")( 10^xg ) / (10^xg)^2, 
	    attr( rawVarFunc(cdsFly,"A")(NA), "size" ) ),
         col="orange", lwd=2 )
      llines( xg, ( 10^xg + alpha * (10^xg)^2 ) / (10^xg)^2,
         col="orange", lwd=2, lty="dashed" )
      for( sf in sizeFactors( cdsFly )[ conditions(cdsFly) == "A" ] )
         llines( xg, 1 / (10^xg * sf), col="purple", lwd=1.5 ) },
   trans = function(x) x^(1/4), 
   inv = function(x) x^4,
   xbins = 80,
   scales = list( x = list( log=TRUE, axs="i" ),
      y = list( log=FALSE, axs="i", tick.number=8 ) ),
   xlab = "mean", ylab="squared coefficient of variation",
   xscale.components = xscale.components.log
))
}

# pdf( "resEcdf_A_fly.pdf", width=8, height=8 )
resEcdf_A_fly <- function ()
{
residualsEcdfPlot( cdsFly, "A" )
}

# pdf( "DE_fly.pdf", width=4, height=4 )
DE_fly <- function ()
{
print(xyplot( 
   log2FoldChange ~ I( baseMean ),
   resFly,  
   pch=16, cex=.3, 
   col=ifelse( resFly$padj < .1, "#FF000040", "#00000040" ),
   panel = function( x, y, col, ... ) {
      above <- (y > 5.8)
      below <- (y < -5.8)
      inside <- !( above | below )
      panel.xyplot( x=x[inside], y=y[inside], col=col[inside], ... ) 
      panel.arrows( x[above], 5.8, x[above], 5.95, col=col[above],
         length=".1", unit="native" )
      panel.arrows( x[below], -5.8, x[below], -5.95, col=col[below],
         length=".1", unit="native" ) },
   axis = function( side, ... ) {
      if( side=="left") {
         panel.axis( side, outside=TRUE, at=seq(-14,14,by=1), labels=FALSE )
         panel.axis( side, outside=TRUE, at=seq(-10,10,by=5), labels=TRUE )
      } 
      if( side=="bottom") {
         panel.axis( side, outside=TRUE, at=seq(-2,10,by=1), rot=0,
            labels = do.call( expression, 
	       lapply( seq(-2,10,by=1), function(a) 
	          substitute( 10^b, list(b=a) ) ) ) )
      } },
   xlab = "mean", ylab = "log2 fold change",
   scales = list( 
      x = list( log=TRUE ), 
      y = list( log=FALSE, limits=c( -6, 6 ) ) ) ))
}

# pdf( "vulcano_fly.pdf", width=4, height=4 )
vulcano_fly <- function ()
{
print(xyplot( -log10( pval ) ~ log2FoldChange, 
   resFly, 
   pch=20, cex=.2,
   col=ifelse( resFly$padj<.1, "#FF000050", "#00000050" ),
   axis = function( side, ... ) {
      if( side=="bottom") {
         panel.axis( side, outside=TRUE, at=seq(-14,14,by=1), labels=FALSE )
         panel.axis( side, outside=TRUE, at=seq(-10,10,by=5), labels=TRUE )
      }
      if( side=="left") {
         panel.axis( side, outside=TRUE, at=seq(0,25,by=1), labels=FALSE )
         panel.axis( side, outside=TRUE, at=seq(0,25,by=5), 
            labels = do.call( expression, 
	       lapply( seq(0,25,by=5), function(a) 
	          substitute( 10^-b, list(b=a) ) ) ) )
      } },
   xlab = "log2 fold change", ylab = "p value",
   scales = list( 
      x = list( limits=c( -6, 6 ) ),
      y = list( limits=c( 0, 25 ) ) ) ))
}

goodRows <- resFly$id %in% rownames(edgerResFlyDT$table)
densAll   <- density( log10(resFly$baseMean)[ goodRows ] )
testrows <- resFly$padj < .1
testrows[is.na(testrows)] <- FALSE
#densHits  <- density( log10(resFly$baseMean)[ goodRows & resFly$padj < .1 ] )
densHits  <- density( log10(resFly$baseMean)[ goodRows & testrows ] )
densedgDT <- density( log10(resFly$baseMean)[ goodRows ][ edgerResFlyDTPadj < .1 ] )
densedgDC <- density( log10(resFly$baseMean)[ goodRows ][ edgerResFlyDCPadj < .1 ] )
densedgNT <- density( log10(resFly$baseMean)[ goodRows ][ edgerResFlyNTPadj < .1 ] )
densedgNC <- density( log10(resFly$baseMean)[ goodRows ][ edgerResFlyNCPadj < .1 ] )

# pdf( "hits_dens_fly_all.pdf", width=6, height=4 )
hits_dens_fly_all <- function ()
{
plot( densAll$x, densAll$y / 7, type="l", ylim=c(0,.05),
   xlab="log10 mean", ylab="density", xaxs="i", yaxs="i", col="gray" )
lines( densHits$x, densHits$y * densHits$n / densAll$n, 
   col="red", lwd=2 )
lines( densedgDT$x, densedgDT$y * densedgDT$n / densAll$n, 
   col="blue", lwd=2, lty="dashed" )
lines( densedgDC$x, densedgDC$y * densedgDC$n / densAll$n, 
   col="lightblue", lwd=2, lty="dashed" )
lines( densedgNT$x, densedgNT$y * densedgNT$n / densAll$n, 
   col="blue", lwd=2, lty="solid" )
lines( densedgNC$x, densedgNC$y * densedgNC$n / densAll$n, 
   col="lightblue", lwd=2, lty="solid" )
text( densAll$x[310]+.2, densAll$y[310]/7, "x7", col="gray" )
}

# pdf( "hits_dens_fly.pdf", width=6, height=4 )
hits_dens_fly <- function ()
{
plot( densAll$x, densAll$y / 7, type="l", ylim=c(0,.05),
   xlab="log10 mean", ylab="density", xaxs="i", yaxs="i", col="gray" )
lines( densHits$x, densHits$y * densHits$n / densAll$n, 
   col="red", lwd=2 )
lines( densedgNT$x, densedgNT$y * densedgNT$n / densAll$n, 
   col="blue", lwd=2, lty="solid" )
lines( densedgNC$x, densedgNC$y * densedgNC$n / densAll$n, 
   col="lightblue", lwd=2, lty="solid" )
text( densAll$x[310]+.2, densAll$y[310]/7, "x7", col="gray" )
}

length(which( resFly$pval < .05 )) / nrow( resFly )
length( which( resFly$padj < .1 ))
c( DT = length( which( edgerResFlyDTPadj < .1 ) ),
   DC = length( which( edgerResFlyDCPadj < .1 ) ),
   NT = length( which( edgerResFlyNTPadj < .1 ) ),
   NC = length( which( edgerResFlyNCPadj < .1 ) ) )

list( 
  table( DESeq = resFly$padj[goodRows] < .1, edgeR_DT = edgerResFlyDTPadj < .1 ),
  table( DESeq = resFly$padj[goodRows] < .1, edgeR_DC = edgerResFlyDCPadj < .1 ),
  table( DESeq = resFly$padj[goodRows] < .1, edgeR_NT = edgerResFlyNTPadj < .1 ),
  table( DESeq = resFly$padj[goodRows] < .1, edgeR_NC = edgerResFlyDCPadj < .1 ) )

sameassame <- function ()
{
a <- counts(cdsFly)[,conditions(cdsFly)=="A"]
cdsA <- newCountDataSet( a[rowSums(a)>0,], c( "A1", "A2" ) )
cdsA <- estimateSizeFactors( cdsA )
cdsA <- estimateVarianceFunctions( cdsA, method="blind" )
resA <- nbinomTest( cdsA, "A1", "A2" )

dglA <- DGEList( counts(cdsA), sizeFactors( cdsA ) * 1e7, c( "A", "A" ) )
# Errors in the following two functions.
# dglA <- estimateCommonDisp( dglA ) 
# dglA <- estimateTagwiseDisp( dglA )
# dglA$samples$group = c( "A1", "A2" )
# eResA <- exactTest( dglA, common.disp=FALSE  ) 
}

# pdf( file="repl_pval_fly_A_ecdf.pdf", width=6, height=6)
repl_pval_fly_A_ecdf <- function ()
{
print(ecdfplot( ~ data | which, 
   make.groups( 
      `DESeq, below 100` = resA$pval[ resA$baseMean < 100 ],
      `DESeq, above 100` = resA$pval[ resA$baseMean > 100 ],   
      `DESeq, all` = resA$pval,
      `edgeR, below 100` = eResA$table$p.value[ 
          resA$baseMean < 100 ],
      `edgeR, above 100` = eResA$table$p.value[ 
          resA$baseMean > 100 ],
      `edgeR, all` = eResA$table$p.value,
      `Poisson, below 100` = pvalPois[ resA$baseMean < 100 ],
      `Poisson, above 100` = pvalPois[ resA$baseMean > 100 ],	
      `Poisson, all` = pvalPois ),
   layout=c(3,3), as.table=TRUE, xlab = "p value",
   panel = function( ...) {
      lsegments( 0, 0, 1, 1, col="darkgray" )
      panel.ecdfplot( ... ) },
   scales = list( tick.number = 4, cex=1 ) )) 
}

deseqneural <- function ()
{
cdsNS <- newCountDataSet( countsTableNS[,-1], condsNS[-1] )
cdsNS <- estimateSizeFactors( cdsNS )
cdsNS <- estimateVarianceFunctions( cdsNS )
resNS <- nbinomTest( cdsNS, "GNS", "NS" )

colGNSL <- which( conditions(cdsNS) == "GNSL" )

dglNS.D <- DGEList( 
   counts = counts(cdsNS)[,-colGNSL], 
   group = factor( conditions(cdsNS)[-colGNSL] ), 
   lib.size = 1e7*sizeFactors(cdsNS)[-colGNSL] )
dglNS.D <- estimateCommonDisp( dglNS.D )
dglNS.D <- estimateTagwiseDisp( dglNS.D )

edgerResNS.DT <- exactTest( dglNS.D, common.disp=FALSE )
edgerResNS.DTPadj <- p.adjust( edgerResNS.DT$table$p.value, method="BH" )

edgerResNS.DC <- exactTest( dglNS.D, common.disp=TRUE )
edgerResNS.DCPadj <- p.adjust( edgerResNS.DC$table$p.value, method="BH" )

libSizesNS <- c( G144=7604834, G166=13625570, G179=12291910,
   CB541=12872125, CB660=10502656 )

dglNS.N <- DGEList( 
   counts = counts(cdsNS)[,-colGNSL], 
   group = factor( conditions(cdsNS)[-colGNSL] ), 
   lib.size = libSizesNS[-colGNSL] )
dglNS.N <- estimateCommonDisp( dglNS.N )
dglNS.N <- estimateTagwiseDisp( dglNS.N )

edgerResNS.NT <- exactTest( dglNS.N, common.disp=FALSE )
edgerResNS.NTPadj <- p.adjust( edgerResNS.NT$table$p.value, method="BH" )

edgerResNS.NC <- exactTest( dglNS.N, common.disp=TRUE )
edgerResNS.NCPadj <- p.adjust( edgerResNS.NC$table$p.value, method="BH" )
}

print(length( which( resNS$padj < .1 ) ))
print(c( DT = length( which( edgerResNS.DTPadj < .1 ) ),
   DC = length( which( edgerResNS.DCPadj < .1 ) ),
   NT = length( which( edgerResNS.NTPadj < .1 ) ),
   NC = length( which( edgerResNS.NCPadj < .1 ) ) ))

# pdf( "cv2_mean.pdf", width=4, height=4 )
neural_cv2_mean <- function ()
{
diagForGNS <- varianceFitDiagnostics( cdsNS, "GNS" )
alpha <- dglNS.D$common.dispersion

print(xyplot(    
   I( baseVar/baseMean^2 ) ~ baseMean, 
   diagForGNS[ diagForGNS$baseMean>0,], 
   panel = function( ... ) { 
      panel.hexbinplot( ... ) 
      llines( 
         xg, 
         adjustScvForBias( 
           baseVarFunc(cdsNS,"GNS")( 10^xg ) / (10^xg)^2, 
           attr( rawVarFunc(cdsNS,"GNS")(NA), "size" ) ),
         col="orange", lwd=2 )
      llines( xg, ( 10^xg + alpha * (10^xg)^2 ) / (10^xg)^2,
         col="orange", lwd=2, lty="dashed" )
      for( sf in sizeFactors( cdsNS )[ conditions(cdsNS) == "GNS" ] )
         llines( xg, 1 / (10^xg * sf), col="purple", lwd=1.5 ) },
   trans = function(x) x^(1/4), 
   inv = function(x) x^4,
   xbins = 80,
   scales = list( x = list( log=TRUE, axs="i" ),
      y = list( log=FALSE, axs="i", tick.number=8 ) ),
   xlab = "mean", ylab="squared coefficient of variation"
))   
}

# pdf( "DE.pdf", width=4, height=4 )
neural_de <- function ()
{
print(xyplot( 
   log2FoldChange ~ I( baseMean ),
   resNS, 
   pch=16, cex=.3, 
   col=ifelse( resNS$padj < .1, "#FF000040", "#00000040" ),
   axis = function( side, ... ) {
      if( side=="left") {
         panel.axis( side, outside=TRUE, at=seq(-14,14,by=1), labels=FALSE )
         panel.axis( side, outside=TRUE, at=seq(-10,10,by=5), labels=TRUE )
      } 
      if( side=="bottom") {
         panel.axis( side, outside=TRUE, at=seq(-2,10,by=1), 
            labels=paste( "10^", seq(-2,10,by=1), sep="" ) )
      } },
   xlab = "mean", ylab = "log2 fold change",
   scales = list( x=list(log=TRUE), y=list(log=FALSE) ) ))
}


# pdf( "hits_dens.pdf", width=6, height=4 )
neural_hits_dens <- function ()
{
goodRowsNS <- resNS$id %in% rownames(edgerResNS.DT$table)

densAll.NS   <- density( log10(resNS$baseMean)[ goodRowsNS ] )
testrowsNS <- resNS$padj < .1
testrowsNS[is.na(testrowsNS)] <- FALSE
densHits.NS  <- density( log10(resNS$baseMean)[ goodRowsNS & testrowsNS ] )
densedgDT.NS <- density( log10(resNS$baseMean)[ goodRowsNS ][ edgerResNS.DTPadj < .1 ] )
densedgDC.NS <- density( log10(resNS$baseMean)[ goodRowsNS ][ edgerResNS.DCPadj < .1 ] )
densedgNT.NS <- density( log10(resNS$baseMean)[ goodRowsNS ][ edgerResNS.NTPadj < .1 ] )
densedgNC.NS <- density( log10(resNS$baseMean)[ goodRowsNS ][ edgerResNS.NCPadj < .1 ] )


plot( densAll.NS$x, densAll.NS$y / 7, type="l", ylim=c(0,.07),
   xlab="log10 mean", ylab="density", xaxs="i", yaxs="i", 
   col="gray" )
lines( densHits.NS$x, densHits.NS$y * densHits.NS$n / densAll.NS$n, 
   col="red", lwd=2 )
lines( densedgDT.NS$x, densedgDT.NS$y * densedgDT.NS$n / densAll.NS$n, 
   col="blue", lwd=2, lty="dashed" )
lines( densedgDC.NS$x, densedgDC.NS$y * densedgDC.NS$n / densAll.NS$n, 
   col="lightblue", lwd=2, lty="dashed" )
lines( densedgNT.NS$x, densedgNT.NS$y * densedgNT.NS$n / densAll.NS$n, 
   col="blue", lwd=2, lty="solid" )
lines( densedgNC.NS$x, densedgNC.NS$y * densedgNC.NS$n / densAll.NS$n, 
   col="lightblue", lwd=2, lty="solid" )
text( densAll.NS$x[310]+.2, densAll.NS$y[310]/7, "x7", col="gray" )
}


# pdf( file="repl_pval_A_ecdf.pdf", width=6, height=4.5)
neural_repl_pval_A_ecdf <- function ()
{
cdsGNS <- newCountDataSet( counts(cdsNS)[,1:2], c( "GNS_A", "GNS_B" ) )
cdsGNS <- estimateSizeFactors( cdsGNS )
cdsGNS <- estimateVarianceFunctions( cdsGNS, method="blind" )
resGNS <- nbinomTest( cdsGNS, "GNS_A", "GNS_B" )
medGNS <- median( getBaseMeansAndVariances( counts(cdsGNS), 
             sizeFactors(cdsGNS) )$baseMean )

dglGNS <- DGEList( counts(cdsGNS), sizeFactors( cdsGNS ) * 1e7, 
   c( "GNS", "GNS" ) )
dglGNS <- estimateCommonDisp( dglGNS )
dglGNS <- estimateTagwiseDisp( dglGNS )
dglGNS$samples$group = c( "GNS_1", "GNS_2" )
eResGNS <- exactTest( dglGNS, common.disp=FALSE ) 
eMedGNS <- median( rowMeans(dglGNS$pseudo.alt) )

print(ecdfplot( ~ data | which, 
   make.groups( 
      `DESeq lower half` = resGNS$pval[ resGNS$baseMean < medGNS ],
      `DESeq upper half` = resGNS$pval[ resGNS$baseMean > medGNS ],   
      `DESeq all` = resGNS$pval,
      `edgeR lower half` = eResGNS$table$p.value[ 
          rowMeans(dglGNS$pseudo.alt) < eMedGNS ],
      `edgeR upper half` = eResGNS$table$p.value[ 
          rowMeans(dglGNS$pseudo.alt) > eMedGNS ],
      `edgeR all` = eResGNS$table$p.value ),
   layout=c(3,2), xlab = "p value", as.table=TRUE, 
   panel = function( ...) {
      lsegments( 0, 0, 1, 1, col="darkgray" )
      panel.ecdfplot( ... ) },
   scales = list( tick.number=4, cex=1 ) ))
}

withreplicate <- function ()
{
cdsFly2 <- cdsFly[ ,c( "A1", "B1" ) ]
cdsFly2 <- estimateVarianceFunctions( cdsFly2, method="blind" )
resFly2 <- nbinomTest( cdsFly2, "A", "B" )
print(addmargins( table( full = resFly$padj < .1, reduced = resFly2$padj < .1 )))
cdsNS2 <- cdsNS[ ,c( "G144", "CB541" ) ]
cdsNS2 <- estimateVarianceFunctions( cdsNS2, method="blind" )
resNS2 <- nbinomTest( cdsNS2, "GNS", "NS" )
print(addmargins( table( full = resNS$padj < .1, reduced = resNS2$padj < .1 ) ))
}

# pdf( "heatmap.pdf" )
heatmap <- function ()
{
cdsNS.all <- newCountDataSet( countsTableNS, rep( "dummy", ncol(countsTableNS) ) )
cdsNS.all <- estimateSizeFactors( cdsNS.all )

vsdNS <- getVarianceStabilizedData( cdsNS.all )
distMatrix <- as.matrix(dist(t( vsdNS )))
rownames( distMatrix ) <- c( "GNS (*)", "GNS (*)", "GNS", "GNS (L)", "NS", "NS" )

heatmap.2( distMatrix, density.info="none", trace="none", revC=TRUE,
   col=colorRampPalette(c("darkblue","yellow"))(100) )
}

yeastreplicateonly <- function ()
{
sfYeast <- estimateSizeFactorsForMatrix( countsYeast )

# A function to get a base variance function:
get_bvf <- function( cl ) {
   rvf <- estimateVarianceFunctionForMatrix( countsYeast[,cl], sfYeast[cl] )
   xim <- sum(1/sfYeast[cl])/length(cl)
   function( q ) rvf( q ) + xim * q
}   

# The variance estimated between the replicates:
bvf_dT_tech <- get_bvf( c("dT_ori","dT_tech") )
bvf_dT_bio  <- get_bvf( c("dT_ori","dT_bio" ) )
bvf_RH_tech <- get_bvf( c("RH_ori","RH_tech") )
bvf_RH_bio  <- get_bvf( c("RH_ori","RH_bio" ) )

# The common-scale means
meansYeast <- getBaseMeansAndVariances( countsYeast, sfYeast )$baseMean
layout( cbind(1:2), heights = c( 1, 3 ) )
par( cex=.75 )

par(mai=c(0.2, 1, 0.5, 0.5))
plot( NULL, main="", xaxs="i", yaxs="i", ylab="density", 
   xlim=c( -.3, 4 ), ylim=c( 0, 1200 ), xaxt="n" )
hist( log10( meansYeast ), breaks=30, add=TRUE )

par(mai=c(1, 1, 0.01, 0.5))
plot( NULL, log="x", ylim=c(0, .65), xaxs="i", yaxs="i",
   xlab="mean", ylab="squared coefficient of variation",
   xlim=c( 10^-.3, 10^4 ) )
xg <- 10^seq( 0, 4, length.out=1000 ) 
rect( xg[-1000], 1/(xg[-1000]*max(sfYeast)), 
   xg[-1], 1/(xg[-1000]*min(sfYeast)), col="#FFD0FF", lty="blank" )
lines( xg, bvf_dT_bio(xg)/xg^2, col="brown" )
lines( xg, bvf_dT_tech(xg)/xg^2, col="darkblue" )
lines( xg, bvf_RH_bio(xg)/xg^2, col="brown", lty="dashed" )
lines( xg, bvf_RH_tech(xg)/xg^2, col="darkblue", lty="dashed" )
lines( xg, 1/xg, col="#FFD0FF" )
}

# I skip ChIP-Seq

sim1 <- function ()
{
# Number of genes to simulate
ngenes <- 20000

# True mean expression value across conditions drawn from an
# exponential distribution
true_q0middle <- rexp( ngenes, 1/250 )

# These values might be typical size factors as could be found
# in a real experiment
sf <- c( .5, 1.7, 1.4, .9 )

# Two replicates for each condition:
conds <- c( "A", "A", "B", "B" )

# Mark around 30% of the genes as differentially expressed:
true_isDE <- runif( ngenes ) < .3

# For the DE genes, draw a log fold change from a normal
true_lfc <- rnorm( ngenes, sd=.7 ) * true_isDE

# Go up and down from the true_q0middle according to these log
# fold changes to get the values q for the conditions
true_q0 <- data.frame(
   A = true_q0middle * 2^(  true_lfc/2 ),
   B = true_q0middle * 2^( -true_lfc/2 ) )
   
# For the variance-mean relation, let's assume that the shape 
# assumed by edgeR is correct, i.e., variance = mean + alpha * mean^2
# and use this alpha:   
alpha <- .015   

# Parametrise the NB to mean and variance:
rnbinomMV <- function( n, mu, v ) 
  rnbinom( n, prob = mu/v, size = mu^2/(v-mu) )

# Now draw counts
countsTable <- 
   sapply( 1:4, function(j)
      sapply( true_q0[[ conds[j] ]], function(q)
         rnbinomMV( 1, sf[j]*q, sf[j]*q + sf[j]^2 * q^2 * alpha ) ) )
	 
# The good rows:
nonzero <- rowSums(countsTable)>0 
# A standard analysis with DESeq:
cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
res <- nbinomTest( cds, "A", "B" )

# The same with edgeR. (edgeR gets the advantage of knowing the true
# library sizes)
dgl <- DGEList( counts=countsTable, group=conds, lib.size=sf*1e7 )
dgl <- estimateCommonDisp(dgl)
edgerRes <- exactTest( dgl )
edgerResPadj <- p.adjust( edgerRes$table$p.value, method="BH" )
}

sim2 <- function ()
{
print(length(which( res$pval[nonzero]<.05 & !true_isDE[nonzero] )) / length(nonzero))
print(length(which( edgerRes$table$p.value<.05 & !true_isDE[nonzero] )) / length(nonzero))
print(length(which( res$padj[nonzero]<.1 & true_isDE[nonzero] )) / 
   length(which( true_isDE[nonzero] )))
print(length(which( edgerResPadj<.1 & true_isDE[nonzero] )) / 
   length(which( true_isDE[nonzero] )))

print(length(which( res$padj[nonzero]<.1 & !true_isDE[nonzero] )) / 
   length(which( res$padj[nonzero]<.1 )))
print(length(which( edgerResPadj<.1 & !true_isDE[nonzero] )) / 
   length(which( edgerResPadj<.1 )))


}

# pdf( "pval_DE_B.pdf", width=4, height=4 )
pval_DE_B <- function ()
{
print(ecdfplot( ~ data | which, 
   make.groups( 
      `DESeq` = res$pval[nonzero][ !true_isDE[nonzero] ],
      `edgeR` = edgerRes$table$p.value[ !true_isDE[nonzero] ] ),
   layout=c(2,1), as.table=TRUE, xlab = "p value",
   panel = function( ...) {
      lsegments( 0, 0, 1, 1, col="darkgray" )
      panel.ecdfplot( ... ) },
   xlim=c( -.005, .1 ), ylim=c( -.005, .1 ),
   scales = list( tick.number = 4, cex=1 ) )) 
}

# pdf( "pval_DE.pdf", width=4, height=4 )
pval_DE <- function ()
{
print(xyplot(  
   res$pval[true_isDE & nonzero] ~ edgerRes$table$p.value[true_isDE[nonzero]],
   panel = function( ... ) {
      panel.hexbinplot( ... )
      llines( c(-35, 0), c(-35, 0), col="#30FF5090" ) },
   trans = function(x) x^(1/4), 
   inv = function(x) x^4,
   xbins = 150,
   scales = list(
      x = list(log=TRUE),
      y = list(log=TRUE), apect=1 ),
   xlim = c( 1e-15, 10 ),
   ylim = c( 1e-15, 10 ),
   aspect = 1,
   xlab="p value from edgeR", ylab="p value from DESeq", asp=1
))
}

alpha <- function(mu) .01 + 9/(mu+100)
sim3 <- function ()
{
ngenes <- 20000
true_q0middle <- rexp( ngenes, 1/250 )
sf <- c( .5, 1.7, 1.4, .9 )
conds <- c( "A", "A", "B", "B" )
true_isDE <- runif( ngenes ) < .3
true_lfc <- rnorm( ngenes, sd=.7 ) * true_isDE
true_q0 <- data.frame(
   A = true_q0middle * 2^(  true_lfc/2 ),
   B = true_q0middle * 2^( -true_lfc/2 ) )

countsTable <- 
   sapply( 1:4, function(j)
      sapply( true_q0[[ conds[j] ]], function(q)
         rnbinomMV( 1, sf[j]*q, sf[j]*q + sf[j]^2 * q^2 * alpha( q ) ) ) )

nonzero <- rowSums(countsTable)>0 

cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds )
res <- nbinomTest( cds, "A", "B" )

dgl <- DGEList( counts=countsTable, group=conds, lib.size=sf*1e7 )
dgl <- estimateCommonDisp(dgl)
dgl <- estimateTagwiseDisp(dgl)
edgerRes <- exactTest( dgl, common.disp=FALSE )
edgerResPadj <- p.adjust( edgerRes$table$p.value, method="BH" )
}

scvSim2 <- function ()
{
scvPlot( cds )
xg <- 10^seq( -2, 4, length.out=1000 )
lines( xg, alpha( xg ), col="purple" )
}

# pdf( file="pvalcompsim2", width=8, height=6 )
pvalcompsim2 <- function ()
{
baseMeans <- t(t(countsTable)/sf)
lowerHalf <- baseMeans < median( baseMeans )

print(ecdfplot( ~ data | which, 
   make.groups( 
      `DESeq, below median` = res$pval[nonzero][ 
         !true_isDE[nonzero] & lowerHalf[nonzero] ],
      `DESeq, above median` = res$pval[nonzero][ 
         !true_isDE[nonzero] & !lowerHalf[nonzero] ],   
      `DESeq, all`          = res$pval[nonzero][ 
         !true_isDE[nonzero] ],
      `edgeR, below median` = edgerRes$table$p.value[ 
         !true_isDE[nonzero] & lowerHalf[nonzero] ],
      `edgeR, above median` = edgerRes$table$p.value[ 
         !true_isDE[nonzero] & !lowerHalf[nonzero] ],
      `edgeR, all` = edgerRes$table$p.value[ 
         !true_isDE[nonzero] ] ),
   layout=c(3,2), as.table=TRUE, xlab = "p value",
   panel = function( ...) {
      lsegments( 0, 0, 1, 1, col="darkgray" )
      panel.ecdfplot( ... ) },
   xlim=c( -.005, .1 ), ylim=c( -.005, .1 ),
   scales = list( tick.number = 4, cex=1 ) )) 
}

# pdf("biasAdjust.pdf", width=4, height=4)
biasAdjust <- function ()
{
plot( NULL, xlim=c(0,2), ylim=c(0,2), 
   xlab = expression( ( hat(sigma)^2 - hat(mu) ) / hat(mu)^2 ), 
   ylab = expression( gamma ) )
xg <- seq( 0, 2, length.out=1000 )
for( i in 2:15 ) 
   lines( xg, adjustScvForBias( xg, i ), 
      col=colorRampPalette(c("red","black"))(15)[i] )
}

# pdf("p1.pdf", width=4, height=4)
p1 <- function ()
{
dnbinomMV <- function( x, mu, v ) 
  dnbinom( x, prob = mu/v, size = mu^2/(v-mu) )

plot( 0:10000, 
   sapply( 0:10000, function(k) 
      dnbinomMV( k, 7000, 7000+.1*7000^2 ) * 
         dnbinomMV( 10000-k, 4000, 4000+.1*4000^2 ) ), 
   type='l', xlab="a", ylab="p(a,b)" )
abline( v = 7000 / (7000+4000) * 10000 )   
}

###################################################################
# From the manual
#library( "DESeq" )
#library( "pasilla" )
#data( "pasillaGenes" )


