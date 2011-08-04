estimateSizeFactors <- function( cds, locfunc=median )
{
   stopifnot( is( cds, "CountDataSet" ) )
   sizeFactors(cds) <- estimateSizeFactorsForMatrix( counts(cds), locfunc )
   cds
}

estimateVarianceFunctions <- function( cds, 
   method = c( "normal", "blind", "pooled" ), pool=NULL, 
   locfit_extra_args=list(), lp_extra_args=list(), modelFrame = NULL )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   if( any( is.na( sizeFactors(cds) ) ) )
      stop( "NAs found in size factors. Have you called already 'estimateSizeFactors'?" )
   
   if( length(method) != 3 & !is.null( pool ) )
      stop( "Do not specify both the 'pool' and the 'method' argument." )      
   if( !is.null( pool ) ) {
      if( pool == FALSE )
         method <- "normal"
      else if( pool == TRUE )
         method <- "blind"
      else
         stop( "Argument 'pool' used incorrectly." )
      warning( "The 'pool' argument to 'estimatevarianceFunction' is deprecated. Use the 'method' argument instead." )
   } else
      method <- match.arg( method )
   
   if( cds@multivariateConditions && ! method %in% c( "blind", "pooled" ) )
      stop( "You have specified multivariate conditions (i.e., passed a data frame with conditions). In this case, you need to specify ethod 'pooled' or 'blind'." )
   
   cds@rawVarFuncs <- new.env( hash=TRUE )
   
   if( method == "blind" ) {
      cds@rawVarFuncs[["_blind"]] <-
         estimateVarianceFunctionForMatrix( counts(cds), 
	    sizeFactors(cds), locfit_extra_args, lp_extra_args )
      if( cds@multivariateConditions )
         rawVarFuncTable(cds) <- c( "_all" = "_blind" )
      else {
         a <- rep( "_blind", length( levels( conditions(cds) ) ) )
	 names(a) <- levels( conditions(cds) )
	 rawVarFuncTable(cds) <- a } }

   else if( method == "normal" ) {
      replicated <- names( which( tapply( conditions(cds), conditions(cds), length ) > 1 ) )
      if( length( replicated ) < 1 )
         stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions." )
      nonreplicated <- names( which( tapply( conditions(cds), conditions(cds), length ) == 1 ) )
      for( cond in replicated )
         cds@rawVarFuncs[[cond]] <- estimateVarianceFunctionForMatrix( 
            counts(cds)[ , conditions(cds)==cond ], sizeFactors(cds)[ conditions(cds)==cond ], 
	       locfit_extra_args, lp_extra_args )
      cds@rawVarFuncs[["_max"]] <- function( q ) {
         a <- lapply( replicated, function(cond) cds@rawVarFuncs[[cond]]( q ) )
         ans <- apply( array( unlist( a, recursive=FALSE ), dim = c( length(q), length(a) ) ), 1, max )
         rownames(ans) <- rownames(q)
         attr( ans, "size" ) <- min( sapply( a, attr, "size" ) )
         ans
      }
         
      rawVarFuncTable(cds) <- sapply( levels(conditions(cds)), function( cond )
            ifelse( cond %in% replicated, cond, "_max" ) ) }

   else if( method == "pooled" ) {
   
      if( cds@multivariateConditions ) {
         if( is.null( modelFrame ) )
            modelFrame <- pData(cds)[ , colnames(pData(cds)) != "sizeFactor" ]
         conds <- modelMatrixToConditionFactor( modelFrame ) }
      else
         conds <- conditions(cds)

      cds@rawVarFuncs[["_pooled"]] <- estimatePooledVarianceFunctionForMatrix( 
         counts(cds), sizeFactors(cds), conds,
         locfit_extra_args, lp_extra_args )
         
      if( cds@multivariateConditions )
         rawVarFuncTable(cds) <- c( "_all" = "_pooled" )
      else {
         a <- rep( "_pooled", length( levels( conditions(cds) ) ) )
	 names(a) <- levels( conditions(cds) )
	 rawVarFuncTable(cds) <- a } }
        
   validObject( cds )
   cds
}

varianceFitDiagnostics <- function( cds, cond=NULL, modelFrame=NULL )
{
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )
   
   if( is.null( cond ) ) {
      if( length(cds@rawVarFuncs) != 1 )
         stop( "There is more than one variance function. 'cond' may not be omitted." )
      cond <- ls( cds@rawVarFuncs )[[1]]
      rvf <- rawVarFunc( cds, cond, byName=TRUE ) 
   } 
   else
      rvf <- rawVarFunc( cds, cond ) 
   
   if( cond == "_pooled" ) {
      if( cds@multivariateConditions ) {
         if( is.null( modelFrame ) )
            modelFrame <- pData(cds)[ , colnames(pData(cds)) != "sizeFactor" ]
         poolconds <- modelMatrixToConditionFactor( modelFrame ) }
      else
         poolconds <- conditions(cds)

      varianceFitDiagnosticsForMatrix(
         counts(cds), sizeFactors(cds), rvf, poolconds )
   }
   
   else if( cond == "_blind" )
      varianceFitDiagnosticsForMatrix(
         counts(cds), sizeFactors(cds), rvf )
   
   else 
      varianceFitDiagnosticsForMatrix(
         counts(cds)[,conditions(cds)==cond], 
         sizeFactors(cds)[conditions(cds)==cond],
         rvf ) 
}

residualsEcdfPlot <- function( cds, condition=NULL, ncuts=7 )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   ensureHasVarFuncs( cds )
   fitdiag <- varianceFitDiagnostics( cds, condition )
   residualsEcdfPlotFromDiagnostics( fitdiag, ncuts,
      sprintf( "Residuals ECDF plot for condition '%s'", condition ) )
}  

nbinomTest <- function( cds, condA, condB, pvals_only=FALSE, eps=1e-4 )
{
   stopifnot( is( cds, "CountDataSet" ) )   
   ensureHasVarFuncs( cds )
   if( cds@multivariateConditions )
      stop( "For CountDataSets with multivariate conditions, only the GLM-based test can be used." )
   stopifnot( condA %in% levels(conditions(cds)) )  
   stopifnot( condB %in% levels(conditions(cds)) )     
   
   colA <- conditions(cds)==condA
   colB <- conditions(cds)==condB

   bmv <- getBaseMeansAndVariances( counts(cds)[,colA|colB], 
      sizeFactors(cds)[colA|colB] )

   rvfA <- rawVarFunc( cds, condA )
   rvfB <- rawVarFunc( cds, condB )
   
   rawScvA <- rvfA( bmv$baseMean ) / bmv$baseMean^2
   rawScvB <- rvfB( bmv$baseMean ) / bmv$baseMean^2
   
   rawScvA <- adjustScvForBias( rawScvA, attr( rawScvA, "size" ) )
   rawScvB <- adjustScvForBias( rawScvB, attr( rawScvB, "size" ) )

   pval <- nbinomTestForMatrices( 
      counts(cds)[,colA], 
      counts(cds)[,colB], 
      sizeFactors(cds)[colA], 
      sizeFactors(cds)[colB], 
      rawScvA, 
      rawScvB,
      eps )
      
   if( pvals_only )
      pval
   else {
      bmvA <- getBaseMeansAndVariances( counts(cds)[,colA], sizeFactors(cds)[colA] )
      bmvB <- getBaseMeansAndVariances( counts(cds)[,colB], sizeFactors(cds)[colB] )
      data.frame( 
         id    = rownames( counts(cds) ),
         baseMean  = bmv$baseMean,
         baseMeanA = bmvA$baseMean,
         baseMeanB = bmvB$baseMean,
         foldChange = bmvB$baseMean / bmvA$baseMean,
         log2FoldChange = log2( bmvB$baseMean / bmvA$baseMean ), 
         pval = pval,
         padj = p.adjust( pval, method="BH" ), 
         resVarA = bmvA$baseVar / ( bmvA$baseMean * sum( 1/sizeFactors(cds)[colA] ) / length(condA) +
            rawVarFunc( cds, condA )( bmv$baseMean ) ),
         resVarB = bmvB$baseVar / ( bmvB$baseMean * sum( 1/sizeFactors(cds)[colB] ) / length(condB) +
            rawVarFunc( cds, condB )( bmv$baseMean ) ),
         stringsAsFactors = FALSE ) }
}

scvPlot <- function( cds, xlim=NULL, ylim=c(0,.8), skipBiasCorrection = FALSE ) {
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )
   
   baseMeans <- getBaseMeansAndVariances( counts(cds), sizeFactors(cds) )$baseMean

   xg <- exp( seq( log( max( min(baseMeans), 2/sum(sizeFactors(cds)) ) ), 
      log( max(baseMeans) ), length.out = 1000 ) )
   if( is.null( xlim ) )
      xlim <- range( xg )
   plot( NULL, xlim = xlim, ylim = ylim, yaxs="i",
      xlab = "base mean", ylab = "squared coefficient of variation", log="x" )
      
   color <- 2
   legendData <- data.frame( text=character(0), col=integer(0), lty=character(0), 
      stringsAsFactors=FALSE ) 
   for( funcName in unique( rawVarFuncTable(cds) ) ) {
   
      rawScv <- rawVarFunc( cds, funcName, byName=TRUE )( xg ) / xg^2      

      lty <- if( funcName != "_max" ) "solid" else "dotted"
      lines( xg, rawScv, col=color, lty=lty )
      legendData <- rbind( legendData, data.frame( text=funcName, col=color, lty=lty,
         stringsAsFactors=FALSE ) )
	 
      for( cond in names(rawVarFuncTable(cds))[ rawVarFuncTable(cds) == funcName ] ) {
         if( cond != "_all" )
	    cols <- which( conditions(cds) == cond )
	 else
	    cols <- 1:ncol(counts(cds))
	 for( j in cols )
	    lines( xg, rawScv + 1 / ( sizeFactors(cds)[[j]] * xg ), col=color, lty="dashed" )
      }

      color <- color + 1 
   }
   
   dens <- density( log(baseMeans) )
   lines( exp(dens$x), .7 * ylim[2] / max(dens$y) * dens$y, col=1, lty="solid" )
   legendData <- rbind( legendData, data.frame( text="base mean density", col=1, lty="solid",
      stringsAsFactors=FALSE ) )   

   legend( "topright", 
      legend = legendData$text,
      col    = legendData$col,
      lty    = legendData$lty )

   invisible( NULL )
}

getVarianceStabilizedData <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( any( is.na( sizeFactors(cds) ) ) )
      stop( "NA found in size factors. Have you called 'estimateSizefactors' yet?" )
   ncounts <- t( t(counts(cds)) / sizeFactors(cds) )
   rvf <- estimateVarianceFunctionForMatrix( counts(cds), sizeFactors(cds) )
   xg <- sinh( seq( asinh(0), asinh(max(ncounts)), length.out=1000 ) )[-1]
   xim <- mean( 1/sizeFactors(cds) )
   baseVarsAtGrid <- xg * xim^2 + pmax( 0, rvf( xg ) )
   integrand <- 1 / sqrt( baseVarsAtGrid )
   splf <- splinefun( 
      asinh( ( xg[-1] + xg[-length(xg)] )/2 ), 
      cumsum( 
         ( xg[-1] - xg[-length(xg)] ) * 
         ( integrand[-1] + integrand[-length(integrand)] )/2 ) )
   tc <- sapply( colnames(counts(cds)), function(clm)
      splf( asinh( ncounts[,clm] ) ) )
   rownames( tc ) <- rownames( counts(cds) )
   tc
}

makeExampleCountDataSet <- function( ) 
{
   ngenes <- 10000
   q0 <- rexp( ngenes, rate=1/250 )
   is_DE <- runif( ngenes ) < .3
   lfc <- rnorm( ngenes, sd=2 )
   q0A <- ifelse( is_DE, q0 * 2^(  lfc/2 ), q0 )
   q0B <- ifelse( is_DE, q0 * 2^( -lfc/2 ), q0 )
   true_sf <- c( 1., 1.3, .7, .9, 1.6 )   
   conds <- c( "A", "A", "B", "B", "B" )
   m <- t( sapply( 1:ngenes, function(i) 
      sapply( 1:5, function( j )
         rnbinom( 1, mu = true_sf[j] * ifelse( conds[j]=="A", q0A[i], q0B[i] ), 
            size = 1/.2 ) ) ) )
   colnames(m) <- c( "A1", "A2", "B1", "B2", "B3" )
   rownames(m) <- paste( "gene", 1:ngenes, 
      ifelse( is_DE, "T", "F" ), sep="_" )
   newCountDataSet( m, conds )
}

nbinomFitGLM <- function( cds, modelFormula, glmControl=list() )
{
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )
   if( is.null( cds@rawVarFuncs[["_pooled"]] ) )
      stop( "No pooled variance function found. Have you called 'estimateVarianceFunctions' with 'method=\"pooled\"'?" )
      
   baseMeans <- colMeans(t(counts(cds))/sizeFactors(cds))
   rawVars <- rawVarFunc( cds, "_pooled", TRUE )( baseMeans )
   rawScv <- adjustScvForBias( rawVars/baseMeans^2, attr( rawVars, "size" ) )

   nbinomGLMsForMatrix( counts(cds), sizeFactors(cds), rawScv, 
      modelFormula, pData(cds), glmControl=glmControl )
}

nbinomGLMTest <- function( resFull, resReduced )
   1 - pchisq( resReduced$deviance - resFull$deviance, 
   attr( resReduced, "df.residual" ) - attr( resFull, "df.residual" ) )

getRawScvForSamplePair <- function( cds, sample1, sample2 ) {
   if( any( is.na( sizeFactors(cds) ) ) )
      stop( "Call 'estimateSizeFactors' first." )
   data <- counts(cds)[ ,c( sample1, sample2 ) ]
   sf <- sizeFactors(cds)[ c( sample1, sample2 ) ]
   vf <- estimateVarianceFunctionForMatrix( data, sf )
   bm <- t( t(data) / sf )
   mean( adjustScvForBias( vf( bm ) / bm^2, 2 ), na.rm=TRUE ) }

getRawScvDistanceMatrix <- function( cds ) {
   res <- matrix( nrow=ncol(cds), ncol=ncol(cds) )
   for( i in 1:ncol(cds) )
      res[ i, i ] <- 0
   for( i in 1:(ncol(cds)-1) )
      for( j in (i+1):ncol(cds) ) {
         res[ i, j ] <- getRawScvForSamplePair( cds, i, j )
         res[ j, i ] <- res[ i, j ] }
   colnames( res ) <- colnames( counts(cds) )
   rownames( res ) <- colnames( counts(cds) )
   res }
