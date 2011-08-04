estimateSizeFactorsForMatrix <- function( counts, locfunc = median )
{
   loggeomeans <- rowMeans( log(counts) ) 
   apply( counts, 2, function(cnts) 
      exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
}


getBaseMeansAndVariances <- function( counts, sizeFactors ) {

   # Devides the counts by sizeFactors and calculates the estimates for
   # base means and variances for each gene.
   
   data.frame(
      baseMean = rowMeans( t( t(counts) / sizeFactors ) ),
      baseVar = rowVars( t( t(counts) / sizeFactors ) ) )
}   

estimateVarianceFunctionForMatrix <- function( counts, sizeFactors, 
         locfit_extra_args=list(), lp_extra_args=list() ) {

   stopifnot( ncol( counts ) == length( sizeFactors ) )
   bmv <- getBaseMeansAndVariances( counts, sizeFactors ) 
   estimateVarianceFunctionFromBaseMeansAndVariances( bmv$baseMean,
      bmv$baseVar, sizeFactors, locfit_extra_args, lp_extra_args )
}      
   
modelMatrixToConditionFactor <- function( modelMatrix ) {

   mmconds <- 1:nrow(modelMatrix)
   for( i in 2:nrow(modelMatrix) )
      for( j in 1:(i-1) )
         if( all( modelMatrix[i,] == modelMatrix[j,] ) ) {
            mmconds[i] = mmconds[j]
            break }
   factor( as.integer( factor( mmconds ) ) )
}


getBaseMeansAndPooledVariances <- function( counts, sizeFactors, conditions ) {

   basecounts <- t( t(counts) / sizeFactors )
   replicated_sample <- conditions %in% names(which(table(conditions)>1))
   df <- sum(replicated_sample) - length( unique( conditions[ replicated_sample ] ) ) 

   data.frame(
      baseMean = rowMeans( basecounts ),
      baseVar =
	 rowSums( 
	    sapply( 
               tapply( 
        	  ( 1:ncol(counts) )[ replicated_sample ], 
        	  factor( conditions[ replicated_sample ] ), 
        	  function(cols) 
        	     rowSums( ( basecounts[,cols] - rowMeans(basecounts[,cols]) )^2 ) ), 
               identity ) ) / df )
}
   
estimatePooledVarianceFunctionForMatrix <- function( counts, sizeFactors, 
      conditions, locfit_extra_args=list(), lp_extra_args=list() ) {
      
   bmv <- getBaseMeansAndPooledVariances( counts, sizeFactors, conditions ) 
   estimateVarianceFunctionFromBaseMeansAndVariances( bmv$baseMean,
      bmv$baseVar, sizeFactors, locfit_extra_args, lp_extra_args )
}      
   
   
estimateVarianceFunctionFromBaseMeansAndVariances <- function( means, 
   variances, sizeFactors, locfit_extra_args=list(), lp_extra_args=list() ) {
   
   variances <- variances[ means > 0 ]
   means <- means[ means > 0 ]
   
   fit <- do.call( "locfit", c( 
      list( 
         variances ~ do.call( "lp", c( list( log(means) ), lp_extra_args ) ),
         family = "gamma" ), 
      locfit_extra_args ) )
   
   rm( means )
   rm( variances )
   xim <- sum( 1/sizeFactors ) / length( sizeFactors )
      
   function( q ) {
      ans <- pmax( safepredict( fit, log(q) ) - xim * q, 1e-8 * q )
      attr( ans, "size" ) <- length( sizeFactors )
      ans }
   # Note: The 'pmax' construct above serves to limit the overdispersion to a minimum
   # of 10^-8, which should be indistinguishable from 0 but ensures numerical stability.
}   
   
safepredict <- function( fit, x )
{
   # A wrapper around predict to avoid the issue that predict.locfit cannot
   # propagate NAs and NaNs properly.

   res <- rep.int( NA_real_, length(x) )
   res[ is.finite(x) ] <- predict( fit, x[is.finite(x)] )
   res   
}

nbinomTestForMatricesRaw <- function( kA, kB, muA, vA, muB, vB, eps=0 )
{
   # Let kA and kB be two count observations from two random variables, for
   # which the null hypothesis assumes negative binomial distributions with
   # means muA, muB and vA and vB. Calculate the probability that kA and kB
   # or as or more extreme counts are observed, conditioned to the sum of the 
   # counts being kA+kB. "As or more extreme" means having conditional probability 
   # at most 
   #                     fNB( kA, muA, vA ) fNB( kA, muA, vA )
   #       -------------------------------------------------------------------
   #        sum of fNB( k, muA, vA ) fNB( kA+kB-k, muA, vA ) for k=0,..,kA+kB
   #
   # 'eps' is a roughly followed guidance on the required presision
   
   if( !all( is.finite( c( kA, kB, muA, vA, muB, vB, eps ) ) ) )
      return( NA )

   pobs <- dnbinom( kA, prob = muA/vA, size = muA^2/(vA-muA) ) * 
           dnbinom( kB, prob = muB/vB, size = muB^2/(vB-muB) )
           
   stopifnot( is.finite( pobs ) )
   
   pobs <- pobs * ( 1 + 1e-7 )
   # This is to avoid rounding errors in checking for p <= pobs

   totals <- .Call( "calc_pvals", as.integer(kA+kB), pobs, muA, vA, muB, vB, eps )
   min( unname( totals[2] / totals[1] ), 1 )
   # The 'min' is to avoid p values slightly exceeding 1 due to
   # approximation errors
}


nbinomTestForMatrices <- function( countsA, countsB, sizeFactorsA, sizeFactorsB, 
   rawScvA, rawScvB, eps=1e-4 )
{
   kAs <- rowSums( cbind(countsA) )
   kBs <- rowSums( cbind(countsB) )
   
   baseMeans <- rowMeans( cbind(      
      t( t( countsA ) / sizeFactorsA ),
      t( t( countsB ) / sizeFactorsB ) ) )      
   muAs <- baseMeans * sum( sizeFactorsA )
   muBs <- baseMeans * sum( sizeFactorsB )

   fullVarA <- pmax( muAs + rawScvA * baseMeans^2 * sum(sizeFactorsA^2), muAs * (1+1e-8) )
   fullVarB <- pmax( muBs + rawScvB * baseMeans^2 * sum(sizeFactorsB^2), muBs * (1+1e-8) )
   
   sapply( 1:nrow(cbind(countsA)), function(i) {
      nbinomTestForMatricesRaw( kAs[i], kBs[i], muAs[i], fullVarA[i], muBs[i], fullVarB[i], eps )
   } )
}


varianceFitDiagnosticsForMatrix <- function( counts, sizeFactors, rawVarFunc, poolingConditions=NULL )
{
   if( is.null( poolingConditions ) )
      res <- getBaseMeansAndVariances( counts, sizeFactors )
   else
      res <- getBaseMeansAndPooledVariances( counts, sizeFactors, poolingConditions )
   res$fittedRawVar <- rawVarFunc( res$baseMean )
   res$fittedBaseVar <- res$fittedRawVar + 
      res$baseMean * sum( 1/sizeFactors ) / length( sizeFactors )
   df <- ncol( cbind(counts) ) - 1
   res$pchisq <- pchisq( df * res$baseVar / res$fittedBaseVar, df = df )
   res
}


multiecdfWithoutLegend <- function( x, ... )
{
   if( all( package_version( packageDescription( "geneplotter", fields="Version" ) ) 
          >= package_version( "1.21.2" ) ) )
      multiecdf( x, ..., legend=NULL )
   else
      multiecdf( x, ... ) 
}

residualsEcdfPlotFromDiagnostics <- function( fitdiag, ncuts=7, 
      plotTitle="Residuals ECDF plot" )
{
   ok <- !is.na(fitdiag$pchisq)
   cols <- colorRampPalette( c("red","blue") )( ncuts )
   cuts <- factor( cut( rank(fitdiag$baseMean[ok]), ncuts ) )

   multiecdfWithoutLegend( 
      fitdiag$pchisq[ok] ~ cuts,
      col = cols,
      xlab = "chi-squared probability of residual",
      ylab = "ECDF",
      main = plotTitle
   )
      
   segments( 0, 0, 1, 1, col="darkgreen" )
   legend( 0, 1, 
      c( sprintf( "%.1e .. %.1e", 
            tapply( fitdiag$baseMean[ok], cuts, min ), 
            tapply( fitdiag$baseMean[ok], cuts, max ) ), 
         "expected" ),
      col = c( cols, "darkgreen" ), lty="solid" )
}  

# Note: The following function is never called; it is here only for
# documentation purposes, as it has been used to produce the data object
# scvBiasCorrectionFits, which is stored in the file
# inst/scvBiasCorrectionFits.rda, gets loadewd by the line after this
# function and is used by the function adjustScvForBias

prepareScvBiasCorrectionFits <- function( maxnrepl=15, mu=100000, ngenes=10000,
      true_raw_scv = c( seq( 0, 2, length.out=100 )[-1], seq( 2, 10, length.out=20 )[-1] ) )
   lapply( 2:maxnrepl, function( m ) {
      est_raw_scv <- sapply( true_raw_scv, function( alpha ) {
         k <- matrix( rnbinom( ngenes*m, mu=mu, size=1/alpha ), ncol=m )
         k <- k[ rowSums(k)>0, ]
         mean( rowVars(k) / rowMeans(k)^2 ) } )
      locfit( true_raw_scv ~ lp( est_raw_scv, nn=.2 ) ) } )

load( system.file ( "extra/scvBiasCorrectionFits.rda", package="DESeq" ) )

adjustScvForBias <- function( scv, nsamples ) {
   stopifnot( nsamples > 1 )
   if( nsamples - 1 > length( scvBiasCorrectionFits ) )
      scv
   else
      pmax( safepredict( scvBiasCorrectionFits[[ nsamples-1 ]], scv ), 1e-8 * scv )
}      

nbkd.sf <- function( r, sf ) {
   fam <- list(
     
      family = sprintf( "nbkd,r=%g", r ),
      link = "log_sf",
      linkfun = function( mu ) log( mu/sf ),
      linkinv = function (eta) pmax(sf*exp(eta), .Machine$double.eps),
      mu.eta = function (eta) pmax(sf*exp(eta), .Machine$double.eps),
      variance = function(mu) mu + mu^2 / r,

      dev.resids = function( y, mu, wt )
          2 * wt * ( ifelse( y > 0, y * log( y / mu ), 0 ) + 
            (r + y) * log( (r+mu)/(r+y) ) ), 
      
      aic = function (y, n, mu, wt, dev) NA,   # what for?
      initialize = expression( {
         n <- rep.int(1, nobs)         # What is n?
         mustart <- y + 0.1
      } ),
      valid.mu <- function(mu) all( mu > 0 ),
      valid.eta <- function(eta) TRUE,
      simulate <- NA
   )
   
   class(fam) <- "family"
   fam }        


nbinomGLMsForMatrix <- function( counts, sizeFactors, rawScv, modelFormula, 
   modelFrame, quiet=FALSE, reportLog2=TRUE, glmControl=list() ) 
{
   stopifnot( length(sizeFactors) == ncol(counts) )
   stopifnot( length(rawScv) == nrow(counts) )
   stopifnot( nrow(modelFrame) == ncol(counts) )

   stopifnot( is( modelFormula, "formula" ) )  
   if( as.character( modelFormula[[1]] ) != "~" )  
      stop( "Formula does not have a '~' as top-level operator." )
   if( as.character( modelFormula[[2]] ) != "count" )  
      stop( "Left-hand side of model formula must be 'count'." )
   
   goodRows <- !is.na( rawScv ) & rowSums(counts) >= 0 
   
   res <- 
   t( sapply( which(goodRows), function(i) {
      if( !quiet & i %% 1000 == 0 )
         cat( '.' ) 
      nbfam <- nbkd.sf( 1 / rawScv[i], sizeFactors )      
      fit <- glm( modelFormula, cbind( count=counts[i,], modelFrame ), 
         family=nbfam, control = glmControl )
      c( 
         coefficients(fit), 
         deviance = deviance(fit), 
         df.residual = fit$df.residual,
         converged = fit$converged ) } ) )
      
   if( !quiet )
      cat( "\n" ) 

   df.residual <- res[ 1, "df.residual" ]
   stopifnot( all( res[ , "df.residual" ] == df.residual ) )

   # Put in the NAs
   res2 <- data.frame(
      row.names = row.names( counts ),
      apply( res, 2, function( col ) {
         a <- rep( NA_real_, nrow(counts) )
         a[goodRows] <- col
         a } ) )
   colnames(res2) <- colnames(res)
      
   if( reportLog2 )
      res2[ , 1:(ncol(res2)-3) ] <- res2[ , 1:(ncol(res2)-3) ] / log(2)

   res2$converged <- as.logical( res2$converged )

   res2 <- res2[ , colnames(res) != "df.residual" ]     
   attr( res2, "df.residual" ) <- df.residual
   res2
}
      
