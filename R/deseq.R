# install.packages("matrixStats")
library(locfit)
library(matrixStats)
library(Biobase)

####################################################################################
# class_and_slots.R: Class, get, and set functions
# Class definition
setClass( "CountDataSet", 
   contains = "eSet",
   representation = representation( 
      rawVarFuncs = "environment",
      rawVarFuncTable = "character",
      multivariateConditions = "logical" ),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), CountDataSet = "1.1.0" ) ) )
)

# Gets counts
counts <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   assayData(cds)[["counts"]]
}   
   
# Gets sizeFactors
sizeFactors <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   sf <- pData(cds)$`sizeFactor`
   names( sf ) <- colnames( counts(cds) )
   sf
}   

# Sets sizeFactors 
`sizeFactors<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   pData(cds)$`sizeFactor` <- value
   validObject( cds )
   cds
}   

# Gets conditions
conditions <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( cds@multivariateConditions )
      stop( "The 'conditions' accessor is only for simple single-factor conditions, but your have specified multivariate conditions. Access them via 'pData'." )
   conds <- pData(cds)$`condition`
   names( conds ) <- colnames( counts(cds) )
   conds
}   

# Sets conditions
`conditions<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( cds@multivariateConditions )
      stop( "The 'conditions' accessor is only for simple single-factor conditions, but your have specified multivariate conditions. Access them via 'pData'." )
   pData(cds)$`condition` <- factor( value )
   validObject( cds )
   cds
}   

# Gets rawVarFunc
rawVarFunc <- function( cds, condOrName=NULL, byName=FALSE ) {
   stopifnot( is( cds, "CountDataSet" ) )
   ensureHasVarFuncs( cds )   
   if( is.null( condOrName ) ) {
      if( length(cds@rawVarFuncs) == 1 )
         return( cds@rawVarFuncs[[ ls(cds@rawVarFuncs)[[1]] ]] )
      else
         stop( "There is more than one variance function. 'condOrName' may not be omitted." )
   }   
   if( byName ) {      
      res <- cds@rawVarFuncs[[ as.character(condOrName) ]]
      if( is.null(res) )
         stop( sprintf( "No raw variance function found with name '%s'.", condOrName ) )
   } else {      
      res <- cds@rawVarFuncs[[ cds@rawVarFuncTable[ as.character(condOrName) ] ]]
      if( is.null(res) )
         stop( sprintf( "No raw variance function found for condition '%s'.", condOrName ) )
   }
   res
}

# Gets rawVarFuncTable
rawVarFuncTable <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   cds@rawVarFuncTable
}   

# Sets rawVarFuncTable
`rawVarFuncTable<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( is.null( names(value) ) )
      names( value ) <- names( rawVarFuncTable(cds) )
   cds@rawVarFuncTable <- value
   validObject( cds )   
   cds
}   

# Checks if VarFuncs are defined
ensureHasVarFuncs <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( length(ls(cds@rawVarFuncs)) == 0 )
      stop( "CountDataSet object does not contain any variance functions. Call 'estimateVarianceFunctions' first." )
   TRUE
}   

varAdjFactors <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   stop( "This function has been removed. Do not use it. See help page." )
}

`varAdjFactors<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   stop( "This function has been removed. Do not use it. See help page." )
}

####################################################################################
# methods.R: 
# Method: sizeFactors
estimateSizeFactors <- function( cds, locfunc=median )
{
   stopifnot( is( cds, "CountDataSet" ) )
   sizeFactors(cds) <- estimateSizeFactorsForMatrix( counts(cds), locfunc )
   cds
}

# Method: variance functions
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
	       rawVarFuncTable(cds) <- a 
      } 
   }
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
         rawVarFuncTable(cds) <- a 
      } 
   }
        
   validObject( cds )
   cds
}

# Method: nbinomTest
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

####################################################################################
# core.R
estimateSizeFactorsForMatrix <- function( counts, locfunc = median )
{
   loggeomeans <- rowMeans( log(counts) ) 
   apply( counts, 2, function(cnts) exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
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

prepareScvBiasCorrectionFits <- function( maxnrepl=15, mu=100000, ngenes=10000,
      true_raw_scv = c( seq( 0, 2, length.out=100 )[-1], seq( 2, 10, length.out=20 )[-1] ) )
   lapply( 2:maxnrepl, function( m ) {
      est_raw_scv <- sapply( true_raw_scv, function( alpha ) {
         k <- matrix( rnbinom( ngenes*m, mu=mu, size=1/alpha ), ncol=m )
         k <- k[ rowSums(k)>0, ]
         mean( rowVars(k) / rowMeans(k)^2 ) } )
      locfit( true_raw_scv ~ lp( est_raw_scv, nn=.2 ) ) } )

load("/Users/goshng/Library/R/2.12/library/DESeq/extra/scvBiasCorrectionFits.rda")

adjustScvForBias <- function( scv, nsamples ) {
   stopifnot( nsamples > 1 )
   if( nsamples - 1 > length( scvBiasCorrectionFits ) )
      scv
   else
      pmax( safepredict( scvBiasCorrectionFits[[ nsamples-1 ]], scv ), 1e-8 * scv )
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

   # totals <- calc_pvals ( as.integer(kA+kB), pobs, muA, vA, muB, vB, eps )
   totals <- .Call( "calc_pvals", as.integer(kA+kB), pobs, muA, vA, muB, vB, eps )
   min( unname( totals[2] / totals[1] ), 1 )
   # The 'min' is to avoid p values slightly exceeding 1 due to
   # approximation errors
}

add_from_both_sides <- function( left, right, kS, pobs, muA, vA, muB, vB, eps)
{
   sizeA = muA*muA / (vA-muA)
   probA = muA/vA
   sizeB = muB*muB / (vB-muB)
   probB = muB/vB

   kl = left
   kr = right
   lval = dnbinom( kl, sizeA, probA ) * dnbinom( kS-kl, sizeB, probB )
   rval = dnbinom( kr, sizeA, probA ) * dnbinom( kS-kr, sizeB, probB )
   prevlval = lval
   prevrval = rval
   total = lval + rval
   esttotalperlength = total/2
   obstotal = 0

   step = 1
   steps = 0
   do_left = TRUE
   if( lval <= pobs )
   {
      obstotal = obstotal + lval
   }
   if( rval <= pobs )
   {
      obstotal = obstotal + rval
   }
   while( kl < kr ) {
      #if (steps %% 1000 == 0)
      #{
        #print (steps)
      #}
      steps <- steps + 1
      if( abs( ( prevrval - rval ) ) / prevrval > .01 )
      {
         do_left = TRUE
      }
      else 
      if( abs( ( prevlval - lval ) ) / prevlval > .01 )
      {
         do_left = FALSE
      }
      else
      {
         do_left = lval > rval
      }
      if( do_left ) {
         prevlval = lval
         if( kl + step > kr )
         {
            step = kr - kl
         }
         kl = kl + step
         lval = dnbinom( kl, sizeA, probA ) * dnbinom( kS-kl, sizeB, probB )
         if( step == 1 )
         {
            total = total + lval
         }
         else
         {
            total = total + min( lval, prevlval ) * step
         }
         if( lval <= pobs ) {
            if( step == 1 )
            {
               obstotal = obstotal + lval
            }
            else {       
               if( prevlval <= pobs )
               {
                  obstotal = obstotal + max( lval, prevlval ) * step
               }
               else
               {
                  obstotal = obstotal + max( lval, prevlval ) * step * abs( (pobs-lval) / (prevlval-lval) )
               }
            }
         }       
         if( abs( prevlval - lval ) / prevlval < eps )
         {
            step = max( step + 1, round(step * 1.5) )
         }
      } else {
         prevrval = rval
         if( kr - step < kl )
         {
            step = kr - kl
         }
         kr = kr - step
         rval = dnbinom( kr, sizeA, probA ) * dnbinom( kS-kr, sizeB, probB )
         if( step == 1 )
         {
            total = total + rval
         }
         else
         {
            total = total + min( rval, prevrval ) * step
         }
         if( rval <= pobs ) {
            if( step == 1 )
            {
               obstotal = obstotal + rval
            }
            else {       
               if( prevrval <= pobs )
               {
                  obstotal = obstotal + max( rval, prevrval ) * step
               }
               else
               {
                  obstotal = obstotal + max ( rval, prevrval ) * step * abs( (pobs-rval) / (prevrval-rval) )
               }
            }
         }       
         if( abs( ( prevrval - rval ) ) / prevrval < eps )
         {
            step = max( step + 1, round(step * 1.5) )
         }
      }
   }
   
   list(total=total, obstotal=obstotal)
}   
calc_pvals <- function ( kS, pobs, muA, vA, muB, vB, eps ) 
{
   kexp <- round(kS * muA / ( muA + muB ))
   vl <- add_from_both_sides( 0,      kexp, kS, pobs, muA, vA, muB, vB, eps)
   vr <- add_from_both_sides( kexp+1, kS,   kS, pobs, muA, vA, muB, vB, eps)

   c(vl$total + vr$total, vl$obstotal + vr$obstotal)
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
   #if (i %% 100 == 0)
     #print (i)
      nbinomTestForMatricesRaw( kAs[i], kBs[i], muAs[i], fullVarA[i], muBs[i], fullVarB[i], eps )
   } )
}

####################################################################################
# Main function
# data file, column names from the second one, two conditions to compare
# Read a table of counts: the first column is gene. Each column represents
# a sample, and the column names does not mean much. You need to specify 
# which columns to compare.
# countsFile <- "/Users/goshng/Documents/Projects/rnaseq-analysis/output/omz/1/bwa/count.txt"
# conds <- c("OMZ175", "OMZHKRR")
# e.g., r1 <- de(countsFile, conds, "OMZ175", "OMZHKRR")
de <- function(countsFile, conditions, condA, condB, method)
{
  countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
  rownames(countsTable) <- countsTable$gene
  countsTable <- countsTable[,-1]

  sizeFactors=NULL
  phenoData = NULL
  featureData = NULL
  countData <- countsTable 
  countData <- as.matrix( countData )
  if( any( round( countData ) != countData ) )
    stop( "The countData is not integer." )
  mode( countData ) <- "integer"

  if ( is.null( sizeFactors ) )
    sizeFactors <- rep( NA_real_, ncol(countData) )
  if ( is.null( phenoData ) )
    phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
  if ( is.null( featureData ) ) 
    featureData <- annotatedDataFrameFrom( countData, byrow=TRUE )

  phenoData$`sizeFactor` <- sizeFactors
  varMetadata( phenoData )[ "sizeFactor", "labelDescription" ] <-
    "size factor (relative estimate of sequencing depth)"
  if( is( conditions, "matrix" ) )
    conditions <- as.data.frame( conditions )

  if( is( conditions, "data.frame" ) || is( conditions, "AnnotatedDataFrame" ) ) {
     stopifnot( nrow( conditions ) == ncol( countData ) )
     conditions <- as( conditions, "AnnotatedDataFrame" )
     dimLabels( conditions ) <- dimLabels( phenoData )
     rownames( pData(conditions) ) <- rownames( pData(phenoData) )
        # TODO: What if the rownames were set?
     phenoData <- combine( phenoData, conditions )
     multivariateConditions <- TRUE
     rvft <- c( `_all` = NA_character_ )
  } else {
     conditions <- factor( conditions )
     stopifnot( length( conditions ) == ncol( countData ) )
     phenoData$`condition` <- factor( conditions )
     varMetadata( phenoData )[ "condition", "labelDescription" ] <-
       "experimental condition, treatment or phenotype"
     multivariateConditions <- FALSE
     rvft <- rep( NA_character_, length(levels(conditions)) )
  }
  cds <- new( "CountDataSet",
    assayData = assayDataNew( "environment", counts=countData ),
    phenoData = phenoData, 
    featureData = featureData,
    multivariateConditions = multivariateConditions,
    rawVarFuncs = new.env( hash=TRUE ),
    rawVarFuncTable = rvft )
  cds <- estimateSizeFactors(cds)
  cds <- estimateVarianceFunctions(cds,method=method)
  res <- nbinomTest(cds, condA, condB)
}

countsFile <- "/Users/goshng/Documents/Projects/rnaseq-analysis/output/omz/1/bwa/count.txt"
conds <- c("OMZ175", "OMZHKRR")
method <- "blind"
r1 <- de(countsFile, conds, "OMZ175", "OMZHKRR", method)


# R: apply's
# http://nsaunders.wordpress.com/2010/08/20/a-brief-introduction-to-apply-in-r/

# for rowVars
# install.packages("matrixStats")
# library(matrixStats)
