# install.packages("matrixStats")

#################################################################################
# Some functions
load("/Users/goshng/Library/R/2.12/library/DESeq/extra/scvBiasCorrectionFits.rda")
adjustScvForBias <- function( scv, nsamples ) {
   stopifnot( nsamples > 1 )
   if( nsamples - 1 > length( scvBiasCorrectionFits ) )
      scv
   else
      pmax( safepredict( scvBiasCorrectionFits[[ nsamples-1 ]], scv ), 1e-8 * scv )
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
   # totals <- calc_pvals ( as.integer(kA+kB), pobs, muA, vA, muB, vB, eps )
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
#
#################################################################################

#################################################################################
# Main procedure
library(locfit)
library(matrixStats)
countsFile <- "/Users/goshng/Documents/Projects/rnaseq-analysis/output/omz/1/bwa/count.txt"
conds <- c("OMZ175", "OMZHKRR")
condA <- "OMZ175"
condB <- "OMZHKRR"
eps <- 1e-8

# Read a read counts file
countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- countsTable$gene
countsTable <- countsTable[,-1]
counts<- as.matrix(countsTable)
conditions <- factor( conds )

# Estimate the size factors
loggeomeans <- rowMeans( log(counts) ) 
sffunc <- function(cnts) exp( median( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) 
sizeFactors <- apply( counts, 2, sffunc)

# Estimate variance functions for method of blind
bmv <- data.frame(
         baseMean = rowMeans( t( t(counts) / sizeFactors ) ),
         baseVar = rowVars( t( t(counts) / sizeFactors ) ) 
       )
variances <- bmv$baseVar[ bmv$baseMean > 0 ]
means <-  bmv$baseMean[ bmv$baseMean > 0 ]
fit <- locfit(variances ~ lp(log(means)),family="gamma")
rm( means )
rm( variances )
xim <- sum( 1/sizeFactors ) / length( sizeFactors )
      
rawVarFunc <- function( q ) {
  x <- log(q)
  res <- rep.int( NA_real_, length(x) )
  res[ is.finite(x) ] <- predict( fit, x[is.finite(x)] )

  ans <- pmax( res - xim * q, 1e-8 * q )
  attr( ans, "size" ) <- length( sizeFactors )
  ans 
}

# Test of the null hypotheses
rvfA <- rawVarFunc
rvfB <- rawVarFunc
colA <- conditions==condA
colB <- conditions==condB
bmvAB <- data.frame(
           baseMean = rowMeans( t( t(counts[,colA|colB]) / sizeFactors[colA|colB] ) ),
           baseVar = rowVars( t( t(counts[,colA|colB]) / sizeFactors[colA|colB] ) ) 
         )
rawScvA <- rvfA( bmvAB$baseMean ) / bmvAB$baseMean^2
rawScvB <- rvfB( bmvAB$baseMean ) / bmvAB$baseMean^2
rawScvA <- adjustScvForBias( rawScvA, attr( rawScvA, "size" ) )
rawScvB <- adjustScvForBias( rawScvB, attr( rawScvB, "size" ) )
countsA <- counts[,colA]
countsB <- counts[,colB]
sizeFactorsA <- sizeFactors[colA]
sizeFactorsB <- sizeFactors[colB]
kAs <- rowSums( cbind(countsA) )
kBs <- rowSums( cbind(countsB) )
baseMeans <- rowMeans( cbind(t( t( countsA ) / sizeFactorsA ),
                             t( t( countsB ) / sizeFactorsB ) ) 
                     )      
muAs <- baseMeans * sum( sizeFactorsA )
muBs <- baseMeans * sum( sizeFactorsB )
fullVarA <- pmax( muAs + rawScvA * baseMeans^2 * sum(sizeFactorsA^2), muAs * (1+1e-8) )
fullVarB <- pmax( muBs + rawScvB * baseMeans^2 * sum(sizeFactorsB^2), muBs * (1+1e-8) )

nbinomTestFuncCall <- function (i)
{
  nbinomTestForMatricesRaw( kAs[i], kBs[i], muAs[i], fullVarA[i], muBs[i], fullVarB[i], eps )
}
pval <- sapply(1:nrow(cbind(countsA)), 
               function (i)
               {
                 nbinomTestForMatricesRaw( kAs[i], kBs[i], muAs[i], fullVarA[i], muBs[i], fullVarB[i], eps )
               }
              )
bmvA <- data.frame(
          baseMean = rowMeans( t( t(counts[,colA]) / sizeFactors[colA] ) ),
          baseVar = rowVars( t( t(counts[,colA]) / sizeFactors[colA] ) ) 
        )
bmvB <- data.frame(
          baseMean = rowMeans( t( t(counts[,colB]) / sizeFactors[colB] ) ),
          baseVar = rowVars( t( t(counts[,colB]) / sizeFactors[colB] ) ) 
        )
result <- data.frame( 
         id    = rownames( counts ),
         baseMean  = bmvAB$baseMean,
         baseMeanA = bmvA$baseMean,
         baseMeanB = bmvB$baseMean,
         foldChange = bmvB$baseMean / bmvA$baseMean,
         log2FoldChange = log2( bmvB$baseMean / bmvA$baseMean ), 
         pval = pval,
         padj = p.adjust( pval, method="BH" ), 
         resVarA = bmvA$baseVar / ( bmvA$baseMean * sum( 1/sizeFactors[colA] ) / length(condA) +
            rvfA( bmvAB$baseMean ) ),
         resVarB = bmvB$baseVar / ( bmvB$baseMean * sum( 1/sizeFactors[colB] ) / length(condB) +
            rvfB( bmvAB$baseMean ) ),
         stringsAsFactors = FALSE )


