
# Step 1: Read a table
newCountDataSet <- function()
{
   sizeFactors=NULL
   phenoData = NULL
   featureData = NULL
   conditions <- conds 
   countData <- countsTable 
   countData <- as.matrix( countData )
   if( any( round( countData ) != countData ) )
      stop( "The countData is not integer." )
   mode( countData ) <- "integer"

   if( is.null( sizeFactors ) )
      sizeFactors <- rep( NA_real_, ncol(countData) )
   if( is.null( phenoData ) )
      phenoData <- annotatedDataFrameFrom( countData, byrow=FALSE )
   if( is.null( featureData ) ) 
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
   cds
}

# size factor
counts <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   assayData(cds)[["counts"]]
}   
 
estimateSizeFactorsForMatrix <- function( counts, locfunc = median )
{
   loggeomeans <- rowMeans( log(counts) ) 
   apply( counts, 2, function(cnts) exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
}

estimateSizeFactors <- function( cds, locfunc=median )
{
   stopifnot( is( cds, "CountDataSet" ) )
   sizeFactors(cds) <- estimateSizeFactorsForMatrix( counts(cds), locfunc )
   cds
}

# Step 2: compute the size factor
estimateSizeFactorsSangChulChoi <- function()
{
   loggeomeans <- rowMeans( log(counts(cds)) ) 
   sizeFactors(cds) <- apply( counts(cds), 2, function(cnts) exp( median( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
}

########################################################################
# Step 3:
getBaseMeansAndVariances <- function( counts, sizeFactors ) {

   # Devides the counts by sizeFactors and calculates the estimates for
   # base means and variances for each gene.
   
   data.frame(
      baseMean = rowMeans( t( t(counts) / sizeFactors ) ),
      baseVar = rowVars( t( t(counts) / sizeFactors ) ) )
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
   
estimateVarianceFunctionsSangChulChoi <- function()
{
   stopifnot( is( cds, "CountDataSet" ) )   
   if( any( is.na( sizeFactors(cds) ) ) )
      stop( "NAs found in size factors. Have you called already 'estimateSizeFactors'?" )
   method="blind"
   c1 <- counts(cds)
   s1 <- sizeFactors(cds)
   baseMean = rowMeans( t( t(c1) / s1 ) )
   baseVar = rowVars( t( t(c1) / s1 ) ) 
   bmv <- data.frame(baseMean, baseVar) 

   # From estimateVarianceFunctionFromBaseMeansAndVariances <- function
   baseVar <- baseVar[ baseMean > 0 ]
   baseMean <- baseMean [ baseMean > 0 ]

   # Test of local regression 
   # data(ethanol, package="locfit")
   # fit <- locfit(NOx ~ E, data=ethanol)
   # fit <- locfit(NOx~lp(E,nn=0.5),data=ethanol)
   # fit <- locfit(NOx~lp(E,C,scale=TRUE),data=ethanol)
   # plot(fit, get.data=TRUE)

   # Return a function
}

# Step 4: nbiomtest
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

   # Functions?
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

commentfunc2 <- function ()
{
library(DESeq)
countsFile <- "/Users/goshng/Documents/Projects/rnaseq/output/omz/1/bwa/count.txt"
countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- countsTable$gene
countsTable <- countsTable[,-1]
conds <- c("OMZ175", "OMZHKRR")
cds <- newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds,method="blind")
res <- nbinomTest(cds, "OMZ175", "OMZHKRR")
}

plotDE <- function( res )
{
   plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.5, 
        col = ifelse( res$padj < .1, "red", "black" ) )
}

commentfunc <- function ()
{
cat ("SCV file: count-de.ps\n")
postscript ("/Users/goshng/Documents/Projects/rnaseq/output/omz/1/bwa/count-de.ps",  width=10, height=10, horizontal = FALSE, onefile = FALSE,
paper = "special")
plotDE(res)
y <- dev.off()

cat ("SCV file: count-scv.ps\n")
postscript ("/Users/goshng/Documents/Projects/rnaseq/output/omz/1/bwa/count-scv.ps",  width=10, height=10, horizontal = FALSE, onefile = FALSE,
paper = "special")
scvPlot(cds)
y <- dev.off()

resSig <- res[res$padj < .1,]
y <- length(countsTable[,1])
cat ("Total number of genes is ", y, "\n",sep="")
y <- length(resSig[,1])
cat ("The number of differentially expressed genes is ", y, "\n",sep="")
cat ("\nList of differentially expressed genes\n")
print(resSig[order(resSig$pval),],with=1000)
cat ("\nThe same table as above but in order of down-regulated first of the differentially expressed genes\n")
print(resSig[order(resSig$foldChange,-resSig$baseMean),])
cat ("\nThe same table as above but in order of up-regulated first of the differentially expressed genes\n")
print(resSig[order(-resSig$foldChange,-resSig$baseMean),])
}
