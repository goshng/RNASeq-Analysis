setClass( "CountDataSet", 
   contains = "eSet",
   representation = representation( 
      rawVarFuncs = "environment",
      rawVarFuncTable = "character",
      multivariateConditions = "logical" ),
   prototype = prototype( new( "VersionedBiobase",
      versions = c( classVersion("eSet"), CountDataSet = "1.1.0" ) ) )
)

newCountDataSet <- function( countData, conditions, sizeFactors=NULL,
      phenoData = NULL, featureData = NULL )
{
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

setValidity( "CountDataSet", function( object ) {
   if( length( object@multivariateConditions ) != 1 )
      return( "multivariateConditions is not scalar." )
   if( ! "sizeFactor"  %in% names(pData(object)) )
      return( "phenoData does not contain a 'sizeFactor' columns.")
   if( ! is( pData(object)$`sizeFactor`, "numeric" ) )
      return( "The 'sizeFactor' column in phenoData is not numeric." )
   if( ! object@multivariateConditions ) {
      if( ! "condition"  %in% names(pData(object)) )
         return( "phenoData does not contain a 'condition' columns." )
      if( ! is( pData(object)$`condition`, "factor" ) )
         return( "The 'condition' column in phenoData is not a factor." )
      if( length(object@rawVarFuncTable) != length(levels(conditions(object))) )
         return( "The rawVarFuncTable does not contain one element per condition." )
      if( any( names(object@rawVarFuncTable) != levels(conditions(object))) )
         return( "The names of the character vector 'rawVarFuncTable' are not identical to the levels of the factor 'conditions'." ) }
   else {
      if( length(object@rawVarFuncTable) != 1 )
         return( "The rawVarFuncTable must be of length 1 if mutivariate conditions are used." )
      if( names(object@rawVarFuncTable) != "_all" )
         return( "The rawVarFuncTable must have one element named '_all' if mutivariate conditions are used." ) }
   for( funcName in object@rawVarFuncTable )
      if( ! is.na( funcName ) )
         if( is.null( object@rawVarFuncs[[funcName]] ) )
            return( sprintf( "rawVarFuncTable mentions a rawVarFunction '%s' which is missing.",
               funcName ) )
   for( rvfName in ls(object@rawVarFuncs) ) {
      if( !is( object@rawVarFuncs[[rvfName]], "function" ) )
         return( sprintf( "rawVarFuncs contains a value, called '%s', which is not a function.", rvfName ) )
      if( length( formals( object@rawVarFuncs[[rvfName]] ) ) != 1 ) 
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does have the right argument list for a raw variance function.", rvfName ) )
      testres <- object@rawVarFuncs[[rvfName]]( c( 1.5, 3.5, 7.3 ) ) 
      if( ! is( testres, "numeric" ) )
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does not return a numeric result.", rvfName ) )
      if( length(testres) != 3 )
         return( sprintf( "rawVarFuncs contains a function, called '%s', which does not return a vector of proper length as result.", rvfName ) )
      if( ! is( attr(testres, "size" ), "numeric" ) )
         return( sprintf( "rawVarFuncs contains a function, called '%s', which return a result with numeric 'size' attribute.", rvfName ) )
   }
   if( !is.integer( counts(object) ) )
      return( "the count data is not in integer mode" )
   if( any( counts(object) < 0 ) )
      return( "the count data contains negative values" )
   TRUE
} )

counts <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   assayData(cds)[["counts"]]
}   
   
sizeFactors <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   sf <- pData(cds)$`sizeFactor`
   names( sf ) <- colnames( counts(cds) )
   sf
}   
   
`sizeFactors<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   pData(cds)$`sizeFactor` <- value
   validObject( cds )
   cds
}   

conditions <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( cds@multivariateConditions )
      stop( "The 'conditions' accessor is only for simple single-factor conditions, but your have specified multivariate conditions. Access them via 'pData'." )
   conds <- pData(cds)$`condition`
   names( conds ) <- colnames( counts(cds) )
   conds
}   
   
`conditions<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( cds@multivariateConditions )
      stop( "The 'conditions' accessor is only for simple single-factor conditions, but your have specified multivariate conditions. Access them via 'pData'." )
   pData(cds)$`condition` <- factor( value )
   validObject( cds )
   cds
}   

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

rawVarFuncTable <- function( cds ) {
   stopifnot( is( cds, "CountDataSet" ) )
   cds@rawVarFuncTable
}   

`rawVarFuncTable<-` <- function( cds, value ) {
   stopifnot( is( cds, "CountDataSet" ) )
   if( is.null( names(value) ) )
      names( value ) <- names( rawVarFuncTable(cds) )
   cds@rawVarFuncTable <- value
   validObject( cds )   
   cds
}   

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
