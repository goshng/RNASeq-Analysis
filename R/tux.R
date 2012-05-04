library(DESeq)
library(smutans)
smomzGenes <- 
  readSmutans (countsFile="output/ua159/1/bwa/count-NC_004350.txt",
               indexFile ="output/ua159/1/bwa/count-NC_004350.txt.index",
               condition = c("UA159", "NotUA159"), 
               firstFactorLabel = c("ua159"), 
               secondFactorLabel = c("UA159", "NotUA159"),
               name="S. mutans UA159")

head( counts(smomzGenes), n=1 )
pData( smomzGenes )

smomzGenes <- estimateSizeFactors( smomzGenes )

# Size Factors
pData( smomzGenes )$sizeFactor

# Estimate dispersion
smomzGenes <- estimateDispersions( smomzGenes )
str( fitInfo(smomzGenes) )
cds <- smomzGenes
max(log10(rowMeans( counts( cds, normalized=TRUE ) )))
xg <- 10^seq( -.5, max(log10(rowMeans( counts( cds, normalized=TRUE ) ))) + 1, length.out=300 )
fitInfo(cds)$dispFun( xg )

# Plot the dispersion estimate
plotDispEsts <- function( cds ) {
  plot( rowMeans( counts( cds, normalized=TRUE ) ), 
        fitInfo(cds)$perGeneDispEsts,
        pch = '.', log="xy" )
  xg <- 10^seq( -.5, max(log10(rowMeans( counts( cds, normalized=TRUE ) ))) + 1, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}
plotDispEsts( cds )
head( fData(cds) )

head(rowMeans( counts( cds, normalized=TRUE ) ))
# How can we use the fitted dispersion function?
# What are estimated for the negative binomial?

# The fitted function
# List of 5
# $ perGeneDispEsts: num [1:1931] 4.33e-05 1.28e-02 1.44e-02 2.56e-02 2.43e-03 ...
# $ dispFunc       :function (q)  
#  ..- attr(*, "coefficients")= Named num [1:2] 0.0269 6.2076
#  .. ..- attr(*, "names")= chr [1:2] "asymptDisp" "extraPois"
#  ..- attr(*, "fitType")= chr "parametric"
# $ fittedDispEsts : num [1:1931] 0.029 0.0293 0.0764 0.0279 0.0319 ...
# $ df             : int 2
# $ sharingMode    : chr "maximum"
