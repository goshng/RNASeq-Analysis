library(DESeq)
library(smutans)
smomzGenes <- 
  readSmutans (countsFile="/Volumes/Elements/Documents/Projects/RNASeq-Analysis/output/tux/1/bwa/count-NC_004350.txt",
               indexFile="/Volumes/Elements/Documents/Projects/RNASeq-Analysis/output/tux/1/bwa/count-NC_004350.txt.index",
               condition = c("UA159", "NotUA159"), 
               firstFactorLabel = c("ua159"), 
               secondFactorLabel = c("UA159", "NotUA159"),
               name="S. mutans UA159")

head( counts(smomzGenes), n=1 )
pData( smomzGenes )

smomzGenes <- estimateSizeFactors( smomzGenes )

# Size Factors
print(pData( smomzGenes )$sizeFactor)

# Estimate dispersion
smomzGenes <- estimateDispersions( smomzGenes )
print(str( fitInfo(smomzGenes) ))
quit("no")
# END

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

#
condA <- "UA159"
colA <- pData(cds)$condition=="UA159"
bmv <- getBaseMeansAndVariances( counts(cds)[,colA],sizeFactors(cds)[colA] )
rawScvA <- fData(cds)[ , paste( "disp", dispTable(cds), sep="_" ) ]

countsA <- counts(cds)[,colA]
sizeFactorsA <- sizeFactors(cds)[colA]
dispsA <- rawScvA

kAs <- rowSums( cbind(countsA) )
# or 
kAs <- rowSums( countsA )
# Computes the average and variance of ???.
mus <- rowMeans( cbind( t( t( countsA ) / sizeFactorsA ) ) )

# fullVarsA is \sigma^2, which is \mu + 
fullVarsA <- pmax( mus * sum( sizeFactorsA ) + dispsA * mus^2 * sum(sizeFactorsA^2), 
                   mus * sum( sizeFactorsA ) * (1+1e-8) )

# From R function of negative binomial, \sigma^2=\mu+\mu^2/size.
# d=1/size, which is sumDispsA, is equal to (\sigma^2 - \mu)/\mu^2, which is the
# following line. So, mus * sum( sizeFactorsA ) is \mu?
sumDispsA <- ( fullVarsA - mus * sum( sizeFactorsA ) ) / ( mus * sum( sizeFactorsA ) )^2




