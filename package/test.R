#library(DESeq)
#source("smutans/R/core.R")
#smutans.prepareData835NPP ()
#q("no")
# library(DESeq)
load("smutans/data/sm835Genes.RData")
source("smutans/R/class_and_slots.R")
source("smutans/R/core.R")
source("smutans/R/methods.R")

head( counts(sm835Genes) )
pData( sm835Genes )

npAndP <- newSmutans( sm835Genes, title="835NP vs. 835P" )
npAndP <- smutans.de2( npAndP, type="ua159", 
                       condA="835np", condB="835p" )
smutans.plotDispersionEstimates(npAndP)

#ua159only <- smutans.de2type( ua159only, type="ua159", 
                              #condA="glucose", condB="galactose" )
#allSample <- smutans.de2( allSample )
#
#padj.ua159only <- smutans.padj( ua159only )
#padj.allSample <- smutans.padj( allSample )
#tab = table( "ua159 only" = padj.ua159only < .1, "all samples" = padj.allSample < .1 )
#addmargins( tab )
#
#pval.ua159only <- smutans.pval( ua159only )
#pval.allSample <- smutans.pval( allSample )
#bottom = function(x, theta=1e-12) pmax(x, theta)
#plot( bottom(pval.ua159only), bottom(pval.allSample), log="xy", pch=20, cex=.3 )
#abline(a=0, b=1, col="blue")
#
#allSampleWithoutType <- newSmutans( smutansGenes )
#allSampleWithoutType <- smutans.de2WithoutType( allSampleWithoutType,
#                                                condA="glucose", 
#                                                condB="galactose" )
#padj.allSampleWithoutType <- smutans.padj( allSampleWithoutType )
#tab2 <- table(
   #`all samples simple` = padj.allSampleWithoutType < 0.1,
   #`all samples GLM`    = padj.allSample < 0.1 )
#addmargins( tab2 )

# smutans.de2TypeHeatmap( ua159only )
# smutans.de2TypeLfc( ua159only, condA="glucose", condB="galactose" )
#smutans.de2Clust( allSample )

