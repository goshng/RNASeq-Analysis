# smutans.de2Clust( npAndUA159 )
# smutans.de2Clust( pAndUA159 )
# smutans.de2TypeHeatmap( npAndP )
# smutans.de2TypeHeatmap( npAndUA159 )
# smutans.de2TypeHeatmap( pAndUA159 )
#smutans.de2Goseq ( npAndP, qval=1e-2, 
#                   feature.gene=smutans.feature.genes,
#                   go.genes=smutans.go.genes,
#                   cat.desc=smutans.cat.desc )


library(GenomicRanges)
#library(DESeq)
source("smutans/R/core.R")
smutans.prepareGoseq( )
#smutans.prepareData835NPP ()
#smutans.prepareData34 ()
#smutans.prepareDataOMZ175( ) 
#smutans.prepareDataUA159CSP( ) 
#smutans.prepareDataSMU86CSP( )
#smutans.prepareDataUA159TW1()
#smutans.prepareTranscript( )
q("no")
library(DESeq)
library(goseq)
load("smutans/data/smutansGenes.RData")
load("smutans/data/smutansGenes2.RData")
load("smutans/data/sm835Genes.RData")
load("smutans/data/sm34Genes.RData")
load("smutans/data/smutans.cat.desc.RData")
load("smutans/data/smutans.genes.criteria.RData")
load("smutans/data/smutans.go.genes.RData")
load("smutans/data/smutansData.txPileup.RData")
source("smutans/R/class_and_slots.R")
source("smutans/R/core.R")
source("smutans/R/methods.R")
# allSample2 <- newSmutans( smutansGenes2 )
# smutans.de2Clust( allSample2 )
#ua159only <- newSmutans( smutansGenes2, title="UA159 Only" )
#ua159only <- smutans.de2( ua159only, type="ua159", 
#                          condA="glucose", condB="galactose" )
#allSample <- newSmutans( smutansGenes2, 
#                         title="All Samples for Comparing Conditions" )
#allSample <- smutans.de2( allSample )
#tabUA159onlyAndAllSample <- smutans.compareResult( ua159only, 
#                                                   allSample )
#addmargins( tabUA159onlyAndAllSample )
allSampleType <- newSmutans( smutansGenes2, title="All samples for type" )
allSampleType <- smutans.de2( allSampleType, compareCondition="no" )
print(
smutans.de2List( allSampleType )[grep("locus",rownames(smutans.de2List(allSampleType))),]
)
#smutans.plotDispersionEstimates(ua159only)
#smutans.plotDiffExp(ua159only)
#smutans.de2List( ua159only )

#smutans.mannwhitney( qval=0.05,
#                   genes.criteria=smutans.genes.criteria,
#                   go.genes=smutans.go.genes,
#                   cat.desc=smutans.cat.desc )

# smutans.de2Clust( allSample )
#allSample <- newSmutans( smutansGenes, title="All Samples for Comparsing Conditions" )
#allSample <- smutans.de2( allSample )
#load("smutans/data/smutans.cat.desc.RData")
#load("smutans/data/smutans.feature.genes.RData")
#load("smutans/data/smutans.go.genes.RData")
#smutans.de2Goseq ( allSample, qval=1e-2, 
#                   feature.gene=smutans.feature.genes,
#                   go.genes=smutans.go.genes,
#                   cat.desc=smutans.cat.desc )


#allSample <- newSmutans( smutansGenes, title="All Samples for Comparsing Conditions" )
#allSample <- smutans.de2( allSample )
#ua159only <- newSmutans( smutansGenes, title="UA159 Only" )
#ua159only <- smutans.de2( ua159only, type="ua159", 
#                          condA="glucose", condB="galactose" )
#tw1only <- newSmutans( smutansGenes, title="TW1 Only" )
#tw1only <- smutans.de2type( tw1only, type="tw1", 
#                            condA="glucose", condB="galactose" )
#npAndP <- newSmutans( sm34Genes, title="835NP vs. 835P" )
#npAndP <- smutans.de2( npAndP, type="ua159", 
#                       condA="835NP", condB="835P" )

#head( counts(sm835Genes) )
#pData( sm835Genes )

#npAndP <- newSmutans( sm835Genes, title="835NP vs. 835P" )
#npAndP <- smutans.de2( npAndP, type="ua159", 
#                       condA="835np", condB="835p" )
#allSample <- newSmutans( smutansGenes )
#allSample <- smutans.de2( allSample )
#npAndUA159 <- newSmutans( sm835Genes, title="835NP vs. UA159" )
#npAndUA159 <- smutans.de2( npAndUA159, type="ua159", 
#                           condA="835np", condB="ua159nocsp" )
#pAndUA159 <- newSmutans( sm835Genes, title="835P vs. UA159" )
#pAndUA159 <- smutans.de2( pAndUA159, type="ua159", 
#                          condA="835p", condB="ua159nocsp" )
# allSample
# ua159only
#smutans.de2Goseq ( ua159only, qval=1e-2, 
#                   feature.gene=smutans.feature.genes,
#                   go.genes=smutans.go.genes,
#                   cat.desc=smutans.cat.desc )

#load("smutans/data/smutans.cat.desc.RData")
#load("smutans/data/smutans.feature.genes.RData")
#load("smutans/data/smutans.go.genes.RData")

#
#feature.genes = smutans.feature.genes
#length.genes = feature.genes$V3 - feature.genes$V2
#assayed.genes = feature.genes$V4
#
##de.genes = scan(file="de.gene.2",what="character",quiet=T)
#de.genes <- smutans.de2Genes ( npAndP, qval=1e-5 )
#gene.vector = as.integer(assayed.genes %in% de.genes)
#names(gene.vector) = assayed.genes
#
#go.genes = smutans.go.genes
#cat.desc = smutans.cat.desc
#rm(feature.genes)
#pwf = nullp( gene.vector, bias.data=length.genes, plot.fit=FALSE )
#
## go.hypergeometric = goseq(pwf,gene2cat=go.genes,method="Hypergeometric") # No length bias correction
#go = goseq(pwf,gene2cat=go.genes,method="Wallenius") # Length bias correction - Approximation
## go.sample = goseq(pwf,gene2cat=go.genes,method="Sampling",repcnt=10000) # Length bias correction - Sampling
#go.fdr = go[p.adjust(go$over_represented_pvalue,method="BH")<.05,]
#
#for (i in go.fdr$category)
#{
#  cat(as.character(cat.desc$V3[cat.desc$V1==i]),"\n")
#}
#
#







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
sm835result <- function () {
  library(smutans)
  data( "sm835Genes" )
  data( "smutans.feature.genes" )
  data( "smutans.go.genes" )
  data( "smutans.cat.desc" )
  npAndP <- newSmutans( sm835Genes, title="835NP vs. 835P" )
  npAndP <- smutans.de2( npAndP, type="ua159", 
                         condA="835np", condB="835p" )
  npAndUA159 <- newSmutans( sm835Genes, title="835NP vs. UA159" )
  npAndUA159 <- smutans.de2( npAndUA159, type="ua159", 
                             condA="835np", condB="ua159nocsp" )
  pAndUA159 <- newSmutans( sm835Genes, title="835P vs. UA159" )
  pAndUA159 <- smutans.de2( pAndUA159, type="ua159", 
                            condA="835p", condB="ua159nocsp" )
  smutans.de2List( npAndP, "npAndP.csv" )
  smutans.de2List( npAndUA159, "npAndUA159.csv" )
  smutans.de2List( pAndUA159, "pAndUA159.csv" )
  smutans.de2Clust( npAndP )
  # smutans.de2Clust( npAndP, margins=c(7,7) ) this is for v0.1.3
}
