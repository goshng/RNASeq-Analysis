#installBaySeq <- function ()
{
source("http://www.bioconductor.org/biocLite.R")
biocLite("baySeq")
}

#loadBaySeq <- function ()
{
library(baySeq)
# library(snow)
# cl <- makeCluster(4, "SOCK")
cl <- NULL
}

#loadSimData <- function ()
{
# countsFile <- "/Users/goshng/Documents/Projects/rnaseq-analysis/output/omz/1/bwa/count.txt"
countsFile <- "/Users/goshng/Documents/Projects/rnaseq/output/omz/1/bwa/count.txt"
countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- countsTable$gene
countsTable <- countsTable[,-1]
simCount <- as.matrix(countsTable)

# simCount[1:10,]
libsizes <- colSums(simCount)
replicates <- c(1,2)
groups <- list(NDE = c(1,1), DE = c(1,2))
CD <- new("countData", data = simCount, replicates = replicates, libsizes = as.integer(libsizes), groups = groups)
#plotMA.CD(CD, samplesA = 1, samplesB = 2, col = c(rep("red", 100), rep("black", 900)))
}

bayseqPoisson <- function ()
{
CDP.Poi <- getPriors.Pois(CD, samplesize = 20, takemean = TRUE, cl = cl)
CDP.Poi@priors
}

bayseqPoisson2 <- function ()
{
CDPost.Poi <- getLikelihoods.Pois(CDP.Poi, pET = "BIC", cl = cl)
CDPost.Poi@estProps
CDPost.Poi@posteriors[1:10,]
CDPost.Poi@posteriors[101:110,]
}

#bayseqNB <- function ()
{
CDP.NBML <- getPriors.NB(CD, samplesize = 10000, estimation = "QL", cl = cl)
CDPost.NBML <- getLikelihoods.NB(CDP.NBML, pET = 'BIC', cl = cl)
#CDPost.NBML@estProps
#CDPost.NBML@posteriors[1:10,]
#CDPost.NBML@posteriors[101:110,]
}


