
# Read a table
# newCountDataSet <- function()
{
   countData <- countsTable 
   countData <- as.matrix( countData )
   if( any( round( countData ) != countData ) )
      stop( "The countData is not integer." )
   mode( countData ) <- "integer"
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
