library(DESeq)
cutoffpval <- 0.05
countsFile <- "output/omz/1/bwa/count.txt"
countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- countsTable$gene
countsTable <- countsTable[,-1]
y <- length(countsTable[,1])
cat ("Total number of genes is ", y, "\n",sep="")

removedCountsTable <- countsTable[countsTable[,1]+countsTable[,2]+countsTable[,3]+countsTable[,4]+countsTable[,5]+countsTable[,6]<=15,]
y <- length(removedCountsTable[,1])
cat ("Number of genes removed is ", y, "\n",sep="")
countsTable <- countsTable[countsTable[,1]+countsTable[,2]+countsTable[,3]+countsTable[,4]+countsTable[,5]+countsTable[,6]>15,]

conds <- c("OMZ175", "OMZHKRR", "OMZ175", "OMZ175", "OMZHKRR", "OMZHKRR")

cds <- newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors(cds)

cds <- estimateVarianceFunctions(cds)
res <- nbinomTest(cds, "OMZ175", "OMZHKRR")
# 1 3 4 vs 2 5 6
# 1  "OMZ175"
# 2  "OMZHKRR"
# 3  "OMZ175"
# 4  "OMZ175"
# 5  "OMZHKRR"
# 6  "OMZHKRR"
      
plotDE <- function( res )
{
   plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.5, 
        col = ifelse( res$padj < cutoffpval, "red", "black" ) )
}

cat ("SCV file: omz-de.ps\n")
postscript ("output/omz/1/bwa/omz-de.ps",  width=10, height=10, horizontal = FALSE, onefile = FALSE,
paper = "special")
plotDE(res)
y <- dev.off()

cat ("SCV file: omz-scv.ps\n")
postscript ("output/omz/1/bwa/omz-scv.ps",  width=10, height=10, horizontal = FALSE, onefile = FALSE,
paper = "special")
scvPlot(cds)
y <- dev.off()

resSig <- res[res$padj < cutoffpval,]
y <- length(countsTable[,1])
cat ("Total number of genes analyzed for differential gene expression is ", y, "\n",sep="")
y <- length(resSig[,1])
cat ("The number of differentially expressed genes is ", y, "\n",sep="")
y <- length(res[res$baseMean < 100,1])
cat ("The number of genes with less than 100 of mean counts is ", y, "\n",sep="")
cat ("\nList of differentially expressed genes\n")
options(width = 1000) 
print(resSig[order(resSig$pval),])
cat ("\nThe same table as above but in order of down-regulated first of the differentially expressed genes\n")
print(resSig[order(resSig$foldChange,-resSig$baseMean),])
cat ("\nThe same table as above but in order of up-regulated first of the differentially expressed genes\n")
print(resSig[order(-resSig$foldChange,-resSig$baseMean),])

cat ("\nList of all of the genes\n")
print(res)
