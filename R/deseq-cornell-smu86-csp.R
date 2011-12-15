library(DESeq)
cutoffpval <- 0.05
countsFile <- "output/cornell/1/bwa/count.txt"
countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- countsTable$gene
countsTable <- countsTable[,-1]
y <- length(countsTable[,1])
cat ("Total number of genes is ", y, "\n",sep="")

removedCountsTable <- countsTable[countsTable[,31]+countsTable[,32]+countsTable[,33]+countsTable[,34]<=10,]
y <- length(removedCountsTable[,1])
cat ("Number of genes removed is ", y, "\n",sep="")
countsTable <- countsTable[countsTable[,31]+countsTable[,32]+countsTable[,33]+countsTable[,34]>10,]

conds <- c("UA159GLU", "TW1GLU", "TW1GAL", "OMZ175", "OMZHKRR", "OMZ175", "OMZ175", "OMZHKRR", "OMZHKRR", "UA159GLU", "UA159GLU", "UA159GAL", "UA159GAL", "TW1GLU", "TW1GLU", "TW1GAL", "TW1GAL", "UA159GAL", "UA159GAL", "835NP", "835NP", "835NP", "835P", "835P", "835P", "UA159noCSP", "UA159noCSP", "UA159CSP", "UA159CSP", "UA159CSP", "Smu86noCSP", "Smu86noCSP", "Smu86CSP", "Smu86CSP")

cds <- newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors(cds)

cds <- estimateVarianceFunctions(cds)
res <- nbinomTest(cds, "Smu86noCSP", "Smu86CSP")
# 31 32 vs 33 34
# 1  "UA159GLU"
# 2  "TW1GLU"
# 3  "TW1GAL"
# 4  "OMZ175"
# 5  "OMZHKRR"
# 6  "OMZ175"
# 7  "OMZ175"
# 8  "OMZHKRR"
# 9  "OMZHKRR"
# 10 "UA159GLU"
# 11 "UA159GLU"
# 12 "UA159GAL"
# 13 "UA159GAL"
# 14 "TW1GLU"
# 15 "TW1GLU"
# 16 "TW1GAL"
# 17 "TW1GAL"
# -------------
# 18 "UA159GAL"
# 19 "UA159GAL"
# 20 "835NP"
# 21 "835NP"
# 22 "835NP"
# 23 "835P"
# 24 "835P"
# 25 "835P"
# 26 "UA159noCSP"
# 27 "UA159noCSP"
# 28 "UA159CSP"
# 29 "UA159CSP"
# 30 "UA159CSP"
# 31 "Smu86noCSP"
# 32 "Smu86noCSP"
# 33 "Smu86CSP"
# 34 "Smu86CSP"
      
plotDE <- function( res )
{
   plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.5, 
        col = ifelse( res$padj < cutoffpval, "red", "black" ) )
}

cat ("SCV file: smu86-csp-de.ps\n")
postscript ("output/cornell/1/bwa/smu86-csp-de.ps",  width=10, height=10, horizontal = FALSE, onefile = FALSE,
paper = "special")
plotDE(res)
y <- dev.off()

cat ("SCV file: smu86-csp-scv.ps\n")
postscript ("output/cornell/1/bwa/smu86-csp-scv.ps",  width=10, height=10, horizontal = FALSE, onefile = FALSE,
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
