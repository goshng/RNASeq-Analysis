# Rscript R/deseq.R count.txt species conditionA conditionB
.Last <- function() {
  cat("\n")
}

plotDispEsts <- function( cds )
{
   plot(
      rowMeans( counts( cds, normalized=TRUE ) ),
      fitInfo(cds)$perGeneDispEsts,
      pch = '.', log="xy" )
   xg <- 10^seq( -.5, 5, length.out=300 )
   lines( xg, fitInfo(cds)$dispFun( xg ), col="red" )
}

plotDE <- function( res )
   plot(
      res$baseMean,
      res$log2FoldChange,
      log="x", pch=20, cex=.3,
      col = ifelse( res$padj < .1, "red", "black" ) )

args <- commandArgs(trailingOnly = TRUE)
#if (length(args) != 2)
#{
#  cat ("Rscript R/deseq-cornell-ua159-tw1.R output/cornell/1/bwa/count.txt 0.01\n")
#  quit("yes")
#}
#countsFile <- args[1]
#cutoffpval <- as.numeric(args[2])
#print(countsFile)
#print(cutoffpval)
countsFile <- "output/cornell/1/bwa/count.txt"
cutoffpval <- 0.01

prepare.input.data <- function ()
{
library(DESeq)
# Count data preparation
conds <- scan(file=paste(countsFile,"index",sep="."), what="character")
subconds <- conds %in% c("UA159GLU", "UA159GAL", "TW1GLU", "TW1GAL")
countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- countsTable$gene
countsTable <- countsTable[,-1]
countsTable <- countsTable[,subconds]
conds <- conds[subconds]
conds.type <- conds %in% c("UA159GLU", "UA159GAL")
factor.type <- rep("ua159",length(conds))
factor.type[!conds.type] <- "tw1"
conds.condition <- conds %in% c("UA159GLU", "TW1GLU")
factor.condition <- rep("glucose",length(conds))
factor.condition[!conds.condition] <- "galactose"
samples <- data.frame(type=factor.type,condition=factor.condition)
rownames(samples) <- colnames(countsTable) 

design <- samples
smutansGenes <- newCountDataSet( countsTable, design )
expdata = new("MIAME", 
   name="S. mutans UA159, TW1, Glucose, and Galactose", 
   lab="University of Florida, and Cornell University", 
   contact="Drs. Robert Burne, Michael Stanhope, and Adam Siepel", 
   title="Streptococcus mutans RNA-Seq Studies", 
   url="http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=XXX", 
   abstract="RNA-seq of 13 biological replicates from Streptococcus mutans")
pubMedIds(expdata) <- "999999999"
experimentData(smutansGenes) <- expdata
save(smutansGenes, file=file.path(".", "smutansGenes.RData"))
quit("no")
}

# library( "smutans" )
# data( "smutansGenes" )
# head( counts(smutansGenes) )
# pData( smutansGenes )

ua159Samples <- pData(smutansGenes)$type == "ua159"
countsTable <- counts(smutansGenes)[ , ua159Samples ]
conds <- pData(smutansGenes)$condition[ ua159Samples ]
cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
# plotDispEsts(cds)
#res <- nbinomTest( cds, "glucose", "galactose" )
# plotDE( res )

design <- pData(smutansGenes)[,c("type","condition")]
# Inference
fullCountsTable <- counts( smutansGenes )
cdsFull <- newCountDataSet( fullCountsTable, design )
cdsFull <- estimateSizeFactors( cdsFull )
cdsFull <- estimateDispersions( cdsFull )
# plotDispEsts( cdsFull )
#fit1 <- fitNbinomGLMs( cdsFull, count ~ type + condition )
#fit0 <- fitNbinomGLMs( cdsFull, count ~ type  )
#pvalsGLM <- nbinomGLMTest( fit1, fit0 )
#padjGLM <- p.adjust( pvalsGLM, method="BH" )

tab = table( "ua159 only" = res$padj < .1, "all samples" = padjGLM < .1 )
addmargins( tab )

#cdsFullB <- newCountDataSet( fullCountsTable, design$condition )
#cdsFullB <- estimateSizeFactors( cdsFullB )
#cdsFullB <- estimateDispersions( cdsFullB )
#resFullB <- nbinomTest( cdsFullB, "glucose", "galactose" )
tab2 <- table(
   `all samples simple` = resFullB$padj < 0.1,
   `all samples GLM`    = padjGLM < 0.1 )
addmargins( tab2 )

cdsBlind <- estimateDispersions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cdsBlind )
mod_lfc <- (rowMeans( vsd[, conditions(cds)=="galactose", drop=FALSE] ) -
            rowMeans( vsd[, conditions(cds)=="glucose", drop=FALSE] ))

lfc <- res$log2FoldChange
finite <- is.finite(lfc)
table(as.character(lfc[!finite]), useNA="always")

largeNumber <- 10
lfc <- ifelse(finite, lfc, sign(lfc) * largeNumber)

plot( lfc, mod_lfc, pch=20, cex=.3, col = ifelse( finite, "#80808040", "red" ) )
abline( a=0, b=1, col="#40404040" )


select <- order(res$pval)[1:40]
colors <- colorRampPalette(c("white","darkblue"))(100)
heatmap( vsd[select,],
         col = colors, scale = "none")

cdsFullBlind <- estimateDispersions( cdsFull, method = "blind" )
vsdFull <- getVarianceStabilizedData( cdsFullBlind )
dists <- dist( t( vsdFull ) )
heatmap( as.matrix( dists ),
   symm=TRUE, scale="none", margins=c(10,10),
   col = colorRampPalette(c("darkblue","white"))(100),
   labRow = paste( pData(cdsFullBlind)$condition, pData(cdsFullBlind)$type ) )
