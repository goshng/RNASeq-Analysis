###################################################
### chunk number 1: 
###################################################
#line 190 "vignettes/DEGseq/inst/doc/DEGseq.Rnw"
library(DEGseq)
step1 <- function()
{
  geneExpFile <- system.file("extdata", "GeneExpExample5000.txt", package="DEGseq")
  geneExpMatrix1 <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(7,9,12,15,18))
  geneExpMatrix2 <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(8,10,11,13,16))
  write.table(geneExpMatrix1[30:31,],row.names=FALSE)
  write.table(geneExpMatrix2[30:31,],row.names=FALSE)
}


###################################################
### chunk number 2: 
###################################################
#line 204 "vignettes/DEGseq/inst/doc/DEGseq.Rnw"
step2 <- function()
{
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow=TRUE))
par(mar=c(2, 2, 2, 2))
DEGexp(geneExpMatrix1=geneExpMatrix1, geneCol1=1, expCol1=c(2,3,4,5,6), groupLabel1="kidney",
       geneExpMatrix2=geneExpMatrix2, geneCol2=1, expCol2=c(2,3,4,5,6), groupLabel2="liver",
       method="MARS")
}

###################################################
### chunk number 3: 
###################################################
#line 223 "vignettes/DEGseq/inst/doc/DEGseq.Rnw"
step3 <- function()
{
layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow=TRUE))
par(mar=c(2, 2, 2, 2))
DEGexp(geneExpMatrix1=geneExpMatrix1, expCol1=2, groupLabel1="kidneyR1L1", 
       geneExpMatrix2=geneExpMatrix2, expCol2=2, groupLabel2="liverR1L2",
       replicateExpMatrix1=geneExpMatrix1, expColR1=3, replicateExpMatrix2=geneExpMatrix1, expColR2=4,
       replicateLabel1="kidneyR1L3", replicateLabel2="kidneyR1L7",
       method="MATR")
}

###################################################
### chunk number 4: 
###################################################
#line 245 "vignettes/DEGseq/inst/doc/DEGseq.Rnw"
step4 <- function()
{
  kidneyR1L1 <- system.file("extdata", "kidneyChr21.bed.txt", package="DEGseq")
  liverR1L2  <- system.file("extdata", "liverChr21.bed.txt", package="DEGseq")
  refFlat    <- system.file("extdata", "refFlatChr21.txt", package="DEGseq")
  mapResultBatch1 <- c(kidneyR1L1)  ## only use the data from kidneyR1L1 and liverR1L2
  mapResultBatch2 <- c(liverR1L2)
  outputDir <- file.path(tempdir(), "DEGseqExample")
  DEGseq(mapResultBatch1, mapResultBatch2, fileFormat="bed", refFlat=refFlat,
         outputDir=outputDir, method="LRT")
}

###################################################
### chunk number 5: 
###################################################
#line 264 "vignettes/DEGseq/inst/doc/DEGseq.Rnw"
step5 <- function()
{
  geneExpFile <- system.file("extdata", "GeneExpExample1000.txt", package="DEGseq")
  set.seed(100)
  geneExpFile1 <- geneExpFile 
  geneExpFile2 <- geneExpFile
  output <- file.path(tempdir(), "samWrapperOut.txt")
  expCol1=c(7,9,12,15,18)
  expCol2=c(8,10,11,13,16)
  measure1=c(-1,-2,-3,-4,-5)
  measure2=c(1,2,3,4,5)
  samWrapper(geneExpFile1=geneExpFile1, geneCol1=1, expCol1=expCol1, measure1=measure1,
             geneExpFile2=geneExpFile2, geneCol2=1, expCol2=expCol2, measure2=measure2,
             nperms=100, min.foldchange=2, max.qValue=1e-4, output=output, paired=TRUE)
}

###################################################
### chunk number 6: 
###################################################
#line 288 "vignettes/DEGseq/inst/doc/DEGseq.Rnw"
step6 <- function()
{
  kidneyR1L1 <- system.file("extdata", "kidneyChr21.bed.txt", package="DEGseq")
  refFlat    <- system.file("extdata", "refFlatChr21.txt", package="DEGseq")
  mapResultBatch <- c(kidneyR1L1)
  output <- file.path(tempdir(), "kidneyChr21.bed.exp")
  exp <- getGeneExp(mapResultBatch, refFlat=refFlat, output=output)
  write.table(exp[30:32,], row.names=FALSE)
}

###################################################
### chunk number 7: 
###################################################
#line 302 "vignettes/DEGseq/inst/doc/DEGseq.Rnw"
step7 <- function()
{
  geneExpFile <- system.file("extdata", "GeneExpExample1000.txt", package="DEGseq")
  exp <- readGeneExp(file=geneExpFile, geneCol=1, valCol=c(7,9,12,15,18,8,10,11,13,16))
  write.table(exp[30:32,], row.names=FALSE)
}

###################################################
### chunk number 8: 
###################################################
#line 310 "vignettes/DEGseq/inst/doc/DEGseq.Rnw"
step8 <- function()
{
  kidneyR1L1 <- system.file("extdata", "kidneyChr21Bowtie.txt", package="DEGseq")
  liverR1L2  <- system.file("extdata", "liverChr21Bowtie.txt", package="DEGseq")
  refFlat    <- system.file("extdata", "refFlatChr21.txt", package="DEGseq")
  kidneyR1L1_aln <- ShortRead::readAligned(dirname(kidneyR1L1), basename(kidneyR1L1), type="Bowtie")
  liverR1L2_aln <- ShortRead::readAligned(dirname(liverR1L2), basename(liverR1L2), type="Bowtie")
  alnBatch1 <- list(kidneyR1L1_aln)  ## only use the data from kidneyR1L1 and liverR1L2
  alnBatch2 <- list(liverR1L2_aln)
  outputDir <- file.path(tempdir(), "DEGseqAlnExample")
  DEGseq.aln(alnBatch1, alnBatch2, refFlat=refFlat, outputDir=outputDir, method="MARS")
  cat("outputDir:", outputDir, "\n")
}

###################################################
### chunk number 9: 
###################################################
#line 323 "vignettes/DEGseq/inst/doc/DEGseq.Rnw"
step9 <- function()
{
  alnBatch <- list(kidneyR1L1_aln, liverR1L2_aln)
  output <- file.path(tempdir(), "kidneyChr21Bowtie.exp")
  exp <- getGeneExp.aln(alnBatch, refFlat=refFlat, output=output)
  write.table(exp[30:32,], row.names=FALSE)
  cat("output: ", output, "\n")
}

