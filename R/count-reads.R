library(Rsamtools)
library(GenomicFeatures)
library(rtracklayer)

bedFile <- "test/pichon/pichon.bed"
x <- import.bed(bedFile)
y <- as(x,"GRanges")
count.table <- data.frame(gene=elementMetadata(y)$name)
for (i in c(1,2,3,4,5,15,16,17,18,19,20,21,22,24,25,26,27,23,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55))
{
  bamFile <- sprintf("output/ua159/1/bwa/FASTQ%03d.sorted.bam",i)
  print(bamFile)
  bv <- readBamGappedAlignments(bamFile,
                                use.names=TRUE,
                                param=ScanBamParam(what=c("mapq")))
  olap <- summarizeOverlaps(y,bv,ignore.strand=FALSE,mode="IntersectionStrict")
  cl <- assays(olap)$counts[,1]

  count.table <- data.frame(count.table,cl)
  colnames(count.table)[ncol(count.table)] <- paste("X",i,sep="")
}
colnames(count.table) <- sub("X","",colnames(count.table))
count.table.file <- sprintf("%s.cl",bedFile)
write.table(count.table,file=count.table.file,quote=FALSE,sep="\t",row.names=FALSE)
print(paste("See cl file",count.table.file))
