library(GenomicFeatures)
library(rtracklayer)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
  cat ("Rscript xbasegff2txdb.R file.gf3\\n")
  quit("no")
}
gffFile <- args[1]
gffFileBasename <- unlist(strsplit(gffFile,"\\."))[1]

sm.gff <- import.gff3(gffFile)
sm.gene <- sm.gff[sm.gff$type=="gene",]
sm.CDS <- sm.gff[sm.gff$type=="CDS",]
sm.source <- sm.gff[sm.gff$type=="region",]

transcripts <- 
  data.frame( tx_id=seq(length(sm.gene$locus_tag)),
              tx_name=sm.gene$locus_tag,
              tx_chrom=sm.gene$space, 
              tx_strand=sm.gene$strand,
              tx_start=start(sm.gene),
              tx_end=end(sm.gene) )

# Compare CDS and gene locus tags
#################################################################
# 1. Using parent information
x <- c()
y <- c()
z <- unlist(sm.CDS$Parent)
# z <- unlist(sm.CDS$locus_tag)
for (i in sm.gene$ID) {
  if (sum(z == i) > 0) {
    if (length(start(sm.CDS)[z == i]) != 1) {               
      print(paste("Check",gffFile))
      print(start(sm.CDS)[z == i])                          
      stop(paste("There are multiple CDS for gene",i))
    }
    x <- c(x,start(sm.CDS)[z == i])
    y <- c(y,end(sm.CDS)[z == i])
  } else {
    x <- c(x,NA)
    y <- c(y,NA)
  }
}
grCDS.start <- x
grCDS.end <- y
rm(i,x,y,z)
##################################################################
# 2. FIXME: Using positions?
x <- c()
y <- c()
# z <- unlist(sm.CDS$Parent)
z <- unlist(sm.CDS$locus_tag)
for (i in sm.gene$ID) {
  if (sum(z == i) > 0) {
    if (length(start(sm.CDS)[z == i]) != 1) {               
      print(paste("Check",gffFile))
      print(start(sm.CDS)[z == i])                          
      stop(paste("There are multiple CDS for",i))
    }
    x <- c(x,start(sm.CDS)[z == i])
    y <- c(y,end(sm.CDS)[z == i])
  } else {
    x <- c(x,NA)
    y <- c(y,NA)
  }
}
grCDS.start <- x
grCDS.end <- y
rm(i,x,y,z)
##################################################################

# 
splicings <- 
  data.frame( tx_id=seq(length(sm.gene$locus_tag)),
              exon_rank=seq(length(sm.gene$locus_tag)),
              exon_start=start(sm.gene),
              exon_end=end(sm.gene),
              exon_name=sm.gene$locus_tag,
              cds_start=grCDS.start,
              cds_end=grCDS.end )

chrominfo <-
  data.frame( chrom=names(sm.gene),
              length=end(sm.source),
              is_circular=FALSE )

txdb <- makeTranscriptDb (transcripts, splicings, chrominfo=chrominfo)
saveFeatures(txdb,file=paste(gffFileBasename,"txdb",sep="."))
print(paste("Check",gffFileBasename,".txdb",sep=""))
q("no")


