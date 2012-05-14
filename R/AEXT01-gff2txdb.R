library(DESeq)
library(ShortRead)
library(rtracklayer)
library(GenomicRanges)
library(VariantAnnotation)
library(GenomicFeatures)
gffFile <- "data/smutans/Streptococcus_agalactiae_FSL_S3-026/AEXT01.gff"
sm.gff <- import.gff3(gffFile)
#sm.source <- sm.gff[1,] #sm.gff=="region",]
# For AEXT01 contig genome
sm.source <- sm.gff[grep("^AEXT01",sm.gff$ID),]
sm.gene <- sm.gff[sm.gff$type=="gene",]
sm.exon <- sm.gff[sm.gff$type=="exon",]
sm.CDS <- sm.gff[sm.gff$type=="CDS",]
sm.CDS <- sm.gff[sm.gff$type=="mRNA",]

transcripts <- 
  data.frame( tx_id=seq(length(sm.gene$locus_tag)),
              tx_name=sm.gene$locus_tag,
              tx_chrom=sm.gene$space, 
              tx_strand=sm.gene$strand,
              tx_start=start(sm.gene),
              tx_end=end(sm.gene) )

# Compare CDS and gene locus tags
x <- c()
y <- c()
# AEXT only
#z <- sub(".t01","",unlist(sm.CDS$Parent))
z <- unlist(sm.CDS$Parent)
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

# 
splicings <- 
  data.frame( tx_id=seq(length(sm.gene$locus_tag)),
              exon_rank=seq(length(sm.gene$locus_tag)),
              exon_start=start(sm.gene),
              exon_end=end(sm.gene),
              exon_name=sm.gene$locus_tag,
              cds_start=grCDS.start,
              cds_end=grCDS.end )

# For AEXT01 contig genome
#              cds_start=start(sm.gene),
#              cds_end=end(sm.gene),

chrominfo <-
  data.frame( chrom=names(sm.gene),
              length=end(sm.source),
              is_circular=FALSE )

TxDb.Smutans.UA159.uid57947.knownGene <- makeTranscriptDb (transcripts, splicings, chrominfo=chrominfo)
txdb <- TxDb.Smutans.UA159.uid57947.knownGene 
saveFeatures(txdb,file="data/smutans/Streptococcus_agalactiae_FSL_S3-026/AEXT01.txdb")
print("Check data/smutans/Streptococcus_agalactiae_A909_uid57935/NC_007432.txdb")

