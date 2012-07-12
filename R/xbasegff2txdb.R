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
sm.source <- sm.gff[sm.gff$type=="DNA",]

transcripts <- 
  data.frame( tx_id=seq(length(sm.gene$locus_tag)), # integer vector
              tx_name=sm.gene$locus_tag,            # string vector or factor
              tx_chrom=sm.gene$space,               # string vector or factor
              tx_strand=sm.gene$strand,             # strand vector
              tx_start=start(sm.gene),              # integer vector
              tx_end=end(sm.gene) )                 # integer vector

# Each transcript contains an exon. The exon may contain a CDS. 
# Compare CDS and gene locus tags
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

# 
splicings <- 
  data.frame( tx_id=seq(length(sm.gene$locus_tag)),     # integer A
              exon_rank=seq(length(sm.gene$locus_tag)), # integer B (unique A-B)
              exon_name=sm.gene$locus_tag,              # character
              exon_start=start(sm.gene),                # integer 
              exon_end=end(sm.gene),                    # integer
              cds_start=grCDS.start,                    # integer
              cds_end=grCDS.end )                       # integer

chrominfo <-
  data.frame( chrom=names(sm.gene),    # character
              length=end(sm.source),   # integer
              is_circular=FALSE ) 

txdb <- makeTranscriptDb (transcripts, splicings, chrominfo=chrominfo)
saveFeatures(txdb,file=paste(gffFileBasename,"txdb",sep="."))
print(paste("Check ",gffFileBasename,".txdb",sep=""))
q("no")
