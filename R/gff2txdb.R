library(rtracklayer)
library(GenomicRanges)
library(VariantAnnotation)
library(GenomicFeatures)
gffFile <- "/Volumes/Elements/Documents/Projects/mauve/bacteria/Streptococcus_mutans_UA159_uid57947/NC_004350.gff"
sm.gff <- import.gff3(gffFile)
sm.source <- sm.gff[sm.gff$type=="source",]
sm.gene <- sm.gff[sm.gff$type=="gene",]
sm.exon <- sm.gff[sm.gff$type=="exon",]
sm.CDS <- sm.gff[sm.gff$type=="CDS",]

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
for (i in sm.gene$locus_tag) {
  if (sum(sm.CDS$locus_tag == i) > 0) {
    x <- c(x,start(sm.CDS)[sm.CDS$locus_tag == i])
    y <- c(y,end(sm.CDS)[sm.CDS$locus_tag == i])
  } else {
    x <- c(x,NA)
    y <- c(y,NA)
  }
}
grCDS.start <- x
grCDS.end <- y
rm(i,x,y)

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
              is_circular=TRUE )

TxDb.Smutans.UA159.uid57947.knownGene <- makeTranscriptDb (transcripts, splicings, chrominfo=chrominfo)
txdb <- TxDb.Smutans.UA159.uid57947.knownGene 

cds(txdb, columns="tx_name")
cds(txdb, columns="exon_name")

# Read short reads
library(ShortRead)
aln <- readAligned(bam.file,type="BAM")
ShortRead::readAligned                   # seems to be good for exploring reads
# AlignedRead object? from which package?
GenomicRanges::summarizeOverlaps         # GappedAlignments is part of GenomicRanges
Rsamtools::readBamGappedAlignments       # seems to be more efficient than readAligned
GenomicRanges::readGappedAlignment # no such function exists.
                                        

# Test of reads
  group_id <- c("A", "B", "C", "C", "D", "D", "E", "F", "G", "H", "H")
  features <- GRanges(
      seqnames = Rle(c("chr1", "chr2", "chr1", "chr1", "chr2", "chr2", 
          "chr1", "chr1", "chr2", "chr1", "chr1")),
      strand = strand(rep("+", length(group_id))),
      ranges = IRanges(
          start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000, 3000, 
              5000, 5400),
          width=c(500, 900, 500, 300, 600, 300, 500, 900, 500, 500, 500)),
     DataFrame(group_id)
  )

  reads <- GappedAlignments(
      names = c("a","b","c","d","e","f","g"),
      seqnames = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
      pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
      cigar = c("500M", "100M", "300M", "500M", "300M", 
          "50M200N50M", "50M150N50M"),
      strand = strand(rep("+", 7)))

solap <- summarizeOverlaps(features, reads)



