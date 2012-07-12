library(DESeq)
library(ShortRead)
library(GenomicRanges)
library(VariantAnnotation)
library(GenomicFeatures)

library(rtracklayer)
gffFile <- "data/NC_004350.gff"
sm.gff <- import.gff3(gffFile)
sm.source <- sm.gff[1,] #sm.gff$type=="region",]

sm.gene <- sm.gff[sm.gff$type=="gene",]
sm.CDS <- sm.gff[sm.gff$type=="CDS",]
sm.source <- sm.gff[sm.gff$type=="DNA",]

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
saveFeatures(txdb,file="data/NC_004350.txdb")
print("Check data/NC_004350.txdb")
q("no")


cds(txdb, columns="tx_name")
cds(txdb, columns="exon_name")

# Read short reads
aln <- readAligned(bam.file,type="BAM")
ShortRead::readAligned                   # seems to be good for exploring reads
# AlignedRead object? from which package?
GenomicRanges::summarizeOverlaps         # GappedAlignments is part of GenomicRanges
Rsamtools::readBamGappedAlignments       # seems to be more efficient than readAligned
GenomicRanges::readGappedAlignment # This is just calling readBamGappedalignments

# Basic statistics of short reads
fq.file <- "/v4scratch/sc2265/rnaseq/output/cornell/1/data/FASTQ051.fq.gz"
bam.file <- "/v4scratch/sc2265/rnaseq/output/cornell/1/bwa/FASTQ051.sorted.bam"

fq.file <- "/Users/goshng/Documents/Projects/RNASeq-Analysis/FASTQ051.fq.gz"
bam.file <- "/Users/goshng/Documents/Projects/RNASeq-Analysis/FASTQ051.sorted.bam"
indexBam(bam.file)
# 1. Total number of short reads in fastq: zcat FASTQ051.fq.gz | wc -l
# 15166944
# scanBamWhat()
# "qname"  "flag"   "rname"  "strand" "pos"    "qwidth" "mapq"   "cigar" "mrnm"   "mpos"   "isize"  "seq"    "qual"
bam.what <- c("qname","flag","pos","mapq","cigar")
bvAll <- scanBam(bam.file, param=ScanBamParam(what=bam.what))
bvAll <- scanBam(bam.file, param=ScanBamParam(what=c("mapq")))
# 1. Number of short reads in the bam file
# > length(bvAll[[1]][["mapq"]])
# [1] 14403971
# $ samtools view FASTQ051.sorted.bam | wc -l
# 14403971
bv <- readBamGappedAlignments(bam.file,use.names=TRUE,param=ScanBamParam(what=c("mapq")))
# 2. Number of short reads mapped 
# > length(bv)
# [1] 14081608
# bv1 <- bvAll[[1]][["qname"]][!bvAll[[1]][["qname"]] %in% names(bv)]
# What are the differences?
# They are not actually aligned. For example,
# samtools view FASTQ051.sorted.bam | grep HWI-ST1085\:57\:D0E15ACXX\:8\:1101\:7072\:1957
# > table(bvAll[[1]][["flag"]])
#       0       4      16 
# 7033720  322363 7047888
#   
# Select mapped reads with mapq greater than 30
bv <- bv[elementMetadata(bv)["mapq"][,1]>30]
# 3. Number of uniquely mapped short reads
# > length(bv)
# [1] 5460982
# Take too much time to read in fastq file: fq <- readFastq(fq.file)
# Map the uniquely mapped reads on CDS regions
feature.cds <- cds(txdb, columns="exon_name")
# Remove strands information
strand(feature.cds) <- '*'
strand(bv) <- '*'
stopifnot(length(runValue(seqnames(feature.cds)))==1)
seqnames(bv) <- rep(runValue(seqnames(feature.cds)),length(bv))
olap <- summarizeOverlaps(feature.cds, bv,mode="IntersectionStrict")
# 4. Number of short reads uniquely mapped on CDS regions
# > sum(assays(olap)$counts)
# [1] 4251781
deseq <- newCountDataSet(countData=assays(olap)$counts, conditions=c("51"))

# Flattened tx
feature.tx.flat <- transcripts(txdb)
strand(feature.tx.flat) <- '*'
feature.tx.flat <- reduce(feature.tx.flat)

feature.nc <- transcripts(txdb)
strand(feature.nc) <- '*'
feature.nc <- gaps(reduce(feature.nc))
feature.nc <- feature.nc[-2:-1]
# sum(width(feature.nc))
stopifnot(length(runValue(seqnames(feature.nc)))==1)
seqnames(bv) <- rep(runValue(seqnames(feature.nc)),length(bv))
olap <- summarizeOverlaps(feature.nc,bv,mode="IntersectionStrict")
# 5. Number of short reads uniquely mapped on non-coding regions
# > sum(assays(olap)$counts)
# [1] 510054
sum(width(feature.nc)) + sum(width(feature.tx.flat)) 

# Simulate (create) short reads
genome.file <- "/v4scratch/sc2265/rnaseq/output/cornell/1/data/NC_004350.fna"
fastq.file <- "test.fq"
file.create(fastq.file)
read.GR <- GRanges()

read.loc <- function(x) { 
  s <- sample(c(TRUE,FALSE),size=1) 
  y <- x + 99
  if (s == FALSE) {
    y <- x - 99
  }
  if (y < 1) {
    y <- 1
  }
  if (y > length(s.mutans)) {
    y <- length(s.mutans)
  }
  z <- strand("+")
  if (x > y) {
    z <- x
    x <- y
    y <- z
    z <- strand("-") 
  }
  list(x,y,z)
}

read.extract <- function (x,y,z,w) {
  s <- subseq(s.mutans,x,y)
  if (z == "-") {
    s <- reverseComplement(s)
  }
  s.quality <- paste((rep("B",times=length(s))),collapse="")
  oneRead <- sprintf("@HWI-ST397:%09d 1:N:0:ACAGTG\n%s\n+\n%s\n", w, s, s.quality)
  cat(oneRead, file=fastq.file, append=TRUE)
}

s.mutans.sequence <- read.DNAStringSet(genome.file)
chrom.list <- names(s.mutans.sequence)

#i <- chrom.list[1]
#i <- chrom.list[2]

for (i in names(s.mutans.sequence)) { 
  cat (i,"\\n")

i <- "NC_004350.1"

  geneIR <- ranges(feature.cds)
  geneIR <- geneIR[width(geneIR)>1000]
  geneIR <- geneIR - 100
  if (length(geneIR) > 0) {
    read.pos <- c(mapply(function(x,y) sample(x:y,size=10,replace=TRUE),start(geneIR),end(geneIR)))
    s.mutans <- s.mutans.sequence[[i]]

    read.start.end.strand <- mapply(read.loc, read.pos)
    
    one.GR <- GRanges( seqnames = Rle(i, length(read.pos)),
                       ranges =
                         IRanges(start=unlist(read.start.end.strand[1,]),
                                 end=unlist(read.start.end.strand[2,]),
                                 names=sprintf("%s.%09d",i,seq(length(read.pos)))),
                        strand = unlist(read.start.end.strand[3,])
                     )
    read.GR <- c(read.GR,one.GR)

    read.simulated <- mapply(read.extract, 
                             unlist(read.start.end.strand[1,]),
                             unlist(read.start.end.strand[2,]), 
                             unlist(read.start.end.strand[3,]),
                             seq(length(unlist(read.start.end.strand[3,]))))
  }
}

