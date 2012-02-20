sample(194:1552,size=2,replace=TRUE)
read.pos <- c(mapply(function(x,y) sample(x:y,size=10,replace=TRUE),start(geneIR),end(geneIR)))
alnIR <- IRanges(start=position(aln2),width=width(aln2))
length(s.mutans)
subseq(s.mutans,10,12)
reverseComplement(subseq(s.mutans,10,12))

library(easyRNASeq)
library(rtracklayer)
gene.range3 <- import.gff3("/v4scratch/sc2265/rnaseq/output/ua159/1/data/NC_004350.gff")
gene.range2 <- gene.range3[gene.range3$type=="gene",]
geneIR <- gene.range2$ranges
s.mutans <- read.DNAStringSet("/v4scratch/sc2265/rnaseq/output/ua159/1/data/NC_004350.fna")
s.mutans <- s.mutans[[1]]

# read.pos <- c(mapply(function(x,y) sample(x:y,size=3,replace=TRUE),start(geneIR[1:2]),end(geneIR[1:2])))
read.pos <- c(mapply(function(x,y) sample(x:y,size=10,replace=TRUE),start(geneIR),end(geneIR)))
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

read.start.end.strand <- mapply(read.loc, read.pos)

read.GR <- GRanges( seqnames = Rle("chr1", length(read.pos)),
                    ranges =
                      IRanges(start=unlist(read.start.end.strand[1,]),
                              end=unlist(read.start.end.strand[2,]),
                              names=sprintf("read%09d",seq(length(read.pos)))),
                     strand = unlist(read.start.end.strand[3,])
                   )
# Now, create a count table and fastQ file.
# @HWI-ST397:202:C08BKACXX:4:1101:1704:2115 1:N:0:ACAGTG
# NTTAAAGCTCGTTTCGCTCGCAATAAGGAAAATAATTCTGAATTTACAGAGCTAAAAAAGATTTGTCTATAAGCAAGTTAATATTCCAAAAGTATTCAGA
# +
# #1=DDDFFHHHHHJJJJJJJJJJJJJJJGGIJJJIIJIJJJJJJJJJJIJGIIJJJIJHFEFFFFFFFEEEEEDDCCDCDDEEEFEEEDDDDACDEEEDD
read.extract <- function (x,y,z,w) {
  s <- subseq(s.mutans,x,y)
  if (z == "-") {
    s <- reverseComplement(s)
  }
  s.quality <- paste((rep("J",times=length(s))),collapse="")
  oneRead <- sprintf("@HWI-ST397:%09d 1:N:0:ACAGTG\n%s\n+\n%s\n", w, s, s.quality)
  cat(oneRead, file="x.fq", append=TRUE)
}

read.simulated <- mapply(read.extract, 
                         unlist(read.start.end.strand[1,]),
                         unlist(read.start.end.strand[2,]), 
                         unlist(read.start.end.strand[3,]),
                         seq(length(unlist(read.start.end.strand[3,]))))
alnIR <- ranges(read.GR)
cl <- countOverlaps(geneIR,alnIR)
save(cl,file="x.cl")

######################################################################

bam.file <- paste(args[1])
cl.file <- paste(args[1], "cl", sep=".")
indexFile <- indexBam(bam.file)
aln <- readAligned(bam.file,type="BAM")
cfilt <- chromosomeFilter('gi|24378532|ref|NC_004350.1|')
aln2 <- aln[cfilt(aln)]
alnIR <- IRanges(start=position(aln2),width=width(aln2))
geneIR <- gene.range2$ranges
cl <- countOverlaps(geneIR,alnIR)
save(cl,file=cl.file)
