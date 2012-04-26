library(seqbias)
library(Rsamtools)
library(ShortRead)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
{
  cat ("Rscript 1.R 8 1\n")
  quit("no")
}

bamFile <- "/v4scratch/sc2265/rnaseq/output/ua159/1/bwa/FASTQ001.sorted.bam"
aln <- readGappedAlignments(bamFile)
aln <- as(aln, "GRanges")

n <- as.integer(args[1]) + 1L
a <- as.integer(seq(1L,length(aln),length.out=n))
a[9] <- a[9] + 1
s.1 <- a[as.integer(args[2])]
s.2 <- a[as.integer(args[2]) + 1] - 1

# print(s.1)
# print(s.2)
# quit("no")

ref_fn <- "/v4scratch/sc2265/rnaseq/output/ua159/1/data/NC_004350.fna"
ref_f <- FaFile( ref_fn )
open.FaFile( ref_f )
ref_seqs <- scanFaIndex( ref_f )
ref_seq <- getSeq(ref_f)
I.all <- GRanges(seqnames=Rle(c("NC_004350.1"),c(2)),
                 ranges=IRanges(c(1,1),width=rep(width(ref_seq),2)),
                 strand=Rle(strand(c("+","-")),c(1,1)))
seqlengths(I.all) <- c(width(ref_seq))

sb5 <- seqbias.load( ref_fn, "/v4scratch/sc2265/rnaseq/output/ua159/1/bwa/FASTQ001-500000.yml")
bias <- seqbias.predict( sb5, I.all )

s1 <- start(aln)
s2 <- end(aln)
s3 <- strand(aln)

# foreach (i=1:length(aln)) %do% {
cvg <- rep(0,width(ref_seq))
# for (i in 1:length(aln)) {
for (i in s.1:s.2) {
  if (runValue(s3[i]) == "+") {
    cvg[s1[i]:s2[i]] <- cvg[s1[i]:s2[i]] + 1/bias[[1]][s1[i]]
  } else {
    cvg[s1[i]:s2[i]] <- cvg[s1[i]:s2[i]] + 1/bias[[2]][s2[i]]
  }
}
cvgFile <- paste("/v4scratch/sc2265/rnaseq/output/ua159/1/bwa/FASTQ001.cvg",args[2],sep=".")
save(cvg,file=cvgFile)
close.FaFile( ref_f )
