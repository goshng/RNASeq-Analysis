############################################################################
# Create a test fastq files by sampling 10 100-bp short reads for each gene.
library(easyRNASeq)
library(rtracklayer)
gff.file <- "/v4scratch/sc2265/rnaseq/output/omz/1/data/SMU109.gff"
genome.file <- "/v4scratch/sc2265/rnaseq/output/omz/1/data/SMU109-OMZ175.fa"
fasta.file <- "SMU109-gene.fa"
file.create(fasta.file)

gene.range3 <- import.gff3(gff.file)
gene.range1 <- gene.range3[gene.range3$type=="gene",]
gene.range2 <- gene.range1[grep("SMU[rt]", gene.range1$locus_tag, invert=TRUE),]

s.mutans.sequence <- read.DNAStringSet(genome.file)
chrom.list <- names(s.mutans.sequence)

for (i in seq(nrow(gene.range2))) {
  j <- as.character(gene.range2[i,]$space)
  x <- start(gene.range2[i,]$ranges)
  y <- end(gene.range2[i,]$ranges)
  s.mutans <- s.mutans.sequence[[j]]
  gene.seq <- subseq(s.mutans,start=x,end=y)
  if (gene.range2[i,]$strand == '-') {
    gene.seq <- reverseComplement(gene.seq)
  }
  oneRead <- sprintf(">%s\n%s\n",gene.range2[i,]$locus_tag,toString(gene.seq))
  cat(oneRead, file=fasta.file, append=TRUE)
}
