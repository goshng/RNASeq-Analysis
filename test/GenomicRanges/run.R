library(GenomicRanges)

# We load a ParseRNAseq output file first. We could load a BED-format file later.
x <- read.table("FASTQ001.parsernaseq1", head=FALSE)
number.tx <- length(rownames(x))
tx.state <- sub("_[[:digit:]]*", "", substring(x$V9,6))
tx.state <- sub("C", "", tx.state)
gr <- GRanges( seqnames = Rle("chr1", number.tx),
               ranges =
                 IRanges(start=x$V4,
                         end=x$V5,
                         names=sprintf("tx%04d",seq(number.tx))),
               strand = Rle(strand("*"),number.tx),
               score=x$V6,
               state=as.integer(tx.state)
             )

# Find the number of transcripts
cat("The number of transcripts predicted:", length(gr), "\n")

# Plot or summarize transcript lengths distribution
cat("The summary of transcript lengths\n")
show(summary(width(gr)))

# Load the pileup file as Rle or a vector object
# rm(x)
x <- read.table("FASTQ001.parsernaseq.pileup", head=FALSE)
pileup.v <- Rle(x$V2)
seqlengths(gr) <- length(pileup.v)
# grl <- split(gr)

# Compute average pilup values of the transcripts
seqselect(pileup.v, ranges(gr)[1])
pileup.x <- lapply(ranges(gr), function(x) seqselect(pileup.v, x))
pileup.y <- mapply(function(x,y) sum(x)/y, x=pileup.x, y=width(gr))
cor(pileup.y,elementMetadata(gr)[,"score"])
# The correlation was 1, so the score is the division of sum of pileup values by
# the length of a transcript.

# Find genes covered by the transcripts.
x <- read.table("feature-genome.out-geneonly", head=FALSE)
number.tx <- length(rownames(x))
grGenes <- GRanges( seqnames = Rle("chr1", number.tx),
               ranges =
                 IRanges(start=x$V2 + 1,
                         end=x$V3,
                         names=x$V4),
               strand = Rle(strand(x$V6))
             )
seqlengths(grGenes) <- length(pileup.v)

mtch <- findOverlaps(grGenes,gr)
# mtch <- grGenes %in% gr

# Plot or summarize number of genes in the transcripts
stem(table(matchMatrix(mtch)[,2]))

# Find the expressed genes
mtch2 <- findOverlaps(grGenes,gr[elementMetadata(gr)$state > 1])
length(matchMatrix(mtch2)[,1])

# Find transcripts without annotated genes as putative non-coding elements
smutansData.txGenes <- grGenes
smutansData.tx <- gr
smutansData.txPileup <- pileup.v
grTx <- smutansData.tx[elementMetadata(smutansData.tx)$state > 1]

grNoncoding <- gaps(smutansData.txGenes)
mtchNoncoding <- findOverlaps(grNoncoding,grTx)
                              
matchMatrix(mtchNoncoding)[,2]
stem(table(matchMatrix(mtchNoncoding)[,2]))


# Find transcripts with genes with conflicted strands
# Get plus strand genes.
grGenesPlus <- smutansData.txGenes[strand(smutansData.txGenes)=='+']
grGenesMnus <- smutansData.txGenes[strand(smutansData.txGenes)=='-']
mtchPlus <- findOverlaps(grGenesPlus,grTx)
mtchMnus <- findOverlaps(grGenesMnus,grTx)
txWithConflictedGenes <- intersect(matchMatrix(mtchPlus)[,2],matchMatrix(mtchMnus)[,2])

# Find the distribution of lengths of expressed genes
mtch2 <- findOverlaps(grGenes,gr[elementMetadata(gr)$state > 1])
length(matchMatrix(mtch2)[,1])
# 1. Grep genes with the expressed genes
# 2. Legnth of the genes
width(grGenes[matchMatrix(mtch2)[,1]])

# Compute average expression levels for genes 
# and correlate two gene expression levels that are adjacent in the same
# transcripts.
# 1. make pairs of genes within the same transcript.
# seqselect(smutansData.txPileup, ranges(smutansData.txGenes)[1])
pileup.x <- lapply(ranges(smutansData.txGenes), function(x) seqselect(smutansData.txPileup, x))
pileup.y <- mapply(function(x,y) sum(x)/y, x=pileup.x, y=width(smutansData.txGenes))
elementMetadata(smutansData.txGenes)["score"] <- pileup.y
mtch3 <- findOverlaps(smutansData.txGenes,grTx2)
x <- rle(matchMatrix(mtch3)[,"subject"])
comp.exp.two.genes <- c()
for (i in x$values[x$lengths > 1]) {
  y <- combn(matchMatrix(mtch3)[,"query"][matchMatrix(mtch3)[,"subject"]==i],2)
  for (j in 1:dim(y)[2]) {
    comp.exp.two.genes <- cbind(comp.exp.two.genes,as.matrix(elementMetadata(smutansData.txGenes)["score"][y[,j],]))
  }
}
plot(log10(comp.exp.two.genes[1,]), log10(comp.exp.two.genes[2,]))
# 2. make pairs of genes across adjacent transcript.
comp.exp.two.genes.between.transcript <- c()
for (i in 1:length(x$values)) {
  if (i == length(x$values)) {
    break
  }
  y1 <- matchMatrix(mtch3)[,"query"][matchMatrix(mtch3)[,"subject"]==x$values[i]]
  i2 <- i + 1
  y2 <- matchMatrix(mtch3)[,"query"][matchMatrix(mtch3)[,"subject"]==x$values[i2]]
  l <- as.matrix(do.call(expand.grid, list(y1,y2)))
  for (j in 1:dim(l)[1]) {
    comp.exp.two.genes.between.transcript <- cbind(comp.exp.two.genes.between.transcript,as.matrix(elementMetadata(smutansData.txGenes)["score"][l[j,],]))
  }
}
plot(log10(comp.exp.two.genes.between.transcript[1,]), log10(comp.exp.two.genes.between.transcript[2,]))
cor(log10(comp.exp.two.genes.between.transcript[1,]), log10(comp.exp.two.genes.between.transcript[2,]))

l <- list(a = 1:2, b = 3:4)



