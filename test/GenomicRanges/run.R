library(GenomicRanges)

# We load a ParseRNAseq output file first. We could load a BED-format file later.
x <- read.table("FASTQ001.parsernaseq1", head=FALSE)
number.tx <- length(rownames(x))
gr <- GRanges( seqnames = Rle("chr1", number.tx),
               ranges =
                 IRanges(start=x$V4,
                         end=x$V5,
                         names=sprintf("tx%04d",seq(number.tx))),
               strand = Rle(strand("*"),number.tx),
               score=x$V6,
               state=substring(x$V9,6)
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

# Find transcripts without annotated genes as putative non-coding elements
grNoncoding <- gaps(grGenes)
mtchNoncoding <- findOverlaps(grNoncoding,gr)
stem(table(matchMatrix(mtchNoncoding)[,2]))

# Find transcripts with genes with conflicted strands
# Get plus strand genes.
grGenesPlus <- grGenes[strand(grGenes)=='+']
grGenesMnus <- grGenes[strand(grGenes)=='-']
mtchPlus <- findOverlaps(grGenesPlus,gr)
mtchMnus <- findOverlaps(grGenesMnus,gr)
txWithConflictedGenes <- intersect(matchMatrix(mtchPlus)[,2],matchMatrix(mtchMnus)[,2])

