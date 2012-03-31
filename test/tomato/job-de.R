library(easyRNASeq)
library(rtracklayer)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1)
{
  cat ("Rscript job-de.R 1\n")
  quit("yes")
}
gene.range3 <- import.gff3("/v4scratch/sc2265/rnaseq/output/tomato/1/data/ITAG2.3_gene_models.gff3.txt.no-negative")
gene.range2 <- gene.range3[gene.range3$type=="exon",]
# gene.range2 <- gene.range1[grep("SMU[rt]", gene.range1$locus_tag, invert=TRUE),]

bam.file <- paste(args[1])
cl.file <- paste(args[1], "cl", sep=".")
indexFile <- indexBam(bam.file)

aln <- readAligned(bam.file,type="BAM")
aln <- aln[!is.na(position(aln))]

cl.sum <- rep(0,times=length(gene.range2$ranges))
for (i in unique(gene.range2$space)) {
  cfilt <- chromosomeFilter(paste("^",as.character(i),"$",sep=""))
  aln2 <- aln[cfilt(aln)]
  if (length(aln2) > 0) {
    alnIR <- IRanges(start=position(aln2),width=width(aln2))
    geneIR <- gene.range2[gene.range2$space==as.character(i), ]$ranges
    cl <- countOverlaps(geneIR,alnIR)
    cl.sum[gene.range2$space==as.character(i)] <- cl.sum[gene.range2$space==as.character(i)] + cl
  }
}
cl <- cl.sum
save(cl,file=cl.file)

