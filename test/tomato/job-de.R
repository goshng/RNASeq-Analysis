library(easyRNASeq)
library(rtracklayer)
args <- commandArgs(trailingOnly = TRUE)
#args <- "BWA-00.bam"
if (length(args) != 1)
{
  cat ("Rscript job-de.R 1\n")
  quit("yes")
}
#gene.range3 <- import.gff3("/v4scratch/sc2265/rnaseq/output/tomato/1/data/ITAG2.3_gene_models.gff3.txt.no-negative")
load("gene.range3.Rd")
gene.range2 <- gene.range3[gene.range3$type=="gene",]
# gene.range2 <- gene.range1[grep("SMU[rt]", gene.range1$locus_tag, invert=TRUE),]

bam.file <- paste(args[1])
cl.file <- paste(args[1], "cl", sep=".")
indexFile <- indexBam(bam.file)

aln <- readAligned(bam.file,type="BAM")
aln <- aln[!is.na(position(aln))]

cl.sum <- rep(0,times=length(gene.range2$ranges))

for (i in unique(gene.range2$space)) {
  for (s in c('+','-')) {
    cfilt <- chromosomeFilter(paste("^",as.character(i),"$",sep=""))
    # aln2 <- aln[cfilt(aln)]
    if (s == '+') {
      sa <- '-'
    } else {
      sa <- '+'
    }
    aln2 <- aln[cfilt(aln) & strand(aln)==sa]
    selected.gene.index <- gene.range2$space==as.character(i) & gene.range2$strand==s
    if (length(aln2) > 0) {
      alnIR <- IRanges(start=position(aln2),width=width(aln2))
      geneIR <- gene.range2[selected.gene.index, ]$ranges
      cl <- countOverlaps(geneIR,alnIR)
      cl.sum[selected.gene.index] <- cl.sum[selected.gene.index] + cl
    }
  }
}
cl <- cl.sum
save(cl,file=cl.file)

# aln.plus <- aln[strand(aln)=='+']

# gene.plus <- gene.range2[gene.range2$strand=='+',]


# bash job-de BWA-00.bam
# bash job-de2 /v4scratch/sc2265/rnaseq/output/tomato/1/bwa/FASTQ001.cl BWA-
