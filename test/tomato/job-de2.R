#library(easyRNASeq)
#library(rtracklayer)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
{
  cat ("Rscript job-de.R 1 BWA-\n")
  quit("yes")
}
# gene.range3 <- import.gff3("/v4scratch/sc2265/rnaseq/output/tomato/1/data/SMU21.gff")
# gene.range1 <- gene.range3[gene.range3$type=="gene",]
# gene.range2 <- gene.range1[grep("SMU[rt]", gene.range1$locus_tag, invert=TRUE),]

# cl.sum <- rep(0,times=length(gene.range2$ranges))

cl.pattern <- paste(args[2],"[[:digit:]]+.bam.cl",sep="")
cl.files <- list.files(pattern=cl.pattern)
load(cl.files[1])
cl.sum <- rep(0,times=length(cl))
for (i in cl.files) {
  load(i)
  cl.sum <- cl.sum + cl
}
cl <- cl.sum
save(cl,file=args[1])
print(paste("file:",args[1]))
