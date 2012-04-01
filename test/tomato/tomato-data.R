# Load libraries
library(DEXSeq)
library(rtracklayer)
library(multicore)
library("pasilla")

gene.range3 <- import.gff3("ITAG2.3_gene_models.gff3.txt.no-negative")
# gene.range3 <- import.gff3("test.gff")
gene.range2 <- gene.range3[gene.range3$type=="exon",]
sa <- data.frame(condition=factor(c(rep("wt",1),rep("il",1))),replicate=c(1,1),type=c("end","end"))
rownames(sa) <- c("wt","il")
sa <- data.frame(condition=factor(c(rep("wt",2),rep("il",2))),replicate=c(1,2,1,2),type=c("end","end","end","end"))
rownames(sa) <- c("wt1","wt2","il1","il2")

exoninfo <- data.frame(chr=gene.range2$space, start=start(gene.range2), end=end(gene.range2), strand=strand(gene.range2)) 
x <- sub("mRNA:","",gene.range2$Parent)
genesrle <- gsub("\\.[0-9]+","",x)
genesrle <- unlist(genesrle)

y <- rep("",times=length(genesrle)) 
for (i in unique(genesrle)) {
  x <- c()
  for (j in 1:sum(genesrle == i)) {
    x <- c(x, sprintf("E%03d", j))
  }
  y[which(genesrle == i)] <- x
}
exons <- y

genexonid <- paste(genesrle,exons,sep=":")

# No alternative splicing
transcripts <- rep("A",length(genexonid))

countsFile <- "/v4scratch/sc2265/rnaseq/output/tomato/1/bwa/count2.txt"
countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- genexonid
#countsTable <- countsTable[,-1]

ecs <- newExonCountSet(countData=countsTable, design=sa, geneIDs=genesrle, exonIDs=exons, exonIntervals=exoninfo, transcripts=transcripts)


pasillaExons <- ecs
pasillaExons <- estimateSizeFactors(pasillaExons)
pasillaExons <- estimateDispersions(pasillaExons)
pasillaExons <- fitDispersionFunction(pasillaExons)
pasillaExons <- testForDEU(pasillaExons)

