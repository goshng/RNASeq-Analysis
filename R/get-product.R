# similar-gene contains, smu21 gene, UA159 gene, and coverage.
# I wish to find gene annotation or gene product of UA159 using NC_004350.gff.
# CDS has the gene product. CDS's parent points to gff's gene. CDS's does not
# have a locus tag. I use locus tag from sm.gene to find parent of CDS, and use
# the CDS to find gene product.
# input files: NC_004350.gff, similar-gene
# output: similar-gene.txt
library(GenomicFeatures)                                                        
library(rtracklayer)

sm.gff <- import.gff3("NC_004350.gff")
sm.gene <- sm.gff[sm.gff$type=="gene",]
sm.CDS <- sm.gff[sm.gff$type=="CDS",]
gene <- read.table("similar-gene")

gene.products <- c()
if (length(unlist(sm.CDS$Parent)) == nrow(sm.CDS)) {
  for (i in gene$V2) {
    print (i)
    gene.id <- sm.gene[sm.gene$locus_tag == i,"ID"]$ID
    gene.product <- sm.CDS[unlist(sm.CDS$Parent)==gene.id,"product"]$product
    gene.products <- c(gene.products, gene.product)
  }
} else {
  print("Error: Multple CDS in genes")
}

x <- data.frame(gene, product=gene.products)
colnames(x) <- c("Smu21","UA159","Coverage","Product")
write.table(x,file="similar-gene.txt",sep="\t",quote=FALSE,row.names=FALSE)
