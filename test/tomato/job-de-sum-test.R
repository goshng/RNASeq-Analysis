#library(rtracklayer)
#args <- commandArgs(trailingOnly = TRUE)

#args <- c("/v4scratch/sc2265/rnaseq/output/tomato/1/bwa", "index.txt")
#gene.range3 <- import.gff3("/v4scratch/sc2265/rnaseq/output/tomato/1/data/ITAG2.3_gene_models.gff3.txt.no-negative")
#gene.range2 <- gene.range3[gene.range3$type=="exon",]


#################################################
# gene names
#x <- sub("mRNA:","",gene.range2$Parent)
#genesrle <- gsub("([\\.*])","",x)

#x <- c() 
#for (i in rle(genesrle)$lengths) {
  #for (j in 1:i) {
    #x <- c(x, sprintf("E%03d", j)) 
  #}
#}
#exons <- x

#genexonid <- paste(genesrle,exons,sep=":")
#################################################
# count.table <- data.frame(gene=gene.range2$locus_tag)
#count.table <- data.frame(gene=genexonid)

fastQIndex <- scan(args[2])
for (i in fastQIndex) {
  cl.file <- sprintf("%s/FASTQ%03d.cl", args[1], i)
  load(cl.file)
  count.table <- data.frame(count.table,cl)
  colnames(count.table)[ncol(count.table)] <- paste("X",i,sep="")
}
colnames(count.table) <- sub("X","",colnames(count.table))
count.table.file <- sprintf("%s/count.txt", args[1])
write.table(count.table,file=count.table.file,quote=FALSE,sep="\t",row.names=FALSE)
