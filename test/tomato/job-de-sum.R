library(rtracklayer)
#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) != 2)
#{
#  cat ("Rscript job-de-sum.R output/ua159/1/subread fastq_index_file\n")
#  quit("yes")
#}
args <- c("/v4scratch/sc2265/rnaseq/output/tomato/1/bwa", "index.txt")
# gene.range3 <- import.gff3("/v4scratch/sc2265/rnaseq/output/tomato/1/data/ITAG2.3_gene_models.gff3.txt.no-negative")
load("gene.range3.Rd")
gene.range2 <- gene.range3[gene.range3$type=="gene",]

#################################################
# gene and exon combination names
if (FALSE) {
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
#################################################
} else {
#################################################
# gene names
x <- sub("gene:","",gene.range2$ID)
genesrle <- gsub("\\.[0-9]","",x)
genesrle <- unlist(genesrle)
genexonid <- genesrle
#################################################
}

# count.table <- data.frame(gene=gene.range2$locus_tag)
count.table <- data.frame(gene=genexonid)

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
print(paste("File:", count.table.file))

