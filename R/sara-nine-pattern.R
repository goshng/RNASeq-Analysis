args <- commandArgs(trailingOnly = TRUE)                                           
if (length(args) != 2)
{
  cat ("Rscript R/find-il-position.R IL1-1_SNP\n")
  quit("no")
}
#args<- c("email/to/sara-palmer/092412/ua159ph-annotation.csv", "email/to/sara-palmer/092412/smu21ph-annotation.csv")

x <- read.csv(args[1])
y <- read.csv(args[2])

# Two-Fold Changes
for (i in c("plus","minus","zero")) {
  for (j in c("plus","minus","zero")) {
    m <- x$padj < 0.1 & y$padj < 0.1
    if (i == "plus" & j == "plus") {
      m <- m & x$log2FoldChange > 1 & y$log2FoldChange > 1
    } else if (i == "plus" & j == "minus") {
      m <- m & x$log2FoldChange > 1 & y$log2FoldChange < -1
    } else if (i == "plus" & j == "zero") {
      m <- m & x$log2FoldChange > 1 & abs(y$log2FoldChange) < 1
    } else if (i == "minus" & j == "plus") {
      m <- m & x$log2FoldChange < -1 & y$log2FoldChange > 1
    } else if (i == "minus" & j == "minus") {
      m <- m & x$log2FoldChange < -1 & y$log2FoldChange < -1
    } else if (i == "minus" & j == "zero") {
      m <- m & x$log2FoldChange < -1 & abs(y$log2FoldChange) < 1
    } else if (i == "zero" & j == "plus") {
      m <- m & abs(x$log2FoldChange) < 1 & y$log2FoldChange > 1
    } else if (i == "zero" & j == "minus") {
      m <- m & abs(x$log2FoldChange) < 1 & y$log2FoldChange < -1
    } else if (i == "zero" & j == "zero") {
      m <- m & abs(x$log2FoldChange) < 1 & abs(y$log2FoldChange) < 1
    }
    m[is.na(m)] <- FALSE
    d <- data.frame(id=x[m,"id"],
                    x.baseMean=x[m,"baseMean"],
                    x.baseMeanA=x[m,"baseMeanA"],
                    x.baseMeanB=x[m,"baseMeanB"],
                    x.foldChange=x[m,"foldChange"],
                    x.log2FoldChange=x[m,"log2FoldChange"],
                    x.padj=x[m,"padj"],
                    y.baseMean=y[m,"baseMean"],
                    y.baseMeanA=y[m,"baseMeanA"],
                    y.baseMeanB=y[m,"baseMeanB"],
                    y.foldChange=y[m,"foldChange"],
                    y.log2FoldChange=y[m,"log2FoldChange"],
                    y.padj=y[m,"padj"],
                    core=y[m,"core"],
                    oralgen.annotation=y[m,"oralgen.annotation"],
                    genename=y[m,"genename"],
                    some.annotation=y[m,"some.annotation"])
    if (nrow(d) > 0) {
      d.file <- paste(dirname(args[1]),"/",sub("-annotation.csv","",basename(args[1])),"-",i,"-",sub("-annotation.csv","",basename(args[2])),"-",j,".csv",sep="")
      write.csv(d,
                file=d.file,
                row.names=FALSE)
    }
  }
}
