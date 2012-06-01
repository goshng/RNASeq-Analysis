library(smutans)

qval <- 0.05
baseDirectory <- "output/email/to/vince-richards/052912"

c.i <- "milk"
c.j <- "std"

for (c.k in c("core2")) {

file.1 <- paste("output/agalactiae/1/bwa/count-",c.k,".cds.txt", sep="")
file.out.1 <- paste(baseDirectory,"/",c.k,"-",c.i,"-",c.j,"-clust.pdf", sep="")
file.out.2 <- paste(baseDirectory,"/",c.k,"-",c.i,"-",c.j,".csv", sep="")
file.out.3 <- paste(baseDirectory,"/",c.k,"-",c.i,"-",c.j,".ps", sep="")

agGenes <- 
  readSmutans(countsFile=file.1,
              indexFile=paste("output/agalactiae/1/run-analysis/count.txt.index",sep=""),
              condition = c("sag82hm", "sag2cm"),
              firstFactorLabel = c("sag82hm"),
              secondFactorLabel = c("sag82hm", "sag2cm"),
              name="S. agalactiae sag82 human and sag2 cow with milk fixed")
ag <- newSmutans( agGenes , title="Human and Cow with milk fixed" )
pdf(file.out.1)
smutans.de2Clust( ag )
dev.off()
ag <- smutans.de2( ag, type="sag82hm", condA="sag82hm", condB="sag2cm")
smutans.de2List( ag, file=file.out.2 )
smutans.plotDiffExp(ag,file=file.out.3,qval=qval)

cowOverExpresion <- ag@res[ag@res$padj < qval & ag@res$foldChange > 1,]
humanOverExpresion <- ag@res[ag@res$padj < qval & ag@res$foldChange < 1,]
print(paste("#Count", c.k, c.i, c.j, length(cowOverExpresion$id), length(humanOverExpresion$id)))

}

