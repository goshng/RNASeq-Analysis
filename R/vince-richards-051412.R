library(smutans)
# REFGENOMEID:NC_007432
# REFGENOMEID:NC_004368
# REFGENOMEID:AEXT01
# REFGENOMEID:single
# count-hcmb-sag56-to-human.index
# 

# 051412, 051512
# for (c.j in c("sag56human","sag56cow")) {
# 051612
# for (c.j in c("56n61vs56n2")) {

qval <- 0.05
baseDirectory <- "output/email/to/vince-richards/052912"
for (c.i in c("milk","broth")) {
for (c.j in c("82n61vs56n2")) {
# for (c.k in c("NC_007432","NC_004368","AEXT01","core2")) {
for (c.k in c("core2")) {
file.1 <- paste("output/agalactiae/1/bwa/count-",c.k,".cds.txt", sep="")
file.out.1 <- paste(baseDirectory,"/",c.k,"-",c.i,"-",c.j,"-clust.pdf", sep="")
file.out.2 <- paste(baseDirectory,"/",c.k,"-",c.i,"-",c.j,".csv", sep="")
file.out.3 <- paste(baseDirectory,"/",c.k,"-",c.i,"-",c.j,".ps", sep="")

if (c.i == "milk") {
agGenes <- 
  readSmutans(countsFile=file.1,
              indexFile=paste("output/agalactiae/1/run-analysis/count-",c.j,".index",sep=""),
              condition = c("hm", "cm"),
              firstFactorLabel = c("hm"),
              secondFactorLabel = c("hm", "cm"),
              name="S. agalactiae human and cow with milk fixed")
ag <- newSmutans( agGenes , title="Human and Cow with milk fixed" )
pdf(file.out.1)
smutans.de2Clust( ag )
dev.off()
ag <- smutans.de2( ag, type="hm", condA="hm", condB="cm")
smutans.de2List( ag, file=file.out.2 )
smutans.plotDiffExp(ag,file=file.out.3,qval=qval)

} else {
agGenes <- 
  readSmutans(countsFile=file.1,
              indexFile=paste("output/agalactiae/1/run-analysis/count-",c.j,".index",sep=""),
              condition = c("hb", "cb"),
              firstFactorLabel = c("hb"),
              secondFactorLabel = c("hb", "cb"),
              name="S. agalactiae human and cow with broth fixed")
ag <- newSmutans( agGenes , title="Human and Cow with broth fixed" )
pdf(file.out.1)
smutans.de2Clust( ag )
dev.off()
ag <- smutans.de2( ag, type="hb", condA="hb", condB="cb")
smutans.de2List( ag, file=file.out.2 )
smutans.plotDiffExp(ag,file=file.out.3,qval=qval)

}

cowOverExpresion <- ag@res[ag@res$padj < qval & ag@res$foldChange > 1,]
humanOverExpresion <- ag@res[ag@res$padj < qval & ag@res$foldChange < 1,]
print(paste("#Count", c.k, c.i, c.j, length(cowOverExpresion$id), length(humanOverExpresion$id)))

}
}
}

