library(smutans)

qval <- 0.05
baseDirectory <- "output/email/to/vince-richards/051412"

c.k <- "NC_004368"
file.1 <- paste("output/agalactiae/1/bwa/count-",c.k,".cds.txt", sep="")
file.out.1 <- paste(baseDirectory,"/",c.k,"-sag56-A-milk-B-broth-clust.pdf", sep="")
file.out.2 <- paste(baseDirectory,"/",c.k,"-sag56-A-milk-B-broth.csv", sep="")
file.out.3 <- paste(baseDirectory,"/",c.k,"-sag56-A-milk-B-broth.ps", sep="")

agGenes <- 
  readSmutans(countsFile=file.1,
              indexFile=paste("output/agalactiae/1/run-analysis/count.txt.index",sep=""),
              condition = c("sag56hm", "sag56hb"),
              firstFactorLabel = c("sag56hm"),
              secondFactorLabel = c("sag56hm", "sag56hb"),
              name="S. agalactiae NEM316 milk and broth")
ag <- newSmutans( agGenes , title="Milk and broth" )
pdf(file.out.1)
smutans.de2Clust( ag )
dev.off()
ag <- smutans.de2( ag, type="sag56hm", condA="sag56hm", condB="sag56hb")
smutans.de2List( ag, file=file.out.2 )
smutans.plotDiffExp(ag,file=file.out.3,qval=qval)
