# This does functional category analyses. 
gobacteriatxt <- "test/RNAz/smutans.gene2go"
gocattxt <- "test/RNAz/smutans.go2ngene"

gocats <- read.table(gobacteriatxt, header=FALSE, col.names=c("gene", "go", "pval"))
gocats <- gocats[gocats$pval <= 1.0e-5,]
descrips <- read.table(gocattxt, sep="\t", header=FALSE, col.names=c("go", "count", "descrip"), quote="")
gocats$go <- as.character(gocats$go)
descrips$go <- as.character(descrips$go)

run.mw <- function(x, gocats, gocats.col="go", descrips=NULL, noisy=FALSE) {
  results <- data.frame() 
  for (gocat in unique(gocats[,gocats.col])) {
    genes <- gocats[gocats[,gocats.col]==gocat, "gene"]
    if (length(genes) >= 10) {
      f <- is.element(x$gene, genes)
      notf <- !f
      if (sum(f) >= 10) {
        wt <- wilcox.test(x[f,"score"], x[!f,"score"], alternative="greater")
        if (noisy) cat(gocat, sum(f), wt$p.value, sep="\t")
        if (!is.null(descrips)) {
          w <- which(descrips[,gocats.col]==gocat)
          if (length(w) != 1L) stop("couldn't find description for ", gocat)
          if (noisy) cat(descrips[w,"descrip"])
          results <- rbind(results, data.frame(gocat=gocat, count=sum(f), p.value=wt$p.value, go.description=descrips[w,"descrip"]))
        } else results <- rbind(results, data.frame(gocat=gocat, count=sum(f), p.value=wt$p.value))
        if (noisy) cat("\n")
      }
    }
  }
  results
}

get.significant <- function(results, p.val) {
  adjust <- p.adjust(results[,p.val], method="BH")
  if (sum(adjust < 0.05) == 0) {
    cat("no significant results\n")
    return(NULL)
  }
  temp <- results[adjust < 0.05,c(p.val, "count", "go.description")]
  format.data.frame(temp[order(temp[,p.val]),], digits=4)
}

# x <- read.table("test/RNAz/target.rnaplex.acc.out", header=TRUE)
x <- read.table("output/cornell/1/bwa/rnaz/target.rnaplex.acc.out", header=TRUE)
mw <- c()
for (i in 2:length(x))
{
  m <- run.mw(data.frame(gene=x[,1], score=x[,i]), gocats, descrips=descrips)
  adjust <- p.adjust(m$p.value, method="BH")
  if (sum(adjust < 0.05) > 0)
  {
    cat ("\n")
    cat (names(x)[i], ": significant category\n")
    cat ("# of categories is ", sum(adjust < 0.05), "\n")
    print(m[adjust < 0.05,])
  }
}

