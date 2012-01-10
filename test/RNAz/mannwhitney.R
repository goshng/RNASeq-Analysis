gocats <- read.table("/usr/projects/strep/SpyMGAS315/SpyMGAS315_go_bacteria.txt", header=FALSE, col.names=c("gene", "go", "pval"))
gocats <- gocats[gocats$pval <= 1.0e-5,]
descrips <- read.table("/usr/projects/strep/SpyMGAS315/SpyMGAS315_go_category_names.txt", sep="\t", header=FALSE, col.names=c("go", "count", "descrip"), quote="")
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



  
x <- read.table("in.gene", header=FALSE)
num.unique.mw <- run.mw(data.frame(gene=x[,1], score=x[,12]), gocats, descips=descrips)
x <- read.table("ri1-refgenome4-map.gene", header=FALSE)
sde.spy.mw <- run.mw(data.frame(gene=x[,1], score=rowSums(x[,10:13])), gocats, descrips=descrips)
spy.sde.mw <- run.mw(data.frame(gene=x[,1], score=rowSums(x[,14:17])), gocats, descrips=descrips)
all.mw <- run.mw(data.frame(gene=x[,1], score=x[,18]), gocats)

# these should all be sorted in the same way but make sure
all.mw <- all.mw[order(all.mw$gocat),]
num.unique.mw <- num.unique.mw[order(num.unique.mw$gocat),]
sde.spy.mw <- sde.spy.mw[order(sde.spy.mw$gocat),]
spy.sde.mw <- spy.sde.mw[order(spy.sde.mw$gocat),]
numrow <- nrow(all.mw)
if (numrow != nrow(spy.sde.mw) || numrow != nrow(sde.spy.mw) || numrow != nrow(num.unique.mw) ||
    sum(all.mw$gocat==spy.sde.mw$gocat) != numrow ||
    sum(all.mw$gocat==sde.spy.mw$gocat) != numrow ||
    sum(all.mw$gocat==num.unique.mw$gocat) != numrow)
  stop("results don't have all the same elements in the same order")

all.results <- data.frame(gocat=all.mw$gocat, count=all.mw$count,
                          p.num.unique=num.unique.mw$p.value,
                          p.all=all.mw$p.value,
                          p.sde.spy=sde.spy.mw$p.value,
                          p.spy.sde=spy.sde.mw$p.value,
                          go.description=all.mw$go.description)

outfile <- "significant.txt"
cat("p-value", "count", "description", file=outfile, sep="\t")
cat("\n", file=outfile, append=TRUE)
for (stat in c("p.num.unique", "p.all", "p.sde.spy", "p.spy.sde")) {
  cat("\n", file=outfile, append=TRUE)
  cat("#", stat, "\n", file="significant.txt", append=TRUE)
  write.table(get.significant(all.results, stat), "significant.txt", row.names=FALSE, quote=FALSE, sep="\t",
              append=TRUE, col.names=FALSE)
}
write.table(format.data.frame(all.results, digits=4), file="mannwhitney_results.txt", quote=FALSE, row.names=FALSE, sep="\t")


# Now check virulence genes
x <- read.table("in.gene", header=FALSE)
gocats <- read.table("/usr/projects/strep/virulent_genes.txt", header=TRUE)
gocats$gene <- as.character(gocats$gene)
gocats <- gocats[is.element(gocats$gene, as.character(x[,1])),]
x <- x[is.element(as.character(x[,1]), gocats$gene),]
run.mw(data.frame(gene=x[,1], score=x[,12]), gocats, gocats.col="virulent")

x <- read.table("ri1-refgenome4-map.gene", header=FALSE)
gocats <- read.table("/usr/projects/strep/virulent_genes.txt", header=TRUE)
gocats$gene <- as.character(gocats$gene)
gocats <- gocats[is.element(gocats$gene, as.character(x[,1])),]
x <- x[is.element(as.character(x[,1]), gocats$gene),]
run.mw(data.frame(gene=x[,1], score=rowSums(x[,10:13])), gocats, gocats.col="virulent")
run.mw(data.frame(gene=x[,1], score=rowSums(x[,14:17])), gocats, gocats.col="virulent")
run.mw(data.frame(gene=x[,1], score=x[,18]), gocats, gocats.col="virulent")
