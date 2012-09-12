smutans.compareResult <- function ( x, y, qval=0.1 )
{
  padj.resA <- smutans.padj( x )
  padj.resB <- smutans.padj( y )
  tab <- table( padj.resA < qval, 
                padj.resB < qval,
                dnn=c(smutans.title(x), smutans.title(y)) )
  tab
}

smutans.prepareData34 <- function ()
{
  # Count data preparation
  countsFile <- "smutans/inst/extdata/count.txt"
  conds <- scan(file=paste(countsFile,"index",sep="."), what="character")

  # subconds <- conds %in% c("UA159GLU", "UA159GAL", "TW1GLU", "TW1GAL")
  countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
  rownames(countsTable) <- countsTable$gene
  countsTable <- countsTable[,-1]
  # countsTable <- countsTable[,subconds]

  # conds <- conds[subconds]
  # conds.type <- conds %in% c("UA159GLU", "UA159GAL")
  factor.type <- rep("ua159",length(conds))
  # factor.type[!conds.type] <- "tw1"
  #conds.condition <- conds %in% c("UA159GLU", "TW1GLU")
  #factor.condition <- rep("glucose",length(conds))
  factor.condition <- conds

  samples <- data.frame(type=factor.type,condition=factor.condition)
  rownames(samples) <- colnames(countsTable) 

  design <- samples
  sm34Genes <- newCountDataSet( countsTable, design )
  expdata = new("MIAME", 
     name="S. mutans UA159, TW1, Glucose, and Galactose", 
     lab="University of Florida, and Cornell University", 
     contact="Drs. Robert Burne, Michael Stanhope, and Adam Siepel", 
     title="Streptococcus mutans RNA-Seq Studies", 
     url="http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=XXX", 
     abstract="RNA-seq of 34 biological replicates from Streptococcus mutans")
  pubMedIds(expdata) <- "999999999"
  experimentData(sm34Genes) <- expdata
  save(sm34Genes, file=file.path("smutans", "data", "sm34Genes.RData"))
}

smutans.prepareDataUA159TW1 <- function (cutoff.count=10)
{
  # Count data preparation
  countsFile <- "smutans/inst/extdata/count.txt"
  countsFile2 <- "smutans/inst/extdata/count2.txt"
  conds <- scan(file=paste(countsFile,"index",sep="."), what="character")
  subconds <- conds %in% c("UA159GLU", "UA159GAL", "TW1GLU", "TW1GAL")

  countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
  rownames(countsTable) <- countsTable$gene
  countsTable <- countsTable[,-1]
  countsTable <- countsTable[,subconds]

  countsTable2 <- read.delim (countsFile2, header=TRUE, stringsAsFactors=TRUE)
  rownames(countsTable2) <- countsTable2$gene
  countsTable2 <- countsTable2[,-1]
  countsTable2 <- countsTable2[,subconds]
  countsTable2 <- countsTable2[rowSums(countsTable2) >=cutoff.count,]

  conds <- conds[subconds]
  conds.type <- conds %in% c("UA159GLU", "UA159GAL")
  factor.type <- rep("ua159",length(conds))
  factor.type[!conds.type] <- "tw1"
  conds.condition <- conds %in% c("UA159GLU", "TW1GLU")
  factor.condition <- rep("glucose",length(conds))
  factor.condition[!conds.condition] <- "galactose"
  samples <- data.frame(type=factor.type,condition=factor.condition)
  rownames(samples) <- colnames(countsTable) 

  design <- samples
  smutansGenes <- newCountDataSet( countsTable, design )
  expdata = new("MIAME", 
     name="S. mutans UA159, TW1, Glucose, and Galactose", 
     lab="University of Florida, and Cornell University", 
     contact="Drs. Robert Burne, Michael Stanhope, and Adam Siepel", 
     title="Streptococcus mutans RNA-Seq Studies", 
     url="http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=XXX", 
     abstract="RNA-seq of 13 biological replicates from Streptococcus mutans")
  pubMedIds(expdata) <- "999999999"
  experimentData(smutansGenes) <- expdata
  save(smutansGenes, file=file.path("smutans", "data", "smutansGenes.RData"))

  smutansGenes2 <- newCountDataSet( countsTable2, design )
  expdata = new("MIAME", 
     name="S. mutans UA159, TW1, Glucose, and Galactose based on transcripts and small RNAs", 
     lab="University of Florida, and Cornell University", 
     contact="Drs. Robert Burne, Michael Stanhope, and Adam Siepel", 
     title="Streptococcus mutans RNA-Seq Studies", 
     url="http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=XXX", 
     abstract="RNA-seq of 13 biological replicates from Streptococcus mutans")
  pubMedIds(expdata) <- "88888888"
  experimentData(smutansGenes2) <- expdata
  save(smutansGenes2, file=file.path("smutans", "data", "smutansGenes2.RData"))
}

smutans.makeData <- function( ngenes=500 ) 
{
  q0 <- rexp( ngenes, rate=1/250 )
  is_DE <- runif( ngenes ) < .3
  lfc <- rnorm( ngenes, sd=2 )
  q0A <- ifelse( is_DE, q0 * 2^(  lfc/2 ), q0 )
  q0B <- ifelse( is_DE, q0 * 2^( -lfc/2 ), q0 )
  true_sf <- c( 1., 1.3, .7, .9, 1.6 )   
  conds <- c( "A", "A", "B", "B", "B" )
  m <- t( sapply( 1:ngenes, function(i) 
     sapply( 1:5, function( j )
        rnbinom( 1, mu = true_sf[j] * ifelse( conds[j]=="A", q0A[i], q0B[i] ), 
           size = 1/.2 ) ) ) )
  colnames(m) <- c( "A1", "A2", "B1", "B2", "B3" )
  rownames(m) <- paste( "gene", 1:ngenes, 
    ifelse( is_DE, "T", "F" ), sep="_" )
  newCountDataSet( m, conds )
}


smutans.prepareGoseq <- function ()
{
  f <- "smutans/inst/extdata/feature-genome.out-geneonly"
  smutans.feature.genes = read.table(file=f,head=F)
  save(smutans.feature.genes, file=file.path("smutans", "data", "smutans.feature.genes.RData"))
  f <- "smutans/inst/extdata/smutans.gene2go"
  smutans.go.genes <- read.table(file=f,head=F)
  save(smutans.go.genes, file=file.path("smutans", "data", "smutans.go.genes.RData"))
  f <- "smutans/inst/extdata/smutans.go2ngene"
  smutans.cat.desc = read.table(file=f,head=F,sep="\t",quote="")
  save(smutans.cat.desc, file=file.path("smutans", "data", "smutans.cat.desc.RData"))
  f <- "smutans/inst/extdata/target.rnaplex.acc.out"
  smutans.genes.criteria = read.table(file=f,head=T,sep="\t",quote="")
  save(smutans.genes.criteria, file=file.path("smutans", "data", "smutans.genes.criteria.RData"))

  f <- "smutans/inst/extdata/smu00001.keg.gene2cat"
  smutans.keg.genes <- read.delim(file=f,head=F)
  save(smutans.keg.genes, file=file.path("smutans", "data", "smutans.keg.genes.RData"))
  f <- "smutans/inst/extdata/smu00001.keg.catdesc"
  smutans.keg.desc <- read.delim(file=f,head=F)
  save(smutans.keg.desc, file=file.path("smutans", "data", "smutans.keg.desc.RData"))
  # smutans.keg.genes <- x[c(4,3)]
}

smutans.run.mw <- function(x, gocats, gocats.col="go", descrips=NULL, noisy=FALSE) {
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

smutans.get.significant <- function(results, p.val) {
  adjust <- p.adjust(results[,p.val], method="BH")
  if (sum(adjust < 0.05) == 0) {
    cat("no significant results\n")
    return(NULL)
  }
  temp <- results[adjust < 0.05,c(p.val, "count", "go.description")]
  format.data.frame(temp[order(temp[,p.val]),], digits=4)
}

smutans.mannwhitney <- function( qval=0.05,genes.criteria,go.genes,cat.desc ) {
  gocats <- go.genes
  colnames(gocats) <- c("gene", "go", "pval")
  gocats <- gocats[gocats$pval <= 1.0e-5,]
  descrips <- cat.desc
  colnames(descrips) <- c("go", "count", "descrip")
  gocats$go <- as.character(gocats$go)
  descrips$go <- as.character(descrips$go)

# This does functional category analyses. 
#gobacteriatxt <- "test/RNAz/smutans.gene2go"
#gocattxt <- "test/RNAz/smutans.go2ngene"
#
#gocats <- read.table(gobacteriatxt, header=FALSE, col.names=c("gene", "go", "pval"))
#gocats <- gocats[gocats$pval <= 1.0e-5,]
#descrips <- read.table(gocattxt, sep="\t", header=FALSE, col.names=c("go", "count", "descrip"), quote="")
#
#gocats$go <- as.character(gocats$go)
#descrips$go <- as.character(descrips$go)
# x <- read.table("test/RNAz/target.rnaplex.acc.out", header=TRUE)

    # x <- read.table("output/cornell/1/bwa/rnaz/target.rnaplex.acc.out", header=TRUE)
  x <- genes.criteria
  mw <- c()
  for (i in 2:length(x))
  {
    m <- smutans.run.mw(data.frame(gene=x[,1], score=x[,i]), gocats, descrips=descrips)
    adjust <- p.adjust(m$p.value, method="BH")
    m <- data.frame(m, q.val=adjust)

    if (sum(adjust < 0.05) > 0)
    {
      #cat ("\n")
      #cat (names(x)[i], ": significant category\n")
      #cat ("# of categories is ", sum(adjust < qval), "\n")
      #print(m[adjust < 0.05,])

      cat (names(x)[i])
      temp <- m[adjust < 0.05,c("q.val", "count", "gocat", "go.description")]
      temp2 <- format.data.frame(temp[order(temp[,"q.val"]),], digits=2)
      for (i in 1:length(temp2$q.val))
      {
        cat (" & ")
        cat (temp2[i,]$q.val, " & ", temp2[i,]$count, " & ")
        cat (temp2[i,]$gocat, " & ", temp2[i,]$go.description, "\\\\\n")
      }
    }
  }
}

smutans.prepareTranscript <- function ()
{
  # ParseRNAseq input pileup or coverage file
  parseRNAseqInFile <- "smutans/inst/extdata/FASTQ001.parsernaseq.pileup"
  x <- read.table(parseRNAseqInFile, head=FALSE)
  pileup.v <- Rle(x$V2)

  # Predicted transcripts
  smutansData.txAll <- NULL
  for (i in c(1,2,3,19,20,21,22,23,24,25,26,27,28)) {
    parseRNAseqOutFile <- paste("smutans/inst/extdata/FASTQ", 
                                sprintf("%03d",i),
                                ".parsernaseq1",sep="")
    x <- read.table(parseRNAseqOutFile, head=FALSE)
    number.tx <- length(rownames(x))
    tx.state <- sub("_[[:digit:]]*", "", substring(x$V9,6))
    tx.state <- sub("C", "", tx.state)
    smutansData.tx <- GRanges( seqnames = Rle("chr1", number.tx),
                                ranges =
                                  IRanges(start=x$V4,
                                          end=x$V5,
                                          names=sprintf("tx%04d",seq(number.tx))),
                                strand = Rle(strand("*"),number.tx),
                                score=x$V6,
                                state=as.integer(tx.state)
                              )
    seqlengths(smutansData.tx) <- length(pileup.v)
    if (i == 1) {
      smutansData.txAll <- list(smutansData.tx)
    } else {
      smutansData.txAll[[length(smutansData.txAll)+1]] <- smutansData.tx
    }
  }
  smutansData.tx <- smutansData.txAll[[1]]

  f <- "smutans/inst/extdata/feature-genome.out-geneonly"
  smutans.feature.genes = read.table(file=f,head=F)

  number.tx <- length(rownames(smutans.feature.genes))
  smutansData.txGenes <- GRanges( seqnames = Rle("chr1", number.tx),
                                   ranges =
                                     IRanges(start=smutans.feature.genes$V2 + 1,
                                             end=smutans.feature.genes$V3,
                                             names=smutans.feature.genes$V4),
                                   strand = Rle(strand(smutans.feature.genes$V6))
                                 )
  seqlengths(smutansData.txGenes) <- length(pileup.v)

  smutansData.txPileup <- pileup.v

  # UA159 annotated entities
  f <- "smutans/inst/extdata/ua159-annotation.txt"
  smutans.feature.annotation = read.table(file=f,head=F)
  number.tx <- length(rownames(smutans.feature.annotation))
  smutansData.txAnnotation <- GRanges( seqnames = Rle("chr1", number.tx),
                                       ranges =
                                         IRanges(start=smutans.feature.annotation$V2 + 1,
                                                 end=smutans.feature.annotation$V3,
                                                 names=smutans.feature.annotation$V4),
                                       strand = Rle(strand(smutans.feature.annotation$V5))
                                     )
  seqlengths(smutansData.txAnnotation) <- length(pileup.v)

  # UA159 RNAz entities
  f <- "smutans/inst/extdata/ua159-rnaz.txt"
  t1 = read.table(file=f,head=F)
  number.tx <- length(rownames(t1))
  smutansData.txRNAz <- GRanges( seqnames = Rle("chr1", number.tx),
                                 ranges =
                                   IRanges(start=t1$V2 + 1,
                                           end=t1$V3,
                                           names=t1$V4),
                                 strand = Rle(strand(t1$V5))
                               )
  seqlengths(smutansData.txRNAz) <- length(pileup.v)

  save(smutans.feature.genes, file=file.path("smutans", "data", "smutans.feature.genes.RData"))
  save(smutansData.txAll, file=file.path("smutans", "data", "smutansData.txAll.RData"))
  save(smutansData.tx, file=file.path("smutans", "data", "smutansData.tx.RData"))
  save(smutansData.txGenes, file=file.path("smutans", "data", "smutansData.txGenes.RData"))
  save(smutansData.txAnnotation, file=file.path("smutans", "data", "smutansData.txAnnotation.RData"))
  save(smutansData.txRNAz, file=file.path("smutans", "data", "smutansData.txRNAz.RData"))
  save(smutansData.txPileup, file=file.path("smutans", "data", "smutansData.txPileup.RData"))
}

readSmutans <- function (countsFile, 
                         indexFile, 
                         condition,
                         firstFactor,
                         firstFactorLabel,
                         secondFactor,
                         secondFactorLabel,
                         name="", lab="", contact="", title="", url="",
                         abstract="", pubmedid="",
                         cutoff.count=10)
{
  # Count data preparation
  conds <- scan(file=indexFile,what="character",quiet=TRUE)
  subconds <- conds %in% condition
  if (length(condition) == 3) {
    firstFactor <- condition
  } else if (length(condition) == 2) {
    firstFactor <- condition
  }

  countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
  rownames(countsTable) <- countsTable$gene
  countsTable <- countsTable[,-1]
  countsTable <- countsTable[,subconds]
  countsTable <- countsTable[rowSums(countsTable) >=cutoff.count,]

  conds <- conds[subconds]
  conds.type <- conds %in% firstFactor
  factor.type <- rep(firstFactorLabel[1],length(conds))

  conds.condition <- c()
  factor.condition <- c()
  if (length(condition) == 4) {
    factor.type[!conds.type] <- firstFactorLabel[2]
    conds.condition <- conds %in% secondFactor
    factor.condition <- rep(secondFactorLabel[1],length(conds))
    factor.condition[!conds.condition] <- secondFactorLabel[2]
  } else if (length(condition) == 3) {
    conds.condition <- conds %in% firstFactor[1]
    factor.condition <- rep(secondFactorLabel[1],length(conds))
    conds.condition <- conds %in% firstFactor[2]
    factor.condition[conds.condition] <- secondFactorLabel[2]
    conds.condition <- conds %in% firstFactor[3]
    factor.condition[conds.condition] <- secondFactorLabel[3]
  } else if (length(condition) == 2) {
    conds.condition <- conds %in% firstFactor[1]
    factor.condition <- rep(secondFactorLabel[1],length(conds))
    conds.condition <- conds %in% firstFactor[2]
    factor.condition[conds.condition] <- secondFactorLabel[2]
  }

  samples <- data.frame(type=factor.type,condition=factor.condition)
  rownames(samples) <- colnames(countsTable) 

  design <- samples
  smutansGenes <- newCountDataSet( countsTable, design )
  expdata = new("MIAME", 
     name=name,
     lab=lab,
     contact=contact,
     title=title,
     url=url,
     abstract=abstract)
  pubMedIds(expdata) <- pubmedid
  experimentData(smutansGenes) <- expdata
  smutansGenes



}


