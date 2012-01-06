smutans.compareResult <- function (resA, resB, qval=0.1)
{
  padj.resA <- smutans.padj( resA )
  padj.resB <- smutans.padj( resB )
  tab <- table( padj.resA < qval, 
                padj.resB < qval,
                dnn=c(smutans.title(resA), smutans.title(resB)) )
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

smutans.prepareDataUA159TW1 <- function ()
{
  # Count data preparation
  countsFile <- "smutans/inst/extdata/count.txt"
  conds <- scan(file=paste(countsFile,"index",sep="."), what="character")
  subconds <- conds %in% c("UA159GLU", "UA159GAL", "TW1GLU", "TW1GAL")
  countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
  rownames(countsTable) <- countsTable$gene
  countsTable <- countsTable[,-1]
  countsTable <- countsTable[,subconds]
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
  save(smutansGenes, file=file.path(".", "smutansGenes.RData"))
}

smutans.prepareData835NPP <- function ()
{
  # Count data preparation
  countsFile <- "smutans/inst/extdata/count.txt"
  conds <- scan(file=paste(countsFile,"index",sep="."), what="character")
  subconds <- conds %in% c("835NP", "835P", "UA159noCSP")
  countsTable <- read.delim (countsFile, header=TRUE, stringsAsFactors=TRUE)
  rownames(countsTable) <- countsTable$gene
  countsTable <- countsTable[,-1]
  countsTable <- countsTable[,subconds]
  conds <- conds[subconds]
  conds.type <- conds %in% c("835NP", "835P", "UA159noCSP")
  factor.type <- rep("ua159",length(conds))
  # All of the same factor: factor.type[!conds.type] <- "tw1"
  conds.condition <- conds %in% c("835NP")
  factor.condition <- rep("835np",length(conds))
  conds.condition <- conds %in% c("835P")
  factor.condition[conds.condition] <- "835p"
  conds.condition <- conds %in% c("UA159noCSP")
  factor.condition[conds.condition] <- "ua159nocsp"
  samples <- data.frame(type=factor.type,condition=factor.condition)
  rownames(samples) <- colnames(countsTable) 

  design <- samples
  sm835Genes <- newCountDataSet( countsTable, design )
  expdata = new("MIAME", 
     name="S. mutans 835NP, 835P, UA159 no CSP", 
     lab="University of Florida, and Cornell University", 
     contact="Drs. Robert Burne, Michael Stanhope, and Adam Siepel", 
     title="Streptococcus mutans RNA-Seq Studies", 
     url="http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=XXX", 
     abstract="RNA-seq of 8 biological replicates from Streptococcus mutans")
  pubMedIds(expdata) <- "999999999"
  experimentData(sm835Genes) <- expdata
  save(sm835Genes, file=file.path("smutans", "data", "sm835Genes.RData"))
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

