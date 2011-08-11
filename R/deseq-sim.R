# DESeq simulations are described in
# http://genomebiology.com/content/supplementary/gb-2010-11-10-r106-s2.pdf

# About FDR
# Another way to look at the difference is that a p-value of 0.05 implies that 5% of all tests will result in false positives. An FDR adjusted p-value (or q-value) of 0.05 implies that 5% of significant tests will result in false positives. The latter is clearly a far smaller quantity.
# http://www.nonlinear.com/support/progenesis/samespots/faq/pq-values.aspx

sim1 <- function ()
{
ngenes <- 2000 # 20000
true_q0middle <- rexp( ngenes, 1/2000 ) # 1/250 )
sf <- c(0.924353, 0.924353) # , 1.081838, 1.081838) # c( .5, 1.7, 1.4, .9 )
conds <- c("A", "B") # c( "A", "A", "B", "B" )
# conds <- c( "A", "A", "B", "B" )
true_isDE <- runif( ngenes ) < 0.025 # < .3
true_lfc <- rnorm( ngenes, sd=.95 ) * true_isDE # .7
true_q0 <- data.frame(
   A = true_q0middle * 2^(  true_lfc/2 ),
   B = true_q0middle * 2^( -true_lfc/2 ) )
alpha <- .015
rnbinomMV <- function( n, mu, v ) 
  rnbinom( n, prob = mu/v, size = mu^2/(v-mu) )
countsTable <- 
   sapply( 1:2, function(j) # 4, function(j)
      sapply( true_q0[[ conds[j] ]], function(q)
         rnbinomMV( 1, sf[j]*q, sf[j]*q + sf[j]^2 * q^2 * alpha ) ) )
nonzero <- rowSums(countsTable)>0 
}
###############################################################
# Analysis
# library(DESeq)
# library(edgeR)
# A standard analysis with DESeq:
constantAlpha1 <- function()
{
cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds , method="blind" )
res <- nbinomTest( cds, "A", "B" )
# The same with edgeR. (edgeR gets the advantage of knowing the true
# library sizes)

#dgl <- DGEList( counts=countsTable, group=conds, lib.size=sf*1e7 )
#dgl <- estimateCommonDisp(dgl)
#edgerRes <- exactTest( dgl )
#edgerResPadj <- p.adjust( edgerRes$table$p.value, method="BH" )
}

checkTypeIerror <- function()
{
deseqTypeI <- length(which( res$pval[nonzero]<.05 & !true_isDE[nonzero] )) / length(nonzero)
cat ("DESeq Type I error is", deseqTypeI, "\n")
deseqPower <- length(which( res$padj[nonzero]<.1 & true_isDE[nonzero] )) / 
   length(which( true_isDE[nonzero] ))
cat ("DESeq power is", deseqPower, "\n")
deseqFDR <- length(which( res$padj[nonzero]<.1 & !true_isDE[nonzero] )) / 
   length(which( res$padj[nonzero]<.1 ))
cat ("DESeq FDR is", deseqFDR, "\n")
}

edgeRcheckTypeIerror <- function()
{
edgeRTypeI <- length(which( edgerRes$table$p.value[nonzero]<.05 & !true_isDE[nonzero] )) / length(nonzero)
cat ("edgeR Type I error is", edgeRTypeI, "\n")
edgeRPower <- length(which( edgerResPadj[nonzero]<.1 & true_isDE[nonzero] )) / 
   length(which( true_isDE[nonzero] ))
cat ("edgeR power is", edgeRPower, "\n")
edgeRFDR <- length(which( edgerResPadj[nonzero]<.1 & !true_isDE[nonzero] )) / 
   length(which( edgerResPadj<.1 ))
cat ("edgeR FDR is", edgeRFDR, "\n")
}

sim2 <- function ()
{
ngenes <- 20000
true_q0middle <- rexp( ngenes, 1/250 )
sf <- c( .5, 1.7, 1.4, .9 )
conds <- c( "A", "A", "B", "B" )
true_isDE <- runif( ngenes ) < .3
true_lfc <- rnorm( ngenes, sd=.7 ) * true_isDE
true_q0 <- data.frame(
   A = true_q0middle * 2^(  true_lfc/2 ),
   B = true_q0middle * 2^( -true_lfc/2 ) )
alpha <- function(mu) .01 + 9/(mu+100)

countsTable <- 
   sapply( 1:4, function(j)
      sapply( true_q0[[ conds[j] ]], function(q)
         rnbinomMV( 1, sf[j]*q, sf[j]*q + sf[j]^2 * q^2 * alpha( q ) ) ) )
}


