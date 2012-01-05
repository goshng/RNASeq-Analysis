# library(edgeR)
# library(limma)



ua159.tw1.glu.gal.test = function()
{
################################################################################
# A Test of multifactor analysis: UA159
# 
# 0 for glucose, and 1 for galactose
# w for UA159, and m for TW1
targets = read.delim(file = "ua159-tw1-glu-gal.txt",sep=" ",stringsAsFactors = FALSE)
sugar <- factor(c(0,0,1,0,0,1,1,1,0,0,1,1,1))
# mutant <- factor(c(w,m,m,w,w,w,w,w,m,m,m,m,w))
d = readDGE(targets,head=F)
d <- calcNormFactors(d)
cpm.sugar <- cpm(d)
print(rownames(d[rowSums(cpm.sugar > 1) < 2,]$counts))
d <- d[rowSums(cpm.sugar > 1) >= 2,]
design <- model.matrix(~sugar + d$samples$group)
rownames(design) <- rownames(d$samples)
colnames(design)[2] <- "sugar"
colnames(design)[3] <- "mutant"
d <- estimateGLMCommonDisp(d,design)
glmfit.sugar <- glmFit(d,design,dispersion=d$common.dispersion)
lrt.sugar <- glmLRT(d,glmfit.sugar,coef=3)
lrt.mutant <- glmLRT(d,glmfit.sugar,coef=2)
# topTags(lrt.sugar)
}

################################################################################
# A Test
np835.proc <- function ()
{
targets = read.delim(file = "np835.txt",sep=" ",stringsAsFactors = FALSE)
d.835 = readDGE(targets,head=F)
d.835 <- calcNormFactors(d.835)
cpm.835 <- cpm(d.835)
print(rownames(d.835[rowSums(cpm.835 > 1) < 2,]$counts))
d.835 <- d.835[rowSums(cpm.835 > 1) >= 2,]
plotMDS(d.835)

# d <- d[rowSums(d$counts) >= 5,] # Filter out those with counts fewer than 5.
# dim(d) # Show the number of genes and the number of samples
# d$samples$lib.size # Show the total number of reads for each of the samples
# colSums(d$counts)
# d$samples$lib.size <- colSums(d$counts)
# d <- calcNormFactors(d) # TMM normalization
# d$samples # Show sample
d.835 <- estimateCommonDisp(d.835)
# names(d) # Show names of d
# d$samples$lib.size # Show library sizes
# d$common.lib.size # Show common library size
# sqrt(d$common.dispersion)
de.835 <- exactTest(d.835)
print(topTags(de.835))
# names(de.com) # Show names
# names(de.com$table) # Show names
}

np835.proc ()
