library(edgeR)
library(limma)
targets = read.delim(file = "Targets.txt",sep=" ",stringsAsFactors = FALSE)
d = readDGE(targets,head=F)
d <- d[rowSums(d$counts) >= 5,] # Filter out those with counts fewer than 5.
# dim(d) # Show the number of genes and the number of samples
# d$samples$lib.size # Show the total number of reads for each of the samples
# colSums(d$counts)
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d) # TMM normalization
# d$samples # Show sample
d <- estimateCommonDisp(d)
# names(d) # Show names of d
# d$samples$lib.size # Show library sizes
# d$common.lib.size # Show common library size
# sqrt(d$common.dispersion)
de.com <- exactTest(d)
# names(de.com) # Show names
# names(de.com$table) # Show names
