library(goseq)

# Measured genes as a named vector (assayed.genes)
# DE genes as a named vector (de.genes)
# Gene lengths for organism, gene identifier
# Association between categories and genes
# Companion packages: GenomicFeatures, Rsamtools
#   for the summarization of mapped reads into a table of counts, or reads per
#   gene
# Length data: the length of each gene (length.genes)
# Category mappings: the mapping between the categories and genes
#   Two-column data frame - gene IDs, and associated gene categories
#   (go.genes)
feature.genes = read.table(file="feature-genome.out-geneonly",head=F)
length.genes = feature.genes$V3 - feature.genes$V2
assayed.genes = feature.genes$V4

# assayed.genes = scan(file="all.gene.2",what="character")
de.genes = scan(file="de.gene.2",what="character",quiet=T)
gene.vector = as.integer(assayed.genes %in% de.genes)
names(gene.vector) = assayed.genes

# supportedGenomes()
# supportedGeneIDs() 
# go_map = getgo(names(genes),"hg18","ensGene")
# go = goseq(pwf,gene2cat=go_map)

go.genes = read.table(file="smutans.gene2go",head=F)
cat.desc = read.table(file="smutans.go2ngene",head=F,sep="\t",quote="")
rm(feature.genes)
pwf = nullp(gene.vector, bias.data=length.genes)


go.hypergeometric = goseq(pwf,gene2cat=go.genes,method="Hypergeometric") # No length bias correction
go = goseq(pwf,gene2cat=go.genes,method="Wallenius") # Length bias correction - Approximation
go.sample = goseq(pwf,gene2cat=go.genes,method="Sampling",repcnt=10000) # Length bias correction - Sampling
go.fdr = go[p.adjust(go$over_represented_pvalue,method="BH")<.05,]

for (i in go.fdr$category)
{
  cat(as.character(cat.desc$V3[cat.desc$V1==i]),"\n")
}

# Plot of p-values for sampling and Wallenius
# plot(log10(go[,2]), log10(go.sample[match(go.sample[,1],go[,1]),2]), xlab="log10(Wallenius p-values)", ylab="log10(Sampling p-values)", xlim=c(-3,0))
# abline(0,1,col=3,lty=2)
# plot(log10(go[,2]), log10(go.hypergeometric[match(go.hypergeometric[,1],go[,1]),2]), xlab="log10(Wallenius p-values)", ylab="log10(Sampling p-values)", xlim=c(-3,0))

