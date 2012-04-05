save.image(file = "tomato.RData")
load("tomato.RData")

library(rtracklayer)
library(GenomicRanges)
library(VariantAnnotation)
library("GenomicFeatures")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
gene.range3 <- import.gff3("/v4scratch/sc2265/rnaseq/output/tomato/1/data/ITAG2.3_gene_models.gff3.txt.no-negative")
#load("gene.range3.Rd")
vcffile <- "FASTQ001.var.flt.vcf.gz"
tomatoFaFile <- FaFile( "/v4scratch/sc2265/rnaseq/output/tomato/1/data/S_lycopersicum_chromosomes.2.40.fna" )
# tomatoFaFile <- FaFile( "S_lycopersicum_chromosomes.2.40.fna" )

# makeTranscriptDb
gr.mRNA <- gene.range3[gene.range3$type=="mRNA",]
gr.exon <- gene.range3[gene.range3$type=="exon",]
gr.cds <- gene.range3[gene.range3$type=="CDS",]

# Find exons with proper CDS.
grExon <- GRanges( seqnames = gr.exon$space, 
               ranges = IRanges(start(gr.exon), end = end(gr.exon), names = gr.exon$ID), 
               strand = strand(gr.exon), 
               ID = gr.exon$ID, Parent = unlist(gr.exon$Parent) )
grCDS <- GRanges( seqnames = gr.cds$space, 
               ranges = IRanges(start(gr.cds), end = end(gr.cds), names = gr.cds$ID), 
               strand = strand(gr.cds) )
ol <- findOverlaps(grCDS,grExon,type="within")
# Unique exons
# Choose unique exons and CDS would be unique.
head(unique(as.matrix(ol)[,2])) 
x <- c()
print(length(unique(as.matrix(ol)[,2])))
j <- 0
for (i in unique(as.matrix(ol)[,2])) {
  j <- j + 1
  x <- c(x, which(as.matrix(ol)[,2]==i)[1])
  if (j %% 1000 == 0) {
    print(j)
  }
}
ol2 <- as.matrix(ol)[x,]
rm (j,i,x)
grExon.unique <- grExon[ol2[,2],]
grCDS.unique <- grCDS[ol2[,1],]

# tx_id of splicings
x <- elementMetadata(grExon.unique)[,"Parent"]
y <- paste("mRNA",gr.mRNA$Name,sep=":")
z1 <- rep(0,length(elementMetadata(grExon.unique)[,"ID"]))
z2 <- rep(0,length(elementMetadata(grExon.unique)[,"ID"]))
j <- 0
for (i in y) {
  j <- j + 1
  z1[which(x %in% i)] <- which(y %in% i)
  z2[which(x %in% i)] <- seq(length(which(x %in% i)))
  if (j %% 1000 == 0) {
    print(j)
  }
}
tx_id <- as.integer(z1)
exon_rank <- z2
rm(x,y,z1,z2,i,j)

splicings <- 
  data.frame( tx_id=tx_id,
              exon_rank=exon_rank,
              exon_start=start(grExon.unique),
              exon_end=end(grExon.unique),
              cds_start=start(grCDS.unique),
              cds_end=end(grCDS.unique))

chrominfo <-
  data.frame( chrom=names(gr.mRNA),
              length=c(21805821L, 90304244L, 49918294L, 64840714L, 64064312L, 65021438L, 46041636L, 65268621L, 63032657L, 67662091L, 64834305L, 53386025L, 65486253L), 
              is_circular=rep(FALSE,length(names(gr.mRNA))) )


# Chromosome names
# names(gr.mRNA)
# Number of Chromosomes
# length(gr.mRNA)

# Create mRNA transcripts
transcripts <- 
  data.frame( tx_id=seq(length(gr.mRNA$Name)),
              tx_name=gr.mRNA$Name,
              tx_chrom=gr.mRNA$space, 
              tx_strand=gr.mRNA$strand,
              tx_start=start(gr.mRNA),
              tx_end=end(gr.mRNA) )

# Crate TranscriptDb object
txdb1 <- makeTranscriptDb (transcripts, splicings, chrominfo=chrominfo)

# Call Variants
vcf <- readVcf(vcffile, "tomato")
rd <- rowData(vcf)
loc <- locateVariants(vcf, txdb1, AllVariants())
aacoding <- predictCoding(vcf, txdb1, seqSource=tomatoFaFile)

write.csv(as.data.frame(loc),file="loc.csv")
elementMetadata(aacoding)[,"proteinLoc"] <- sapply(elementMetadata(aacoding)[,"proteinLoc"],paste,collapse=",") 

names(aacoding)[names(aacoding) =="SL2.40ch08:54030169"] <- c("SL2.40ch08:54030169-1","SL2.40ch08:54030169-2")
names(aacoding)[names(aacoding) =="SL2.40ch08:54030177"] <- c("SL2.40ch08:54030177-1","SL2.40ch08:54030177-2")

write.csv(as.data.frame(aacoding),file="aacoding.csv")
