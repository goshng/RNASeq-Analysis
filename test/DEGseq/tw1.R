tw1galactose <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/FASTQ03.bed.degseq"
tw1glucose <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/FASTQ11.bed.degseq"
refFlat <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/run-analysis/refFlat.txt"
mapResultBatch1 <- c(tw1glucose)
mapResultBatch2 <- c(tw1galactose)
DEGseq(mapResultBatch1, mapResultBatch2, fileFormat="bed", refFlat=refFlat, outputDir="tw1glucosegalactose", method="LRT")

