library(DEGseq)
# TW1-Glucose vs. TW1-Galactose or 11 vs. 03
tw1glu2gal <- function ()
{
  tw1galactose <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/FASTQ03.bed.degseq"
  tw1glucose <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/FASTQ11.bed.degseq"
  refFlat <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/run-analysis/refFlat.txt"
  mapResultBatch1 <- c(tw1glucose)
  mapResultBatch2 <- c(tw1galactose)
  DEGseq(mapResultBatch1, mapResultBatch2, fileFormat="bed", refFlat=refFlat, outputDir="tw1glucosegalactose", method="LRT")
}

# UA159 vs. TW1-Glucose or 01 vs. 11
wt2tw <- function ()
{
  wt <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/FASTQ01.bed.degseq"
  tw <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/FASTQ11.bed.degseq"
  refFlat <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/run-analysis/refFlat.txt"
  mapResultBatch1 <- c(wt)
  mapResultBatch2 <- c(tw)
  DEGseq(mapResultBatch1, mapResultBatch2, fileFormat="bed", refFlat=refFlat, outputDir="wt2tw", method="LRT")
}

# OMZ175 vs. OMZ/HKRR or 09 vs. 07
omz <- function ()
{
  omz175 <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/FASTQ09.bed.degseq"
  omzhkrr <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/data/FASTQ07.bed.degseq"
  refFlat <- "/Users/goshng/Documents/Projects/rnaseq/output/smutans12/run-analysis/refFlat.txt"
  mapResultBatch1 <- c(omz175)
  mapResultBatch2 <- c(omzhkrr)
  DEGseq(mapResultBatch1, mapResultBatch2, fileFormat="bed", refFlat=refFlat, outputDir="omz", method="LRT")
}

# I do not have UA159-Galactose strain to compare UA159-Glucose.
# tw1glu2gal()
# wt2tw()
omz()
