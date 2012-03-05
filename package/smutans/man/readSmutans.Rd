\name{readSmutans}
\Rdversion{1.1}
\alias{readSmutans}

\title{
  Read a count table
}
\description{
  This function creates CountDataSet object by reading a count table. See DESeq
  for the class CountDataSet.
}
\usage{
readSmutans (countsFile, 
             indexFile, 
             condition, 
             firstFactor,
             firstFactorLabel,
             secondFactor,
             secondFactorLabel,
             name="", lab="", contact="", title="", url="",
             abstract="", pubmedid="",
             cutoff.count=10)
}
\arguments{
  \item{countsFile}{
    a string for a count table file.
}
  \item{indexFile}{
    a string for a count table index file. The file contains a single text line
    of as many factors as columns in the count table. 
}
  \item{condition}{
    a vector of strings for a subset of the indices in the count table index.
    The vector can be of size 2, 3, or 4. It is used to subset columns of the
    count table file. 
}
  \item{firstFactor}{
    a vector of two strings for grouping the same first condition.
}
  \item{firstFactorLabel}{
    a vector of two strings for denoting the two groups of the first condition.
}
  \item{secondFactor}{
    a vector of two strings for grouping the same second condition.
}
  \item{secondFactorLabel}{
    a vector of two strings for denoting the two groups of the second condition.
}
  \item{name}{
    [optional] a string for the name of the comparison experiment.
}
  \item{lab}{
    [optional] a string for the names of laboratories.
}
  \item{contact}{
    [optional] a string for the names of contact persons.
}
  \item{title}{
    [optional] a string for the title of a journal article.
}
  \item{url}{
    [optional] a string for the address of the source of the data.
}
  \item{abstract}{
    [optional] a text of the abstract.
}
  \item{pubmedid}{
    [optional] a string for the PubMed ID.
}
  \item{cutoff.count}{
    [optional] a minimum value for the number of short reads for a gene that is used in the
    comparison. Genes with short reads less than this number are filtered out
    before detemining differentially expressed genes.
}
}
\details{
  A CountDataSet object must be prepared with two factors using newCountDataSet
  function provided by R package DESeq. The title of the analysis would be
  useful to distinguish different analyses later. 

  If condition is a vector of size 4, then a two-factor count table is created:
  e.g., condition = c("UA159GLU", "UA159GAL", "TW1GLU", "TW1GAL"). One factor is
  whether the samples are from ua159 type or tw1 type. Group first
  by this factor by specifying the two levels of the first factor:
  firstFactorLabel = c("ua159", "tw1"). Choose two conditions that belong to
  ua159 or the first element of firstFactorLabel:
  firstFactor = c("UA159GLU", "UA159GAL"). Secondly, group the four conditions
  whether they beong to glucose or galactose:
  secondFactorLabel = c("glucose","galactose"). Choose two conditions that
  belong to glucose or the first element of secondFactorLabel:
  secondFactor = c("UA159GLU", "TW1GLU"). 
 
  If condition is a vector of size 3, then a one-factor count table is craeted.
  For example, we have condition = c("835NP", "835P", "UA159noCSP"). First group
  all of the three conditions in one type ua159:
  firstFactor = c("835NP", "835P", "UA159noCSP"), and
  firstFactorLabel = c("ua159").
  Because we do not have any other factor to group in the second factor:
  secondFactor = c(), and secondFactorLabel = c("835np","835p", "ua159nocsp").

  If condition is a vector of size 2, then a one-factor count table is craeted.
}
\value{
An CountDataSet object.
}
\author{
   Sang Chul Choi, schoi@cornell.edu
}
\examples{
# To read in a two-factor data
inDir <- system.file("extdata", package="smutans", mustWork=TRUE)
smutansGenes <- 
  readSmutans (countsFile=file.path(inDir, "count-ua159-bwa.txt"),
               indexFile =file.path(inDir, "count-ua159.txt.index"),
               condition = c("UA159GLU", "UA159GAL", "TW1GLU", "TW1GAL"),
               firstFactor = c("UA159GLU", "UA159GAL"),
               firstFactorLabel = c("ua159", "tw1"),
               secondFactor = c("UA159GLU", "TW1GLU"),
               secondFactorLabel = c("glucose","galactose"))
}
