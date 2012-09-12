\name{smutans.keg.genes}
\docType{data}
\alias{smutans.keg.genes}
\title{Association table of gene ontology terms and genes}
\description{
\code{smutans.keg.genes} is a table of three columns of gene name, functional
category, and a p-value: e.g., SMU.01 GO:0000166 0.

Use smutans.keg.genes <- smutans.keg.genes[c(4,3)] for KEGG functional category
analysis.
}
\usage{
data("smutans.keg.genes")
}

\format{ExonCountSet, CountDataSet}

\source{Processed data from NCBI Gene Expression Omnibus under accession numbers GSM461176 to GSM461181.}
\references{Brooks et al. Conservation of an RNA regulatory map between Drosophila and mammals. Genome Research, 2010}
\keyword{datasets}
