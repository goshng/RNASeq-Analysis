\name{smutans.de2Genes}
\Rdversion{1.1}
\alias{smutans.de2Genes}

\title{
  Returns a list of differentially expressed genes
}
\description{
  A 
  This function tests for differences between the base means of two
  conditions using two-factor data sets. One factor is named type, and another
  condition. The function behaves differently depending arguments. We could
  compare two conditions given one-type. We could also compare two types given
  one-condition. We could compare two conditions considering both types. As you
  expect, we can compare two types consdering both conditions.
}
\usage{
smutans.de2Genes(object, qval=0.05)
}
\arguments{
  \item{object}{
    a Smutans with an experimental design of two factors: one is type, and
    another is condition. See function \code{smutans.prepareDataUA159TW1} for
    example.
}
  \item{qval}{
    a cutoff adjusted p-value for differentially expressed genes.
}
}
\details{
  We use this to find a list of differentially expressed genes using a cutoff
  adjusted p-value. This list could be used with R package goseq.
}
\value{
A list of gene names.
}
\author{
   Sang Chul Choi, schoi@cornell.edu
}
\examples{
%data( "smutansGenes" )
%x <- newSmutans(smutansGenes)
%x <- smutans.de2Genes(x, type="ua159", condA="glucose", condB="galactose")
}
