\name{smutans.de2Goseq}
\Rdversion{1.1}
\alias{smutans.de2Goseq}

\title{
   Test for differences between the base means for two conditions using
   a two-factor data set
}
\description{
   This function tests for differences between the base means of two
   conditions using two-factor data sets. One factor is named type, and another
   condition. The function behaves differently depending arguments.
}
\usage{
smutans.de2Goseq(object,qval=0.05, feature.genes, go.genes, cat.desc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
      a Smutans with an experimental design
}
  \item{qval}{
      a cutoff value for determining differentially expressed genes.
}
  \item{feature.genes}{
      a table of six columns of chrmosome, start, end, gene name, 0, and strand.
}
  \item{go.genes}{
      a table of three columns of gene name, functional category, and a p-value:
      e.g., SMU.01 GO:0000166  0
}
  \item{cat.desc}{
      a table of three columns of functional category, number of genes, and
      category description: e.g., GO:0007205 1 alcohol metabolic process
}
}
\details{
   See \code{\link{smutans.de2Goseq}} for more technical informations
}
\value{
An Smutans object.
}
\author{
   Sang Chul Choi, schoi@cornell.edu
}
\examples{
%data( "smutansGenes" )
%x <- newSmutans(smutansGenes)
%x <- smutans.de2Goseq(x, type="ua159", condA="glucose", condB="galactose")
%smutans.de2Goseq ( ua159only, qval=1e-2,
%                   feature.gene=smutans.feature.genes,
%                   go.genes = smutans.go.genes,
%                   cat.desc = smutans.cat.desc )
}
