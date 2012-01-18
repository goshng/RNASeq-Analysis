\name{smutans.mannwhitney}
\Rdversion{1.1}
\alias{smutans.mannwhitney}

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
smutans.mannwhitney(qval=0.05,genes.criteria,go.genes,cat.desc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qval}{
      a cutoff value for determining differentially expressed genes.
}
  \item{genes.criteria}{
      a table of values for a gene and one of criteria; columns are criteria.
      See e.g., smutans/inst/extdata/target.rnaplex.acc.out
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
   See \code{\link{smutans.mannwhitney}} for more technical informations
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
%x <- smutans.mannwhitney(x, type="ua159", condA="glucose", condB="galactose")
%smutans.mannwhitney ( ua159only, qval=1e-2,
%                   feature.gene=smutans.feature.genes,
%                   go.genes = smutans.go.genes,
%                   cat.desc = smutans.cat.desc )
}
