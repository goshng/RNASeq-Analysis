\name{smutans.de2type}
\Rdversion{1.1}
\alias{smutans.de2type}

\title{
   Test for differences between the base means for two conditions using
   a two-factor data set
}
\description{
   This function tests for differences between the base means of two
   conditions using two-factor data sets. One factor is named type, and another
   condition. We consider only one of the two types to compare samples based on
   experimental condition. 
}
\usage{
smutans.de2type(x, type="ua159", condA="glucose", condB="galactose")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
      a Smutans with an experimental design
}
  \item{type}{
      one of the types in 'x'
}
  \item{condA}{
      one of the conditions in 'x'
}
  \item{condB}{
      another one of the conditions in 'x'
}
}
\details{
   See \code{\link{smutans.de2typeForMatrices}} for more technical informations
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
%x <- smutans.de2type(x, type="ua159", condA="glucose", condB="galactose")
}
