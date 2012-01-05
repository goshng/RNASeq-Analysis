\name{smutans.de2}
\Rdversion{1.1}
\alias{smutans.de2}

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
x <- smutans.de2(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
      a Smutans with an experimental design
}
}
\details{
   See \code{\link{smutans.de2}} for more technical informations
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
