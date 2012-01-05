\name{smutans.plotDiffExp}
\Rdversion{1.1}
\alias{smutans.plotDiffExp}

\title{
   Function for plotting estimated pooled dispersion values with regression
   curve. 
}
\description{
   This function plots dispersion values with respect to normalized mean counts
   over all of the sample.
}
\usage{
smutans.plotDiffExp(x, file="1.ps")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
      a Smutans with an experimental design
}
  \item{file}{
      a postscript file name
}
}
\details{
   See \code{\link{smutans.plotDiffExpForMatrices}} for more technical informations
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
%smutans.plotDiffExp(x, file="1.ps")
}
