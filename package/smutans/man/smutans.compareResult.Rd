\name{smutans.compareResult}
\Rdversion{1.1}
\alias{smutans.compareResult}

\title{
   Compares two analyses in their power to detect differentiall expressed genes
}
\description{
   This function compare two analyses. 
}
\usage{
smutans.compareResult(x,y,qval=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
      a Smutans with an experimental design
}
  \item{y}{
      another Smutans with an experimental design
}
  \item{qval}{
      a numeric value for q-value cutoff
}
}
\details{
   See \code{\link{smutans.compareResult}} for more technical informations
   A resulting table is a 2-by-2 matrix. If a variable is set to the output of
   the function, say tab, then we can say the following: the analyses find
   tab[2,2] genes in common. Of the length(tab) genes, tab[1,2] genes are found
   in analysis x or the first argument of the function. tab[2,1] genes are found
   in analysis y or the second argument of the function. tab[1,1] genes are not
   found to be differentially expressed in either analysis.
}
\value{
A 2-by-2 contigency table.
}
\author{
   Sang Chul Choi, schoi@cornell.edu
}
\examples{
%data( "smutansGenes" )
%x <- newSmutans(smutansGenes)
%
%tab <- smutans.compareResult(x, y)
%addmargins(tab)
}
