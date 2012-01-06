\name{smutans.de2TypeHeatmap}
\Rdversion{1.1}
\alias{smutans.de2TypeHeatmap}

\title{
   Draw a heatmap to compare samples
}
\description{
   Draw a heatmap to compare samples.
}
\usage{
smutans.de2TypeHeatmap(x, file="default", nselect=40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
      a Smutans with an experimental design
}
  \item{file}{
      one of the types in 'x'
}
  \item{nselect}{
      one of the conditions in 'x'
}
}
\details{
   See \code{\link{smutans.de2TypeHeatmapForMatrices}} for more technical informations
}
\value{
A plot.
}
\author{
   Sang Chul Choi, schoi@cornell.edu
}
\examples{
%data( "smutansGenes" )
%x <- newSmutans(smutansGenes)
%x <- smutans.de2TypeHeatmap(x, type="ua159", condA="glucose", condB="galactose")
}
