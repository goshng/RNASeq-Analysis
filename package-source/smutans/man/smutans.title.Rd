\name{smutans.title}
\Rdversion{1.1}
\alias{smutans.title}

\title{
  Returns a list of adjusted p-values
}
\description{
  A list of multiple testing adjusted p-values.
}
\usage{
smutans.title(object)
}
\arguments{
  \item{object}{
    a Smutans with an experimental design of two factors: one is type, and
    another is condition. See function \code{smutans.prepareDataUA159TW1} for
    example.
}
}
\details{
  A list of p-values.
}
\value{
A list of p-values.
}
\author{
   Sang Chul Choi, schoi@cornell.edu
}
\examples{
%data( "smutansGenes" )
%x <- newSmutans(smutansGenes)
%x <- smutans.title(x, type="ua159", condA="glucose", condB="galactose")
}
