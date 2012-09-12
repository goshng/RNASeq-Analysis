\name{newSmutans}
\Rdversion{1.1}
\alias{newSmutans}

\title{
  Creates an Smutans class object
}
\description{
  This function creates an object of Smutans class type. This is the main type
  of smutans R package.
}
\usage{
newSmutans( countData, title = "default" )
}
\arguments{
  \item{countData}{
    a CoundDataSet with an experimental design of two factors: one is type, and
    another is condition. See function \code{smutans.prepareDataUA159TW1} for
    example.
}
  \item{title}{
    a string for the analysis.
}
}
\details{
  A CountDataSet object must be prepared with two factors using newCountDataSet
  function provided by R package DESeq. The title of the analysis would be
  useful to distinguish different analyses later.
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
%x <- newSmutans(x, type="ua159", condA="glucose", condB="galactose")
}
