\name{smutans.de2WithoutType}
\Rdversion{1.1}
\alias{smutans.de2WithoutType}

\title{
  Test for differences between the base means for two conditions using
  a two-factor data set
}
\description{
  This function tests for differences between the base means of two
  conditions using two-factor data sets. One factor is named type, and another
  condition. The function behaves differently depending arguments. We could
  compare two conditions given one-type. We could also compare two types given
  one-condition. We could compare two conditions considering both types. As you
  expect, we can compare two types consdering both conditions.
}
\usage{
smutans.de2WithoutType(object, condA, condB)
}
\arguments{
  \item{object}{
    a Smutans with an experimental design of two factors: one is type, and
    another is condition. See function \code{smutans.prepareDataUA159TW1} for
    example.
}
  \item{condA}{
    a string for one condition
}
  \item{condB}{
    a string for another condition
}
}
\details{
  We assume that the input data is of a two-factor. First, a two-factor data is
  created using newCountDataSet function provided by DESeq package. The
  CountDataSet class object is used to create a Smutans class object using
  newSmutans. Then, use smutans.de2WithoutType function to do a two-factor analysis to
  compare conditions using all of the data by giving the Smutans object to the
  function. The retuned object is can be used to see some results of the
  analysis.

  When type or condition is given, the function performs one-factor analysis to
  compare condA and condB.
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
%x <- smutans.de2WithoutType(x, type="ua159", condA="glucose", condB="galactose")
}
