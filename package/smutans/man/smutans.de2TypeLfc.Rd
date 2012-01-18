\name{smutans.de2TypeLfc}
\Rdversion{1.1}
\alias{smutans.de2TypeLfc}

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
smutans.de2TypeLfc(object, condA, condB, file="default")
}
\arguments{
  \item{object}{
    a Smutans with an experimental design of two factors: one is type, and
    another is condition. See function \code{smutans.prepareDataUA159TW1} for
    example.
}
  \item{condA}{
    a string for requesting comparison of conditions in a two-factor analysis.
    This is ignored when type is not empty string or condition is not empty
    string. This is used only when type, condition, condA, and condB are all
    set to an empty string. If compareCondition is yes (default), then we compare
    conditions. Otherwise, the function compares types.
}
  \item{condB}{
    a string for the type, if this is given, then condA and condB must be
    given because a single-factor analysis will be done
}
  \item{file}{
    a string for the type, if this is given, then condA and condB must be
    given because a single-factor analysis will be done
}
}
\details{
  We assume that the input data is of a two-factor. First, a two-factor data is
  created using newCountDataSet function provided by DESeq package. The
  CountDataSet class object is used to create a Smutans class object using
  newSmutans. Then, use smutans.de2TypeLfc function to do a two-factor analysis to
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
%x <- smutans.de2TypeLfc(x, type="ua159", condA="glucose", condB="galactose")
}
