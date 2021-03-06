\name{smutans.de2}
\Rdversion{1.1}
\alias{smutans.de2}

\title{
  Test differences between the base means for two conditions using
  a two-factor data set
}
\description{
  This function tests for differences between the base means of two
  conditions using two-factor data sets. One factor is named type, and another
  condition. The function behaves differently depending arguments. We could
  compare two conditions given one-type. We can also compare two types given
  one-condition. We can compare two conditions considering both types. As you
  expect, we can compare two types consdering both conditions.
}
\usage{
smutans.de2(object,compareCondition="yes",
            type="",condition="",condA="",condB="",
            method = "pooled")
}
\arguments{
  \item{object}{
    a Smutans with an experimental design of two factors: one is type, and
    another is condition. See function \code{newSmutans} for
    example and \code{readSmutans}.
}
  \item{compareCondition}{
    a string for requesting comparison of conditions in a two-factor analysis.
    This is ignored when type is not empty string or condition is not empty
    string. This is used only when type, condition, condA, and condB are all
    set to an empty string. If compareCondition is yes (default), then we compare
    conditions. Otherwise, the function compares types.
}
  \item{type}{
    a string for the type, if this is given, then condA and condB must be
    given because a single-factor analysis will be done
}
  \item{condition}{
    a string for the type, if this is given, then condA and condB must be
    given because a single-factor analysis will be done
}
  \item{condA}{
    [required] a string for one condition
}
  \item{condB}{
    [required] a string for another condition
}
  \item{method}{
    a string for another condition. See R package DESeq. 
}
}
\details{
  We assume that the input data is of a two-factor. First, a two-factor data is
  created using newCountDataSet function provided by DESeq package. The
  CountDataSet class object is used to create a Smutans class object using
  newSmutans. Then, use smutans.de2 function to do a two-factor analysis to
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
data( "smutansGenes" )
x <- newSmutans(smutansGenes)
x <- smutans.de2(x, type="ua159", condA="glucose", condB="galactose")
}
