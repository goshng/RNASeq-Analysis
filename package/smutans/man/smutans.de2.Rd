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
   condition. The function behaves differently depending arguments.
}
\usage{
smutans.de2(x,compareCondition="yes",type="",condition="",condA="",condB="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
      a Smutans with an experimental design
}
  \item{compareCondition}{
      a string for requesting comparison of conditions in a two-factor analysis.
      This is ignored when type is not empty string or condition is not empty
      string. This is used only when type, condition, condA, and condB are all
      set to an empty string. If compareCondition is yes, then we compare
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
      a string for one condition
}
  \item{condB}{
      a string for another condition
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
%x <- smutans.de2(x, type="ua159", condA="glucose", condB="galactose")
}
