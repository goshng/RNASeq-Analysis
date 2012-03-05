\name{smutans.de2List}
\Rdversion{1.1}
\alias{smutans.de2List}

\title{
  Writes an output file in CSV format
}
\description{
  Users can use a spreadsheet program to open a CSV format file.
}
\usage{
smutans.de2List(object,file="",qval=0.1)
}
\arguments{
  \item{object}{
    a Smutans with an experimental design of two factors: one is type, and
    another is condition. See function \code{newSmutans} for
    example and \code{readSmutans}.
}
  \item{file}{
    a file name. When file is given, it prints values to a CSV file. Otherwise,
    it prints values to a standard output in R.
}
  \item{qval}{
    a numeric value of cutoff when printing to a standard output in R. This is
    not used when printing values in a file.
}
}
\value{
A list of gene names or a CSV file.
}
\author{
   Sang Chul Choi, schoi@cornell.edu
}
\examples{
data( "smutansGenes" )
x <- newSmutans(smutansGenes)
x <- smutans.de2(x, type="ua159", condA="glucose", condB="galactose")
smutans.de2List( x, "x.csv" )
}
