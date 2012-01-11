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
smutans.de2List(object, file="default",qval=0.1)
}
\arguments{
  \item{object}{
    a Smutans with an experimental design of two factors: one is type, and
    another is condition. See function \code{smutans.prepareDataUA159TW1} for
    example.
}
  \item{file}{
    a file name.
}
  \item{qval}{
    a file name.
}
}
\details{
  A CSV file will be created using this function.
}
\value{
A list of gene names.
}
\author{
   Sang Chul Choi, schoi@cornell.edu
}
\examples{
%data( "smutansGenes" )
%x <- newSmutans(smutansGenes)
%x <- smutans.de2List(x, type="ua159", condA="glucose", condB="galactose")
}
