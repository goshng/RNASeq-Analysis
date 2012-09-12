\name{Smutans-class}
\Rdversion{1.1}
\docType{class}
\alias{Smutans-class}

\title{Class "Smutans" -- a container for S. mutans RNA-seq count data}
\description{	This is the main class for the present package. }
\section{Objects from the Class}{
Objects should be created with calls to \code{\link{newSmutans}} (q.v.).
}
\section{Extends}{
Class \code{CountDataSet} (package 'DESeq'), directly.
}
\note{ 
Note: This is a summary for reference. For an explanation of the actual usage, see the vignette.

An Smutans object stores RNA-seq count data from S. mutans RNA samples.

CountDataSet object created using DESeq should be given to the function
\code{\link{newSmutans}} to create Smutans object.
After creation with \code{\link{newSmutans}}, an Smutans typically contains 
a count table, i.e., a matrix of integer data, that is accessible with the accessor function \code{\link{counts}}.
Each row of the matrix corresponds to a gene (or binding region, or the like), and
each column to an experimental sample. The experimental conditions of the samples
are stored in a factor (with one element for each row of the counts matrix), which
can be read with the accessor function \code{\link{conditions}}.

}

\examples{
# See the vignette
}

