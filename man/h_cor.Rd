\name{h_cor}
\alias{h_cor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Heterozygosity-Heterozygosity Correlation }
\description{
Read genotype data from a file and get a sample of heterozygosity-heterozygosity correlations
}
\usage{
h_cor(input, na.string, n, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{ Path and name of the tab or space delimited input file }
  \item{na.string}{ String denoting missing values in the input file }
  \item{n}{ Number of randomizations }
  \item{method}{ Multilocus heterozygosity measure to use: "hl" for homozygosity by loci, "ir" for internal relatedness and "sh"
  for stardardized heterozygosity }
}
\value{
Returns an array of heterozygosity-heterozygosity correlations computed by randomly dividing the markers into two groups and
calculating the chosen estimate of individual multilocus heterozygosity for the two set of loci.
}
\note{
The genotype data is automatically checked for problems with \code{\link{chkdata}} when this function is used.
}
\seealso{ \code{\link{hh}} }
