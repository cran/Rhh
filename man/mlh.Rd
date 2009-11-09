\name{mlh}
\alias{mlh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Multilocus Heterozygosity Estimates }
\description{
Read genotype data from a file, calculate homozygosity by loci, internal relatedness and stardardized heterozygosity, and save
the results to an output file.
}
\usage{
mlh(input, output, na.string, n.digits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{ Path and name of the tab or space delimited input file }
  \item{output}{ Path and name of the output file }
  \item{na.string}{ String denoting missing values in the input file }
  \item{n.digits}{ Number of decimal places to be used in the output file }  
}
\value{
Returns a three-column array with individual estimates of homozygosity by loci, internal relatedness and stardardized
heterozygosity. This data is also saved to the output file.
}
\note{
The genotype data is automatically checked for problems with \code{\link{chkdata}} when this function is used.
}
\seealso{ \code{\link{hl}}, \code{\link{ir}} and \code{\link{sh}} }
