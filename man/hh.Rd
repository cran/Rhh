\name{hh}
\alias{hh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Heterozygosity-Heterozygosity Correlation }
\description{
Get a sample of heterozygosity-heterozygosity correlations
}
\usage{
hh(genotypes, n, method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{An array with each row representing an individual and every two columns representing genotype data for two
  alleles in a loci}
  \item{n}{Number of randomizations}
  \item{method}{Multilocus heterozygosity measure to use: "hl" for homozygosity by loci, "ir" for internal relatedness and "sh"
  for stardardized heterozygosity}
}
\value{
Returns an array of heterozygosity-heterozygosity correlations computed by randomly dividing the markers into two groups and
calculating the chosen estimate of individual multilocus heterozygosity for the two set of loci.
}
\note{
The genotype data should be checked for problems with \code{\link{chkdata}} before this function is used.
}
\seealso{ \code{\link{chkdata}}, \code{\link{hl}}, \code{\link{ir}}, \code{\link{sh}}, \code{\link{mlh}} and \code{\link{h_cor}} }
