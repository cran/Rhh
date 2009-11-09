\name{chkdata}
\alias{chkdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Check Genotype Data}
\description{
Check genotype data for problems.
}
\usage{
chkdata(genotypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{An array with each row representing an individual and every two columns representing genotype data for two
  alleles in a loci.}
}
\details{
This function checks that 1) there is even number of columns, i.e. two alleles per individual for each locus, in the genotype
data; 2) each loci is polymorphic, i.e. has more than one allele in the data; and 3) that each individual has either two null
alleles or two non-null alleles in a given locus.

Error occurs if the first two checks fail. If the last check fails, there will be a warning.
}
\value{
None
}
\note{
This function is meant to be run with genotype data before using functions \code{\link{hl}}, \code{\link{ir}}, \code{\link{sh}}
and \code{\link{hh}}. It is called automatically by functions \code{\link{mlh}} and \code{\link{h_cor}}.
}
\seealso{ \code{\link{hl}}, \code{\link{ir}}, \code{\link{sh}}, \code{\link{hh}}, \code{\link{mlh}} and \code{\link{h_cor}}}
