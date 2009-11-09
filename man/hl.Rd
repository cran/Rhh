\name{hl}
\alias{hl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Homozygosity By Loci}
\description{
Calculate homozygosity by loci.
}
\usage{
hl(genotypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{An array with each row representing an individual and every two columns representing genotype data for two
  alleles in a loci.}
}
\details{
Calculates homozygosity by loci, a multilocus heterozygosity measure developed by Aparicio et al. (2007), for each individual in
the genotype data.
}
\value{
An array with homozygosity by loci estimate for each individual. If the estimate could not be calculated for a given individual,
the corresponding value is \code{\link{NaN}}.
}
\note{
The genotype data should be checked for problems with \code{\link{chkdata}} before this function is used.
}
\seealso{\code{\link{ir}}, \code{\link{sh}}, \code{\link{chkdata}} and \code{\link{mlh}}}
\references{
Aparicio JM, Ortego J, Cordero PJ (2006) What should we weigh to estimate heterozygosity, alleles or loci? Mol Ecol 15:4659-4665
}
