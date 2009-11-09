\name{sh}
\alias{sh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standardized Heterozygosity}
\description{
Calculate standardized heterozygosity.
}
\usage{
sh(genotypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{An array with each row representing an individual and every two columns representing genotype data for two
  alleles in a loci.}
}
\details{
Calculates standardized heterozygosity, a multilocus heterozygosity measure developed by Coltman et al. (1999), for each
individual in the genotype data.
}
\value{
An array with standardized heterozygosity estimate for each individual. If the estimate could not be calculated for a given
individual, the corresponding value is \code{\link{NaN}}.
}
\note{
The genotype data should be checked for problems with \code{\link{chkdata}} before this function is used.
}
\seealso{\code{\link{hl}}, \code{\link{ir}}, \code{\link{chkdata}} and \code{\link{mlh}}}
\references{
Coltman DW, Pilkington JG, Smith JA et al (1999) Parasite-mediated selection against inbred Soay sheep in a free-living island
population. Evolution 53:1259-1267
}
