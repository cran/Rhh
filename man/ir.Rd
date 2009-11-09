\name{ir}
\alias{ir}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal Relatedness}
\description{
Calculate internal relatedness.
}
\usage{
ir(genotypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{An array with each row representing an individual and every two columns representing genotype data for two
  alleles in a loci.}
}
\details{
Calculates internal relatedness, a multilocus heterozygosity measure developed by Amos et al. (2001), for each individual in the
genotype data.
}
\value{
An array with internal relatedness estimate for each individual. If the estimate could not be calculated for a given individual,
the corresponding value is \code{\link{NaN}}.
}
\note{
The genotype data should be checked for problems with \code{\link{chkdata}} before this function is used.
}
\seealso{\code{\link{hl}}, \code{\link{sh}}, \code{\link{chkdata}} and \code{\link{mlh}}}
\references{
Amos W, Worthington Wilmer J, Fullard K et al (2001) The influence of parental relatedness on reproductive success. Proc R Soc
Lond B 268:2021-2027
}
