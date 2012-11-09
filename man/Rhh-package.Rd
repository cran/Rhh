\name{Rhh-package}
\alias{Rhh-package}
\alias{Rhh}
\docType{package}
\title{
Calculating Multilocus Heterozygosity And Heterozygosity-Heterozygosity Correlation
}
\description{
Rhh is a free, open source extension package for calculating multilocus heterozygosity measures and the
heterozygosity-heterozygosity correlation in molecular ecological and evolutionary studies. 
}
\details{
Estimates of individual multilocus heterozygosity have been used as proxy measures of inbreeding, but these estimates are
expected to correlate only weakly with inbreeding coefficients calculated from the pedigree. Therefore, an inbreeding signal
of a set of microsatellite markers should be tested before using multilocus heterozygosity measures in studies of inbreeding.
This can be done by calculating the heterozygosity-heterozygosity correlation - i.e. repeatedly and randomly dividing the
loci into two groups, calculating an estimate of individual multilocus heterozygosity for the two set of loci, and obtaining
the mean correlation between two groups, which should be significant and positive if an inbreeding signal is present.

Rhh includes functions to estimate heterozygosity-heterozygosity correlation and to calculate three estimates of individual
multilocus heterozygosity: homozygosity by loci, internal relatedness, and standardized heterozygosity.

The package includes seven functions: \code{\link{chkdata}} to check data for problems; \code{\link{hl}}, \code{\link{ir}}
and \code{\link{sh}} to calculate estimates of multilocus heterozygosity; \code{\link{hh}} to get a sample of
heterozygosity-heterozygosity correlations; \code{\link{mlh}} to read data from a text file, calculate all three estimates of
multilocus heterozygosity and save them to an output file; and \code{\link{h_cor}} to read data from a text file and to get a
sample of heterozygosity-heterozygosity correlations.
}
\author{
Jussi Alho, Kaisa Välimäki

Maintainer: Jussi Alho \email{jussi.s.alho@iki.fi}

\url{http://www.helsinki.fi/biosci/egru/research/software}
}

\seealso{ \code{\link{chkdata}}, \code{\link{hl}}, \code{\link{ir}}, \code{\link{sh}}, \code{\link{hh}}, \code{\link{mlh}} and \code{\link{h_cor}}}
