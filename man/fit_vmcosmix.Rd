% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_vmcos_fns.R
\name{fit_vmcosmix}
\alias{fit_vmcosmix}
\title{Fitting bivariate von Mises cosine model mixtures using MCMC}
\usage{
fit_vmcosmix(...)
}
\arguments{
\item{...}{arguments (other than \code{model}) passed to \link{fit_angmix}}
}
\description{
Fitting bivariate von Mises cosine model mixtures using MCMC
}
\details{
Wrapper for \link{fit_angmix} with \code{model = "vmcos"}.
}
\examples{
# illustration only - more iterations needed for convergence
fit.vmcos.10 <- fit_vmcosmix(tim8, ncomp = 3, n.iter =  10,
                             n.chains = 1)
fit.vmcos.10

}
