% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_wnorm_fns.R
\name{rwnorm}
\alias{rwnorm}
\alias{dwnorm}
\title{The univariate Wrapped Normal distribution}
\usage{
rwnorm(n = 1, kappa = 1, mu = 0)

dwnorm(x, kappa = 1, mu = 0, int.displ, log = FALSE)
}
\arguments{
\item{n}{number of observations. Ignored if at least one of the other parameters have length k > 1, in which
case, all the parameters are recycled to length k to produce k random variates.}

\item{kappa}{vector of concentration (inverse-variance) parameters; \code{kappa} > 0.}

\item{mu}{vector of means.}

\item{x}{vector of angles (in radians) where the densities are to be evaluated.}

\item{int.displ}{integer displacement. If \code{int.displ =} M, then the infinite sum in the
density is approximated by a sum over 2*M + 1 elements. (See Details.) The allowed values are 1, 2, 3, 4 and 5. Default is 3.}

\item{log}{logical. Should the log density be returned instead?}
}
\value{
\code{dwnorm} gives the density  and \code{rwnorm} generates random deviates.
}
\description{
The univariate Wrapped Normal distribution
}
\details{
If \code{mu} and \code{kappa} are not specified they assume the default values of \code{0} and \code{1} respectively.

The univariate wrapped normal distribution has density
\deqn{f(x) = \sqrt(\kappa/(2\pi)) \sum \exp(-\kappa/2 (x - \mu(2\pi\omega))^2)}
where the sum extends over all integers \eqn{\omega},
and is approximated by a sum over \eqn{\omega} in \eqn{\{-M, -M+1, ..., M-1, M \}} if \code{int.displ = } \eqn{M}.
}
\examples{

kappa <- 1:3
mu <- 0:2
x <- 1:10
n <- 10


# when x and both parameters are scalars, dwnorm returns a single density
dwnorm(x[1], kappa[1], mu[1])

# when x is a vector but both the parameters are scalars, dmv returns a vector of
# densities calculated at each entry of x with the same parameters
dwnorm(x, kappa[1], mu[1])

# if x is scalar and at least one of the two paraemters is a vector, both parameters are
# recycled to the same length, and dwnorm returns a vector of with ith element being the
# density evaluated at x with parameter values kappa[i] and mu[i]
dwnorm(x[1], kappa, mu)

# if x and at least one of the two paraemters is a vector, x and the two parameters are
# recycled to the same length, and dwnorm returns a vector of with ith element being the
# density at ith element of the (recycled) x with parameter values kappa[i] and mu[i]
dwnorm(x, kappa, mu)

# when parameters are all scalars, number of observations generated by rwnorm is n
rwnorm(n, kappa[1], mu[1])

# when at least one of the two parameters is a vector, both are recycled to the same length,
# n is ignored, and the number of observations generated by rwnorm is the same as the length
# of the recycled vectors
rwnorm(n, kappa, mu)

}
