% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_model_select.R
\name{bestmodel}
\alias{bestmodel}
\alias{bestcriterion}
\title{Convenience function for extracting angmcmc object, and the value of the model
selection criterion corresponding to the best fitted model in stepwise fits}
\usage{
bestmodel(step_object)

bestcriterion(step_object)
}
\arguments{
\item{step_object}{stepwise fitted object obtained from \link{fit_incremental_angmix}.}
}
\value{
\code{bestmodel} returns an \code{angmcmc} object, and
\code{bestcriterion} returns the  corresponding value of model selection criterion  for the best fitted model in \code{step_object}.
}
\description{
Convenience function for extracting angmcmc object, and the value of the model
selection criterion corresponding to the best fitted model in stepwise fits
}
\details{
These are convenience functions; the best fitted model and the corresponding value of model selection criterion
can also be directly obtained by
extracting the elements \code{"fit.best"} and \code{"crit.best"} from \code{step_object} respectively.
Note that \verb{bestcriterion\} returns: (a) a scalar number (class = }numeric\verb{) if }crit\verb{used in original}fit_incremental_angmix\verb{call is}'AIC'\verb{, }'BIC'\code{or}'DIC'\verb{, (b) an element of class }bridge\verb{from package}bridgesampling\code{if}crit\code{is}LOGML\verb{, (c) an element of class }c("waic", "loo")\code{if}crit = 'WAIC'\verb{, and (d) an element of class }c("psis_loo", "loo")\code{if}crit = "LOOIC"`. See documentations of these model
selection criteria for more details.
}
\examples{
# illustration only - more iterations needed for convergence
set.seed(1)
fit.vmsin.step.15 <- fit_incremental_angmix("vmsin", tim8, start_ncomp = 1,
                                            max_ncomp = 3, n.iter = 15,
                                            n.chains = 1,
                                            crit = "WAIC")
fit.vmsin.best.15 <- bestmodel(fit.vmsin.step.15)
fit.vmsin.best.15

crit.best <- bestcriterion(fit.vmsin.step.15)
crit.best

}
