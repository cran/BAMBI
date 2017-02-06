#' Component size selection via stepwise incremented univariate mixture models
#' @inheritParams pointest
#' @param data vector of observations. If outside, the values are transformed into the scale \eqn{[0, 2\pi)}.
#' @param model univariate angular mixture model to be fitted.  Must be one of \code{"vm"} or \code{"wnorm"}.
#' @param start_ncomp starting component size. A single component model is fitted if \code{start_ncomp} is equal to one.
#' @param max_ncomp maximum number of components allowed in the mixture model.
#' @param crit criteria for model comparison, one of \code{"AIC", "BIC", "DIC"} or \code{"WAIC"}. Default is
#' \code{"WAIC"}.
#' @param fn function to evaluate on MCMC samples to estimate parameters.
#' Defaults to \code{mean}, which computes the estimated posterior means. Ignored if \code{prev.par} is \code{FALSE}.
#' @param prev.par logical. Should the final parameters from the model with \code{ncomp = K} be used in the model
#' with \code{ncomp = K+1} as the starting parameters?
#' @param form form of crit to be used. Available choices are 1 and 2. Used only if crit is \code{"WAIC"} or
#' \code{"DIC"} and ignored otherwise.
#' @param ... additional arguments passed to \link{fit_vmmix} or \link{fit_wnormmix}, depending on \code{model}.
#' @usage
#' fit_stepwise_univariate(data, model, fn = mean,  start_ncomp = 1,
#'                         max_ncomp = 10, crit = "WAIC", form = 1,
#'                         prev.par = FALSE, burnin = 1/3, thin = 1, ...)
#'
#' @details The default \code{method} is \code{"hmc"}. Each fit uses the number of MCMC iterations specified by
#' the argument \code{n.iter} (passed to \link{fit_vmmix} or \link{fit_wnormmix}, as specified
#' through the argument \code{model}).
#'
#' @return Returns a named list with the following seven elements:
#'
#' \code{fit.all} - a list all angmcmc objects created at each component size;
#'
#' \code{fit.best} - angmcmc object corresponding to the optimum component size;
#'
#' \code{ncomp.best} - optimum component size (integer);
#'
#' \code{crit} - which model comparison criterion used (one of \code{"AIC", "BIC", "DIC"} or \code{"WAIC"});
#'
#' \code{crit.all} - all \code{crit} values calculated (for all component sizes);
#'
#' \code{crit.best} - \code{crit} value for the optimum component size; and
#'
#' \code{check_min} - logical; is the optimum component size less than \code{max_ncomp}?
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vm.step.15 <- fit_stepwise_univariate(wind, "vm", start_ncomp = 1,
#'                                           max_ncomp = 3, n.iter = 15,
#'                                           ncores = 1)
#' (fit.vm.best.15 <- bestmodel(fit.vm.step.15))
#' densityplot1d(fit.vm.best.15)
#'
#' @export


fit_stepwise_univariate <- function(data, model, fn = mean,  start_ncomp=1, max_ncomp=10,
                                    crit = "WAIC", form = 1, prev.par = FALSE, burnin = 1/3, thin = 1, ...)
{
  if(missing(model)) stop("argument \"model\" is missing, with no default")
  if(!crit %in% c("AIC", "BIC", "DIC", "WAIC")) stop("non-compatible criterion")
  if(!model %in% c("vm", "wnorm"))
    stop(paste(model, "is not a supported model. Supported models are \'vm\' and \'wnorm\' for univariate data."))

  all_ncomp <- start_ncomp:max_ncomp
  all_fit <- list()
  all_input <- list("data" = data, ...)
  all_crit <- c()
  if(!form %in% 1:2) form <- 1
  check_min <- FALSE
  for(j in seq_len(length(all_ncomp))) {
    all_input$ncomp <- all_ncomp[j]
    if(j == 1) {
      all_input$start_par <- NULL
    } else if(prev.par) {
      all_par <- pointest(all_fit[[j-1]], fn = fn, burnin = burnin, thin = thin)
      if(is.null(dim(all_par))) all_par <- t(as.matrix(all_par))
      all_par <- rbind(all_par, c(min(all_par[, "pmix"]), 1, 0))
      all_par[, "pmix"] <- all_par[, "pmix"] / sum(all_par[, "pmix"])
      start_par <- lapply(seq_len(ncol(all_par)), function(i) all_par[,i])
      names(start_par) <- all_fit[[j-1]]$par.name
      all_input$start_par <- start_par
    } else {
      all_input$start_par <- NULL
    }

    # all_input$show.progress <- FALSE

    cat("\n")

    all_fit[[j]] <- do.call(paste0("fit_", model, "mix"), all_input)

    crit_input <- list("object" = all_fit[[j]],
                       "burnin" = burnin, "thin" = thin)
    if(crit %in% c("WAIC", "DIC")) crit_input$form <- form
    all_crit[j] <- do.call(crit, crit_input)
    # if(j > start_ncomp) cat("\n")

    cat(paste("\t", "ncomp = ", all_ncomp[j], ",\t", crit, "=", round(all_crit[j], 4)))

    if(j > 1 && all_crit[j] > all_crit[j-1]) {
      check_min <- TRUE
      j <- j-1
      break
    }

  }

  result <- list("fit.all" = all_fit, "fit.best" = all_fit[[j]], "ncomp.best" = all_ncomp[j], "crit" = crit,
                 "crit.all" = all_crit, "crit.best" = all_crit[j], "check_min" = check_min)
  class(result) <- "stepfit"

  result

}



#' Component size selection via stepwise incremented bivariate mixture models
#' @inheritParams fit_stepwise_univariate
#' @inheritParams pointest
#' @param fn function to evaluate on MCMC samples to estimate parameters.
#' Defaults to \code{mean}, which computes the estimated posterior means. Ignored if \code{prev.par} is \code{FALSE}.
#' @param data two column matrix (each row being a bivariate vector) of observations. If outside, the values
#' are transformed into the scale \eqn{[0, 2\pi)}.
#' @param model bivariate angular mixture model to be fitted.  One of \code{"vmsin"}, \code{"vmcos"}, or \code{"wnorm2"}.
#' @param ... additional arguments passed to \link{fit_vmsinmix} or \link{fit_vmcosmix} or \link{fit_wnorm2mix}, depending on \code{model}.
#' @usage
#' fit_stepwise_bivariate(data, model, fn = mean,  start_ncomp = 1,
#'                        max_ncomp = 10, crit = "WAIC", form = 1,
#'                        prev.par = FALSE, burnin = 1/3, thin = 1, ...)
#'
#' @return Returns a named list with the following seven elements:
#'
#' \code{fit.all} - a list all angmcmc objects created at each component size;
#'
#' \code{fit.best} - angmcmc object corresponding to the optimum component size;
#'
#' \code{ncomp.best} - optimum component size (integer);
#'
#' \code{crit} - which model comparison criterion used (one of \code{"AIC", "BIC", "DIC"} or \code{"WAIC"});
#'
#' \code{crit.all} - all \code{crit} values calculated (for all component sizes);
#'
#' \code{crit.best} - \code{crit} value for the optimum component size; and
#'
#' \code{check_min} - logical; is the optimum component size less than \code{max_ncomp}?
#'
#' @details The default \code{method} is \code{"hmc"}. Each fit uses the number of MCMC iterations specified by
#' the argument \code{n.iter} (passed to \link{fit_vmsinmix} or \link{fit_vmcosmix} or \link{fit_wnorm2mix}, as specified
#' through the argument \code{model}).
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.step.15 <- fit_stepwise_bivariate(tim8, "vmsin", start_ncomp = 3,
#'                                             max_ncomp = 5, n.iter = 15,
#'                                             ncores = 1)
#' (fit.vmsin.best.15 <- bestmodel(fit.vmsin.step.15))
#' contour(fit.vmsin.best.15)
#'
#' @export

fit_stepwise_bivariate <- function(data, model, fn = mean, start_ncomp = 1, max_ncomp = 10,
                                   crit = "WAIC", form = 1, prev.par = FALSE, burnin = 1/3, thin = 1, ...)
{
  if(missing(model)) stop("argument \"model\" is missing, with no default")
  if(!crit %in% c("AIC", "BIC", "DIC", "WAIC")) stop("non-compatible criterion")
  if(!model %in% c("vmsin", "vmcos", "wnorm2"))
    stop(paste(model, "is not a supported model. Supported models are \'vmsin\', \'vmcos\' and \'wnorm2\' for bivariate data."))

  all_ncomp <- start_ncomp:max_ncomp
  all_fit <- list()
  all_input <- list("data" = data, ...)
  if(!form %in% 1:2) form <- 1
  all_crit <- c()
  check_min <- FALSE
  for(j in seq_len(length(all_ncomp))) {
    all_input$ncomp <- all_ncomp[j]
    if(j == 1) {
      all_input$start_par <- NULL
    } else if(prev.par){
      all_par <- pointest(all_fit[[j-1]], fn = fn, burnin = burnin, thin = thin)
      if(is.null(dim(all_par))) all_par <- t(as.matrix(all_par))
      all_par <- rbind(all_par, c(min(all_par[, "pmix"]), 1, 1, 0, 0, 0))
      all_par[, "pmix"] <- all_par[, "pmix"] / sum(all_par[, "pmix"])
      start_par <- lapply(seq_len(ncol(all_par)), function(i) all_par[,i])
      names(start_par) <- colnames(all_par)
      all_input$start_par <- start_par
    } else {
      all_input$start_par <- NULL
    }

    # all_input$show.progress <- FALSE

    cat("\n")

    all_fit[[j]] <- do.call(paste0("fit_", model, "mix"), all_input)

    crit_input <- list("object" = all_fit[[j]],
                       "burnin" = burnin, "thin" = thin)
    if(crit %in% c("WAIC", "DIC")) crit_input$form <- form
    all_crit[j] <- do.call(crit, crit_input)

    # if(j > start_ncomp) cat("\n")

    cat(paste("\t", "ncomp = ", all_ncomp[j], ",\t", crit, "=", round(all_crit[j], 4)))

    if(j > 1 && all_crit[j] > all_crit[j-1]) {
      check_min <- TRUE
      j <- j-1
      break
    }

  }

  result <- list("fit.all" = all_fit, "fit.best" = all_fit[[j]], "ncomp.best" = all_ncomp[j], "crit" = crit,
                 "crit.all" = all_crit, "crit.best" = all_crit[j], "check_min" = check_min)
  class(result) <- "stepfit"

  result
}


#' Extracting angmcmc object corresponding to the best fitted model in stepwise fits
#'
#' @param step_object stepwise fitted object (output of \code{\link{fit_stepwise_univariate}}
#' or \code{\link{fit_stepwise_bivariate}}).
#'
#' @return Returns an angmcmc object corresponding to the the best fitted model in step_object.
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.step.15 <- fit_stepwise_bivariate(tim8, "vmsin", start_ncomp = 3,
#'                                             max_ncomp = 5, n.iter = 15,
#'                                             ncores = 1)
#' fit.vmsin.best.15 <- bestmodel(fit.vmsin.step.15)
#' fit.vmsin.best.15
#'
#' @export

bestmodel <- function(step_object) {
  if(class(step_object) != "stepfit") stop("\'step_object\' is not a stepwise fitted object")
  step_object$fit.best
}
