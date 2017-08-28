#' Fix label switching in angmcmc objects
#' @inheritParams pointest
#' @param method method to use for fixing label switching. Available choices are 1 for \link{dataBased}
#' and 2 for \link{ecr.iterative.1}.
#'
#' @details \code{fix_label} is a wrapper for \link{dataBased} or \link{ecr.iterative.1} (depending on
#' \code{method}) from the \code{label.switching} package for \code{angmcmc} objects.
#'
#' @return Returns another \code{angmcmc} object, with the parameter values (after burn-in and thin)
#' re-ordered according to the resulting permutation from \link{dataBased} or \link{ecr.iterative.1}.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' # now apply fix_label
#' fit.vmsin.20.fix <- fix_label(fit.vmsin.20)
#'
#' @export

fix_label <- function(object, burnin = 1/3, thin = 1, method = 1) {
  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  if(class(object) != "angmcmc") stop("\"object\" must be an angmcmc object")

  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  lab_all <- matrix(1:object$ncomp, nrow = object$n.iter+1, ncol = object$ncomp, byrow = TRUE)
  if(method == 1) {
    lab_switch <- dataBased(x = object$data, K = object$ncomp, z = t(object$clus.ind[ , final_iter]))$per
  } else {
    lab_switch <- ecr.iterative.1(z = t(object$clus.ind[ , final_iter]), K = object$ncomp)$per
  }
  lab_all[final_iter, ] <- lab_switch
  par_val_adj <- par_mat_permute(object$par.value, lab_all)
  dimnames(par_val_adj) <- dimnames(object$par.value)

  object_fixed <- object
  object_fixed$par.value <- par_val_adj
  object_fixed$fixed.label <- TRUE
  object_fixed$burnin <- burnin
  object_fixed$thin <- thin
  object_fixed
}



#' Point estimates for parameters from an angmcmc object
#' @param object angular MCMC object.
#' @param fn function to evaluate on MCMC samples to estimate parameters.  Defaults to \code{mean}, which computes the estimated posterior mean.
#' @param par.name vector of names of parameters for which point estimates are to be computed.  If \code{NULL}, results for all parameters are provided.
#' @param comp.label vector of component labels (positive integers, e.g., \code{1, 2, ...}) for which point estimates are to be computed.
#' If \code{NULL}, results for all components are provided.
#' @param burnin initial fraction of the MCMC samples to be discarded as burn-in. Must be a value in [0, 1).
#' @param thin positive integer. If \code{thin =} \eqn{n}, only every \eqn{n}-th realizations of the Markov chain is kept.
#'
#' @return Returns a matrix of point estimates, or vector of point estimates if \code{length(par.name)==1} or \code{length(comp.label)==1}.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' # estimate parameters by sample mean
#' (est_mean <- pointest(fit.vmsin.20))
#' # estimate parameters by sample median
#' (est_median <- pointest(fit.vmsin.20, fn = median))
#' @export

pointest <- function(object, fn = mean, par.name = NULL,
                     comp.label = NULL, burnin = 1/3, thin = 1)
{
  if(class(object) != "angmcmc") stop("\'object\' must be an angmcmc object")
  fn <- match.fun(fn)

  if(!is.null(object$fixed.label)) {
    if(missing(burnin)) burnin <- object$burnin
    if(missing(thin)) thin <- object$thin
  }

  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")
  if(missing(par.name)) par.name <- object$par.name
  else if(!all(par.name %in% object$par.name)) stop("invalid parameter name")


  if(missing(comp.label)) comp.label <- 1:object$ncomp
  else if(!all(comp.label %in% 1:object$ncomp)) stop("invalid component label")

  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  if(all(c(length(comp.label), length(par.name)) == 1 )) {
    res <- fn(object$par.value[par.name, comp.label, final_iter])
  } else if(all(c(length(comp.label), length(par.name)) > 1 )) {
    if(identical(fn, mean)) {
      res <- sapply(par.name, function(name) rowMeans(object$par.value[name, comp.label, final_iter]))
    } else {
      res <- sapply(par.name, function(name) apply(object$par.value[name, comp.label, final_iter], 1, fn))
    }
  } else  {
    res <- apply(object$par.value[par.name, comp.label, final_iter], 1, fn)
  }

  res
}


#' Quantile estimates for parameters from an angmcmc object
#' @inheritParams pointest
#' @param x angmcmc object
#' @inheritParams stats::quantile
#' @return Returns a three dimensional array of quantiles, or a matrix (vector) of quantiles
#'  if one (or two) among \code{par.name},  \code{comp.label}, \code{probs} has length 1.
#' @param ... further arguments to pass to \code{quantile}.  In particular, \code{probs = seq(0, 1, 0.25)}
#' is the default vector of quantiles computed for each parameter.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' # 0.025th quantiles
#' (quant_025 <- quantile(fit.vmsin.20, prob = 0.025))
#' # 0.975th quantiles
#' (quant_975 <- quantile(fit.vmsin.20, prob = 0.975))
#' # default quantiles
#' (quant_def <- quantile(fit.vmsin.20))
#'
#' @export


quantile.angmcmc <- function(x, par.name = NULL, comp.label=NULL,
                             burnin = 1/3, thin = 1, ...)
{
  if(!is.null(x$fixed.label)) {
    if(missing(burnin)) burnin <- x$burnin
    if(missing(thin)) thin <- x$thin
  }

  if(missing(par.name)) par.name <- x$par.name
  else if(!all(par.name %in% x$par.name)) stop("invalid parameter name")
  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  if(missing(comp.label)) comp.label <- 1:x$ncomp
  else if(!all(comp.label %in% 1:x$ncomp)) stop("invalid component label")


  burnin_iter <- 1:floor((x$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:x$n.iter)[-burnin_iter][thin_filter]

  if(all(c(length(comp.label), length(par.name)) == 1 )) {
    res <- quantile(x$par.value[par.name, comp.label, final_iter], ...)
  } else if(all(c(length(comp.label), length(par.name)) > 1 ) ) {
    res <- apply(x$par.value[par.name, comp.label, final_iter], c(2, 1), function(x) quantile(x, ...))
  } else {
    res <- apply(x$par.value[par.name, comp.label, final_iter], 1, function(x) quantile(x, ...))
  }

  res
}


#' Extract MCMC samples for parameters from an angmcmc object
#' @inheritParams pointest
#' @param object angular MCMC object
#' @details The default for both \code{par.name} and \code{comp.label} are the all possible choices
#' available in \code{object}.
#' @return
#' Returns a  matrix (vector) if one (both) of \code{par.name} and \code{comp.label} is of length 1, and a three dimensional array with
#' third dimension corresponding to iterations if both \code{par.name} and \code{comp.label} have length \code{>1}.
#'
#' @usage
#' extractsamples(object, par.name = NULL, comp.label = NULL,
#'                burnin = 1/3, thin = 1)
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' # extract Markov chain realizations for kappa1 from component 1
#' extr_kappa1_1 <- extractsamples(fit.vmsin.20, "kappa1", 1)
#' # for kappa1 from component from all components
#' extr_kappa1 <- extractsamples(fit.vmsin.20, "kappa1")
#' # for all parameters in component 1
#' extr_1 <- extractsamples(fit.vmsin.20, comp.label = 1)
#'
#' @export

extractsamples <- function(object, par.name = NULL, comp.label = NULL,
                           burnin = 1/3, thin = 1)
{
  if(class(object) != "angmcmc") stop("object must be an angmcmc object")

  if(missing(par.name)) par.name <- object$par.name
  else if(!all(par.name %in% object$par.name)) stop("invalid parameter name")

  if(missing(comp.label)) comp.label <- 1:object$ncomp
  else if(!all(comp.label %in% 1:object$ncomp)) stop("invalid component label")

  if(!is.null(object$fixed.label)) {
    if(missing(burnin)) burnin <- object$burnin
    if(missing(thin)) thin <- object$thin
  }

  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")


  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  object$par.value[par.name, comp.label, final_iter]
}


#' Summary statistics for parameters from an angmcmc object
#' @inheritParams base::summary
#' @inheritParams pointest
#' @param object angular MCMC object.
#' @details Computes (after thinning and discarding burn-in) point estimates with 95\% posterior credible sets for all components and all parameters,
#' together with the sample averages of log likelihood and log posterior density.
#' @return Returns a list with elements \code{estimate, lower, upper, llik} and \code{lpd}.
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' summary(fit.vmsin.20)
#'
#' @export

summary.angmcmc <- function(object, burnin = 1/3, thin = 1, ...)
{
  if(!is.null(object$fixed.label)) {
    if(missing(burnin)) burnin <- object$burnin
    if(missing(thin)) thin <- object$thin
  }
  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")


  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  if(object$ncomp == 1) {
    est <- rowMeans(object$par.value[, , final_iter])[-1]
    upp <- apply(object$par.value[ , , final_iter], 1, function(x) quantile(x, probs = 0.975))[-1]
    low <- apply(object$par.value[ , , final_iter], 1, function(x) quantile(x, probs = 0.025))[-1]

  } else {
    est <- sapply(object$par.name, function(name) rowMeans(object$par.value[name, , final_iter]))
    upp <- apply(object$par.value[ , , final_iter], c(2, 1), function(x) quantile(x, probs = 0.975))
    low <- apply(object$par.value[ , , final_iter], c(2, 1), function(x) quantile(x, probs = 0.025))
  }


  llik <- mean(object$llik[final_iter])
  lpd <- mean(object$lpd[final_iter])

  res <- list("estimate" = est, "upper" = upp, "lower" = low,
              "llik" = llik, "lpd" = lpd)

  res
}


#' AIC and BIC for angmcmc objects
#' @inheritParams stats::AIC
#' @param object an angmcmc object.
#' @inheritParams summary.angmcmc
#'
#' @return AIC computes the AIC and BIC computes BIC for \code{angmcmc} objects.
#'
#' @details
#' Let \eqn{\hat{L}} be the maximum value of the likelihood function for the model, \eqn{m} be the number of
#' estimated parameters in the model and \eqn{n} be the number of data points. Then AIC and BIC are defined as
#' \eqn{AIC = -2 \log \hat{L} + mk} and \eqn{BIC = -2 \log \hat{L} + m \log(n)}.
#'
#' \eqn{\hat{L}} is estimated by the sample maximum obtained from the MCMC realizations.
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' AIC(fit.vmsin.20)
#' BIC(fit.vmsin.20)
#'
#' @export

AIC.angmcmc <- function(object, burnin = 1/3, thin = 1, k = 2, ...)
{
  if(!is.null(object$fixed.label)) {
    if(missing(burnin)) burnin <- object$burnin
    if(missing(thin)) thin <- object$thin
  }

  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  llik <- max(object$llik[final_iter])
  npar <- prod(dim(object$par.value)[1:2]) - 1
  aic <- k * npar - 2 * llik
  aic
}


#' @rdname AIC.angmcmc
#' @export

BIC.angmcmc <- function(object, burnin = 1/3, thin = 1, ...)
{
  if(!is.null(object$fixed.label)) {
    if(missing(burnin)) burnin <- object$burnin
    if(missing(thin)) thin <- object$thin
  }

  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  llik <- max(object$llik[final_iter])
  npar <- prod(dim(object$par.value)[1:2]) - 1
  bic <- log(object$n.data) * npar - 2 * llik
  bic
}



#' Deviance Information Criterion (DIC) for angmcmc objects
#' @inheritParams WAIC
#' @inheritParams pointest
#' @param ... additional model specific arguments to be passed to \code{DIC}. For example, \code{int.displ}
#' specifies integer dispacement in wnorm and wnorm2 models. See \link{fit_wnormmix} and
#' \link{fit_wnorm2mix} for more details.
#' @param form form of DIC to use. Available choices are 1 (default) and 2. See details.
#' @return Computes the DIC for a given angmcmc object
#' @details Given a deviance function \eqn{D(\theta) = -2 log(p(y|\theta))}, and an estimate
#' \eqn{\theta* = (\sum \theta_i) / N} of the posterior mean
#' \eqn{E(\theta|y)}, where \eqn{y} denote the data, \eqn{\theta} are the unknown
#' parameters of the model, \eqn{\theta_1, ..., \theta_N} are MCMC samples from the posterior
#' distribution of \eqn{\theta} given \eqn{y} and \eqn{p(y|\theta)} is the likelihood function,
#' the (form 1 of) Deviance Infomation Criterion (DIC) is defined as
#' \deqn{DIC = 2 ( (\sum_{s=1}^N D(\theta_s)) / N - D(\theta*) )}
#' The second form for DIC is given by
#' \deqn{DIC = D(\theta*) - 4 \hat{var} \log p(y|\theta_s)}
#' where for \eqn{i = 1, ..., n}, \eqn{\hat{var} \log p(y|\theta)} denotes the estimated variance
#' of the log likelihood based on the realizations \eqn{\theta_1, ..., \theta_N}.
#'
#' Like AIC and BIC, DIC is an asymptotic approximation for large samples, and
#' is only valid when the posterior distribution is approximately normal.
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' DIC(fit.vmsin.20)
#'
#' @export

DIC <- function(object, form = 1, burnin = 1/3, thin = 1, ...)
{
  if(!is.null(object$fixed.label)) {
    if(missing(burnin)) burnin <- object$burnin
    if(missing(thin)) thin <- object$thin
  }

  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  if(!form %in% 1:2) form <- 1

  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  all_llik <- object$llik[final_iter]

  par_est <- pointest(object, fn = mean, burnin = burnin, thin = thin)

  if(object$ncomp > 1) {
    all_input <- list("par" = t(par_est[,-1]), "pi" = par_est[,1])
  } else{
    all_input <- list("par" = matrix(par_est[-1], ncol = 1), "pi" = par_est[1])
  }

  ell <- list(...)

  if(object$model == "wnorm") {
    if(is.null(ell$int.displ) || ell$int.displ < 1 || ell$int.displ > 5) int.displ <- 3
    displ <- as.integer(int.displ)
    all_input$omega_2pi_1d <- omega.2pi.1d <- (-displ):displ * (2*pi) # 2pi * 1d integer displacements
  } else if(object$model == "wnorm2") {
    if(is.null(ell$int.displ) || ell$int.displ < 1 || ell$int.displ > 5) int.displ <- 3
    displ <- as.integer(int.displ)
    all_input$omega_2pi <- as.matrix(expand.grid(-displ:displ,-displ:displ) * (2*pi))
  }

  all_input$data <- object$data

  modelname <- object$model
  if(object$type == "uni") modelname <- paste0("uni", modelname)

  const_calc_input <- list("par_mat" = all_input$par)
  if(object$model == "vmcos") const_calc_input$uni_rand <- matrix(runif(2e4), ncol = 2)
  all_input$log_c <- do.call(paste0("log_const_", modelname, "_all"), const_calc_input)

  llik_estpar <- do.call(paste0("llik_", modelname, "_full"), all_input)

  if(form == 1) {
    D_bar <- -2 * mean(all_llik)
    D_estpar <- -2 * llik_estpar
    dic <- 2 * D_bar - D_estpar
  } else {
    dic <- -2 * llik_estpar + 4 * var(all_llik)
  }
  dic
}



#' Watanabe-Akaike Information Criterion (WAIC) for angmcmc objects
#' @inheritParams pointest
#' @param form the form of p_W to use. Available choices are 1 (default) and 2. See details.
#' @param ... additional model specific arguments to be passed to \code{WAIC}. For example, \code{int.displ}
#' specifies integer dispacement in wnorm and wnorm2 models. See \link{fit_wnormmix} and
#' \link{fit_wnorm2mix} for more details.
#' @return Computes the WAIC for a given angmcmc object.
#' @details
#' Given a deviance function \eqn{D(\eta) = -2 \log(p(y|\eta))}, and an estimate
#' \eqn{\eta* = (\sum \eta_i) / n} of the posterior mean
#' \eqn{E(\eta|y)}, where \eqn{y = (y_1, ..., y_n)} denote the data, \eqn{\eta} is the unknown
#' parameter vector of the model, \eqn{\eta_1, ..., \eta_N} are MCMC samples from the posterior
#' distribution of \eqn{\eta} given \eqn{y} and \eqn{p(y|\eta)} is the likelihood function,
#' the Watanabe-Akaike Information Criterion (WAIC) is defined as
#' \deqn{WAIC = LPPD - p_W}
#' where
#' \deqn{LPPD  = \sum_{i=1}^n \log (N^{-1} \sum_{s=1}^N p(y_i|\eta_s) )}
#' and (form 1 of)
#' \deqn{p_W =  2 \sum_{i=1}^n [ \log (N^{-1} \sum_{s=1}^N p(y_i|\eta_s) ) - N^{-1} \sum_{s=1}^N \log \:p(y_i|\eta_s) ].}
#' An alternative form (form 2) for \eqn{p_W} is given by
#' \deqn{p_W = \sum_{i=1}^n \hat{var} \log p(y_i|\eta)}
#' where for \eqn{i = 1, ..., n}, \eqn{\hat{var} \log p(y_i|\eta)} denotes the estimated variance
#' of \eqn{\log p(y_i|\eta)} based on the realizations \eqn{\eta_1, ..., \eta_N}.
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' WAIC(fit.vmsin.20)
#'
#' @export


WAIC <- function(object, form = 1, burnin = 1/3, thin = 1, ...)
{
  if(missing(burnin)) {
    if(is.null(object$fixed.label)) burnin <- 1/3
    else burnin <- object$burnin
  }

  if(missing(thin)) {
    if(is.null(object$fixed.label)) thin <- 1
    else thin <- object$thin
  }

  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  if(!form %in% 1:2) form <- 1

  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  ncomp <- object$ncomp

  ell <- list(...)

  all_par <- extractsamples(object, burnin = burnin, thin = thin)

  if(!object$model %in% c("wnorm", "wnorm2")) {
    if(ncomp > 1) {
      all_den_mat <- suppressWarnings(sapply(1:length(final_iter),
                                             function(iter) do.call(paste0("d", object$model, "mix"),
                                                                    addtolist(list_by_row(all_par[ , , iter], 1:ncomp),
                                                                              x = object$data)) ))
    } else {
      all_den_mat <- suppressWarnings(sapply(1:length(final_iter),
                                             function(iter) do.call(paste0("d", object$model),
                                                                    addtolist(as.list(all_par[-1 , iter]), x = object$data) )))
    }

  } else {

    int.displ <- ell$int.displ
    if(is.null(int.displ) || int.displ < 1 || int.displ > 5) int.displ <- object$int.displ
    int.displ <- as.integer(int.displ)

    if(ncomp > 1) {
      all_den_mat <- (sapply(1:length(final_iter),
                             function(iter) do.call(paste0("d", object$model, "mix"),
                                                    addtolist(list_by_row(all_par[ , , iter], 1:ncomp),
                                                              x = object$data, int.displ = int.displ) )))
    } else {
      all_den_mat <- (sapply(1:length(final_iter),
                             function(iter) do.call(paste0("d", object$model),
                                                    addtolist(as.list(all_par[-1 , iter]),
                                                              x = object$data, int.displ = int.displ) )))
    }

  }

  lppd <- sum(log(rowMeans(all_den_mat)))
  if(form == 1) {
    p_waic <- 2 * (lppd - sum(rowMeans(log(all_den_mat))))
  } else {
    p_waic <- sum(rowVars(log(all_den_mat)))
  }

  waic <- 2 * (p_waic - lppd)

  waic
}





#' Density and random deviates from an angmcmc object
#' @inheritParams pointest
#' @param object angular MCMC object. The dimension of the model must match with \code{x}.
#' @param x vector (if univariate) or a two column matrix (if bivariate, with each row a 2-D vector) of points where the
#' densities are to be computed.
#' @param n number of observations to be generated.
#' @return \code{d_fitted} gives a vector the densities computed at the given points  and \code{r_fitted}
#' creates a vector (if univariate) or a matrix (if bivariate) with each row being a 2-D point, of random deviates.
#'
#' @details
#' To estimate the mixture density, first the parameter vector \eqn{\eta} is estimated
#' by applying \code{fn} on the MCMC samples (using the function \link{pointest}), yielding the (consistent) Bayes estimate \eqn{\hat{\eta}}. Then the mixture density
#' \eqn{f(x|\eta)} at any point \eqn{x} is (consistently) estimated by \eqn{f(x|\hat{\eta})}.
#'
#' The random deviates are generated from the estimated mixture density \eqn{f(x|\hat{\eta})}.
#'
#'
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' d_fitted(c(0,0), fit.vmsin.20)
#' r_fitted(10, fit.vmsin.20)
#' @export


d_fitted <- function(x, object, fn = mean, burnin = 1/3, thin = 1)
{
  if(class(object) != "angmcmc") stop("object must be an angmcmc object")
  if(object$type == "bi") {
    if((length(dim(x)) < 2 && length(x) != 2) || (length(dim(x)) == 2 && tail(dim(x), 1) != 2)
       || (length(dim(x)) > 2)) stop("x must either be a bivariate vector or a two-column matrix")
  }
  if(missing(burnin)) {
    if(is.null(object$fixed.label)) burnin <- 1/3
    else burnin <- object$burnin
  }

  if(missing(thin)) {
    if(is.null(object$fixed.label)) thin <- 1
    else thin <- object$thin
  }

  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  if(object$ncomp == 1) {
    est <- pointest(object, burnin = burnin, thin = thin, fn = fn)
    allinput <- as.list(est)
    allinput$pmix <- NULL
    allinput$x <- x
    suppressWarnings(do.call(what = paste0("d", object$model), args = allinput))

  } else {
    est <- pointest(object, burnin = burnin, thin = thin, fn = fn)
    allinput <- lapply(seq_len(ncol(est)), function(i) est[,i])
    names(allinput) <- object$par.name
    allinput$x <- x
    suppressWarnings(do.call(what = paste0("d", object$model, "mix"), args = allinput))
  }
}



#' @rdname d_fitted
#' @export
r_fitted <- function(n, object, fn = mean, burnin = 1/3, thin = 1)
{
  if(class(object) != "angmcmc") stop("\"object\" must be an angmcmc object")

  if(!is.null(object$fixed.label)) {
    if(missing(burnin)) burnin <- object$burnin
    if(missing(thin)) thin <- object$thin
  }

  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  if(object$ncomp == 1) {
    est <- pointest(object, burnin = burnin, thin = thin, fn = fn)
    allinput <- as.list(est)
    allinput$pmix <- NULL
    allinput$n <- n
    suppressWarnings(do.call(what = paste0("r", object$model), args = allinput))

  } else {
    est <- pointest(object, burnin = burnin, thin = thin, fn = fn)
    allinput <- lapply(seq_len(ncol(est)), function(i) est[,i])
    names(allinput) <- object$par.name
    allinput$n <- n
    suppressWarnings(do.call(what = paste0("r", object$model, "mix"), args = allinput))
  }

}


#' Extract Log-Likelihood from angmcmc objects
#' @inheritParams pointest
#' @param fn function to evaluate on MCMC samples to estimate parameters.  Defaults to \code{mean}, which computes the estimated posterior mean. Used for parameter estimation.
#' @param burnin initial fraction of the MCMC samples to be discarded as burn-in. Must be a value in [0, 1). Used for parameter estimation.
#' @param thin positive integer. If \code{thin =} \eqn{n}, only every \eqn{n}-th realizations of the Markov chain is kept. Used for parameter estimation.
#' @inheritParams stats::logLik
#'
#' @details In order to estimate the log likelihood for the model, first the parameter vector is estimated using \link{pointest},
#' and then log the likelihood is calculated on the basis of the estimated parameter.
#'
#' The degrees of the likelihood function is the total number of free parameters estimated in the mixture models,
#' which is equal to \eqn{6K - 1} for bivariate models (vmsin, vmcos and wnorm2), or \eqn{3K - 1} for univariate
#' models (vm and wnorm), where \eqn{K} denotes the number of components in the mixture model.
#'
#' @return Returns an object of class \link{logLik}. This is a number (the estimated log likelihood) with attributes "df"
#' (degrees of freedom) and "nobs" (number of observations).
#'
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' logLik(fit.vmsin.20)
#' @export


logLik.angmcmc <- function(object, fn = mean, burnin = 1/3, thin = 1, ...)
{
  if(class(object) != "angmcmc") stop("\"object\" must be an angmcmc object")
  if(!is.null(object$fixed.label)) {
    if(missing(burnin)) burnin <- object$burnin
    if(missing(thin)) thin <- object$thin
  }
  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  all_den <- d_fitted(object$data, object, fn = fn, burnin = burnin, thin = thin)
  llik <- sum(log(all_den))

  if(object$type == "uni") df <- 3*object$ncomp - 1
  else df <- 6*object$ncomp - 1

  object_nobs <- object$n.data
  result <- llik
  attributes(result) <- list("df" = df, "nobs" = object_nobs)
  class(result) <- "logLik"
  result
}
