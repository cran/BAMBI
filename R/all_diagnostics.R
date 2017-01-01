#' Contour plot for angmcmc objects with bivariate data
#'
#' @inheritParams pointest
#' @param x angular MCMC object (with bivariate data).
#' @param show.data logical. Should the data points be added to the contour plot? Ignored if \code{object} is NOT supplied.
#' @param nlevels	number of contour levels desired \strong{if} levels is not supplied;
#' passed to the \link{contour} function in graphics.
#' @param levels numeric vector of levels at which to draw contour lines;
#' passed to the \link{contour} function in graphics.
#' @param cex,col graphical parameters passed to \code{\link{points}} in graphics for plotting the data points.
#' Ignored if {show.data == FALSE}.
#' @param ... additional arguments to be passed to the function \code{\link{contour}}.
#'
#' @details
#' \code{contour.angmcmc} is an S3 function for angmcmc objects that calls \code{\link{contour}} from graphics.
#'
#' To estimate the mixture density required to construct the contour plot, first the parameter vector \eqn{\eta} is estimated
#' by applying \code{fn} on the MCMC samples, yielding the (consistent) Bayes estimate \eqn{\hat{\eta}}. Then the mixture density
#' \eqn{f(x|\eta)} at any point \eqn{x} is (consistently) estimated by \eqn{f(x|\hat{\eta})}.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' # now create a contour plot with the default first 1/3 as burn-in and thin = 1
#' contour(fit.vmsin.20)
#'
#' @export

contour.angmcmc <-  function(x, fn = mean, show.data = TRUE, nlevels = 20, levels, burnin = 1/3,
                             thin = 1, cex = 1, col = "red", ...)
{
  if(class(x) != "angmcmc" || x$type != "bi") stop("\"x\" is not a bivariate angmcmc object")
  if(missing(burnin)) {
    if(is.null(x$fixed.label)) burnin <- 1/3
    else burnin <- x$burnin
  }
  if(missing(thin)) {
    if(is.null(x$fixed.label)) thin <- 1
    else thin <- x$thin
  }
  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  if(missing(levels)) {
    levels <- exp(seq(-20,2, length.out = nlevels))
  } else {
    levels <- NULL
  }

  if(x$ncomp == 1) {
    all.par <- pointest(x, fn = fn, burnin = burnin, thin=thin)
    par.mat <- as.matrix(all.par[-1])
    pi.mix  <- all.par[1]
  } else {
    all.par <- pointest(x, fn = fn, burnin = burnin, thin=thin)
    par.mat <- t(all.par[, -1])
    pi.mix  <- all.par[, 1]
  }

  all_input <- list("par.mat" = par.mat, "pi.mix" = pi.mix, "levels" = levels, ...)
  if(x$model == "wnorm2") all_input$omega.2pi.mat = x$omega.2pi

  colnames_data <- colnames(x$data)

  if(is.null(colnames_data)) {
    xlab <- ylab <- ""
  } else {
    xlab <- colnames_data[1]
    ylab <- colnames_data[2]
  }

  if(x$ncomp > 1) {
    main <- paste("contour plot for fitted", x$ncomp, "component", x$model, "mixtures")
  } else {
    main <- paste("contour plot for fitted (single component)", x$model)
  }

  do.call(paste0("contour_", x$model), all_input)
  if(show.data) points(x$data, col = col, cex = cex)

  graphics::title(main = main, xlab = xlab, ylab = ylab)
}


#' Density surface for angmcmc objects with bivariate data
#'
#' @inheritParams pointest
#' @param object angular MCMC object (with bivariate data).
#' @param log.density logical. Should log density be used for the plot?
#' @param theta,phi,shade,expand arguments passed to \code{\link{persp}} from \code{\link{graphics}}.
#' @param ... additional arguments passed to \code{\link{persp}} from \code{\link{graphics}}.
#'
#' @details
#' \code{densityplot2d} is a wrapper for \code{\link{persp}} from graphics applied on angmcmc objects.
#'
#' To estimate the mixture density, first the parameter vector \eqn{\eta} is estimated
#' by applying \code{fn} on the MCMC samples, yielding the (consistent) Bayes estimate \eqn{\hat{\eta}}. Then the mixture density
#' \eqn{f(x|\eta)} at any point \eqn{x} is (consistently) estimated by \eqn{f(x|\hat{\eta})}.
#'
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' # now create density surface with the default first 1/3 as burn-in and thin = 1
#' densityplot2d(fit.vmsin.20)
#' # the viewing angles can be changed through the arguments theta and phi
#' # (passed to persp from graphics)
#' densityplot2d(fit.vmsin.20, theta = 45, phi = 45)
#'
#' @export

densityplot2d <- function(object, fn = mean,  burnin = 1/3,
                          thin = 1, log.density = FALSE,
                          theta = 30, phi = 30, shade = 0.01, expand = 0.5, ...)
{
  if(class(object) != "angmcmc" || object$type != "bi") {
    stop("\"object\" is not a bivariate angmcmc object")
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


  est <- pointest(object, fn = fn, burnin = burnin, thin = thin)
  estlist <- list_by_row(t(est))

  x <- seq(0, 2*pi, length.out = 100)
  y <- seq(0, 2*pi, length.out = 100)
  coord <- as.matrix(expand.grid(x,y))

  inputpar <- addtolist(estlist, x = coord)
  if(object$model == "wnorm2") inputpar$int.displ <- object$int.displ

  if(object$ncomp  > 1) {
    den <- do.call(paste0("d", object$model, "mix"), args = inputpar)
  } else {
    den <- do.call(paste0("d", object$model), args = inputpar)
  }

  if(log.density) {
    denmat <- matrix(log(den), nrow = 100)
  } else {
    denmat <- matrix(den, nrow = 100)
  }
  nrden <- nrow(denmat)
  ncden <- ncol(denmat)

  if(object$ncomp > 1) {
    main <- paste("Density surface for fitted", object$ncomp, "component",
                  object$model, "mixtures")
  } else {
    main <- paste("Density surface for fitted (single component)", object$model)
  }

  # Create a function interpolating colors in the range of specified colors
  jet.colors <- grDevices::colorRampPalette( c("blue", "green",
                                               "yellow", "orange", "red") )
  # Generate the desired number of colors from this palette
  nbcol <- 500
  color <- jet.colors(nbcol)

  denfacet <- denmat[-1, -1] + denmat[-1, -ncden] +
    denmat[-nrden, -1] + denmat[-nrden, -ncden]
  # Recode facet z-values into color indices
  facetcol <- cut(denfacet, nbcol)

  graphics::persp(x, y, z=denmat, theta = theta, phi = phi, expand = expand, col = color[facetcol],
                  ltheta = 120, shade = shade, ticktype = "detailed",
                  xlab = "phi", ylab = "psi", zlab = "Density",
                  main = main) -> res
}



#' Density curve for angmcmc object with univariate data
#'
#' @inheritParams pointest
#' @param object angular MCMC object (with univariate data).
#' @param show.hist logical. Should a histogram for the data points be added to the plot?
#' @param ... other arguments to be passed to the function \code{\link{hist}} from graphics. Ignored if
#' \code{show.hist == FALSE}.
#'
#' @details
#' To estimate the mixture density, first the parameter vector \eqn{\eta} is estimated
#' by applying \code{fn} on the MCMC samples, yielding the (consistent) Bayes estimate \eqn{\hat{\eta}}. Then the mixture density
#' \eqn{f(x|\eta)} at any point \eqn{x} is (consistently) estimated by \eqn{f(x|\hat{\eta})}.
#'
#' @examples
#' # first fit a vm mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vm.20 <- fit_vmmix(wind, ncomp = 2, n.iter =  20,
#'                        ncores = 1)
#' # now create density curve with the default first 1/3 as burn-in and thin = 1
#' densityplot1d(fit.vm.20)
#'
#' @export


densityplot1d <- function(object, fn = mean, show.hist = TRUE, burnin = 1/3, thin = 1, ...)
{
  if(class(object) != "angmcmc" || object$type != "uni"){
    stop("\"object\" is not a univariate angmcmc object")
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

  if(object$ncomp == 1) {
    all.par <- pointest(object, fn = fn, burnin = burnin, thin=thin)
    par.mat <- as.matrix(all.par[-1])
    pi.mix  <- all.par[1]
  } else {
    all.par <- pointest(object, fn = fn, burnin = burnin, thin=thin)
    par.mat <- t(all.par[, -1])
    pi.mix  <- all.par[, 1]
  }

  den.arg <- seq(0, 2*pi, length.out = 100)
  l.const <- do.call(paste0("log_const_uni", object$model, "_all"), list("par" = par.mat))

  all_input <- list("x" = den.arg, "par" = par.mat, "pi" = pi.mix, "log_c" = l.const)
  if(object$model == "wnorm") all_input$omega_2pi_1d = object$omega.2pi.1d

  den <- do.call(paste0("uni", object$model, "mix_manyx"), all_input)

  if(show.hist){
    histplot <-  hist(object$data, plot = FALSE, ...)
  } else {
    histplot <- NULL
  }
  y_max <- 1.1* max(den, histplot$density)
  graphics::plot(den.arg,den, type = "l", lty = 2, ylim = c(0, y_max), ylab = "Density", xlab = "Angles in Radian")
  if(show.hist) plot(histplot, freq = FALSE, add = TRUE, ...)

  if(object$ncomp > 1) {
    main <- paste("density plot for fitted", object$ncomp, "component", object$model, "mixtures")
  } else {
    main <- paste("density plot for fitted (single component)", object$model)
  }

  graphics::title(main = main)
}







#' Trace plot for parameters from an angmcmc object
#' @inheritParams pointest
#' @param object angular MCMC object.
#' @param par parameter for which trace plot is to be created.
#'
#' @return
#' Returns a single plot if a single \code{par} and a single \code{comp.label} is supplied.
#' Otherwise, a series of plots is produced.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' # trace plot for kappa1 in component 1
#' paramtrace(fit.vmsin.20, "kappa1", 1)
#' # for kappa1 in all components
#' paramtrace(fit.vmsin.20, "kappa1")
#' # for all parameters in component 1
#' paramtrace(fit.vmsin.20, comp.label = 1)
#'
#' @export

paramtrace <- function(object, par, comp.label, burnin = 1/3, thin = 1)
{
  if(class(object) != "angmcmc") stop("paramtrace can only be used for \'angmcmc\' objects")
  if(missing(par)) par <- object$par.name
  if(missing(comp.label)) comp.label <- 1:object$ncomp
  if(any(c(length(par), length(comp.label)) > 1)) singleplot <- FALSE
  else singleplot <- TRUE
  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  if(singleplot) {
    val <- extractsamples(object, par, comp.label, burnin, thin)
    plot(val, type = "l", ylab="", xlab = "Iteration")
    if(object$ncomp > 1) {
      ylab <- paste(par, "for component", comp.label)
      main <- paste("Traceplot for", ylab, "in", object$ncomp, "component", object$model, "mixtures")
    } else {
      ylab <- par
      main <- paste("Traceplot for", ylab, "in (single component)", object$model)
    }

    graphics::title(main = main, ylab = ylab)
  } else {
    nplots <- length(par) * length(comp.label)
    currplotno <- 1L
    for(par.curr in par) {
      for(comp.label.curr in comp.label) {
        val <- extractsamples(object, par.curr, comp.label.curr, burnin, thin)
        plot(val, type = "l", ylab="", xlab = "Iteration")
        if(object$ncomp > 1) {
          ylab <- paste(par.curr, "for component", comp.label.curr)
          main <- paste("Traceplot for", ylab, "in ", object$ncomp, "component", object$model, "mixtures")
        } else {
          ylab <- par.curr
          main <- paste("Traceplot for", ylab, "in  (single component)", object$model)
        }

        graphics::title(main = main, ylab = ylab)
        if(currplotno < nplots) {
          press_enter()
          graphics::frame()
        }
        currplotno <- currplotno + 1
      }
    }
  }
}


#' Trace plot of log posterior density or log likelihood from an angmcmc object
#' @inheritParams pointest
#' @param object angular MCMC object.
#' @param use.llik logical. Should log likelihood be plotted instead of log posterior? Set
#' to \code{FALSE} by default.
#'
#' @examples
#' # first fit a vmsin mixture model
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' # log posterior density trace
#' lpdtrace(fit.vmsin.20)
#' # log likelihood trace
#' lpdtrace(fit.vmsin.20, use.llik = TRUE)
#'
#' @export

lpdtrace <- function(object, use.llik = FALSE, burnin = 1/3, thin = 1)
{
  if(class(object) != "angmcmc") stop("lpdtrace can only be used for \'angmcmc\' objects")
  if(burnin < 0 || burnin >= 1) stop("\"burnin\" must be in [0, 1)")
  if(thin <= 0  || thin != as.integer(thin) ) stop("\"thin\" must be a positive integer")

  burnin_iter <- 1:floor((object$n.iter)*burnin)
  thin_filter <- c(TRUE, rep(FALSE, thin-1))
  final_iter <- (1:object$n.iter)[-burnin_iter][thin_filter]

  if(use.llik){
    val <- object$llik[final_iter]
    plot(val, type = "l", ylab="Log Likelihood", xlab = "Iteration")
    if(object$ncomp > 1) {
      main <- paste("Log Likelihood traceplot for ", object$ncomp, "component", object$model, "mixtures")
    } else {
      main <- paste("Log Likelihood traceplot for (single component)", object$model)
    }
    graphics::title(main = main)
  } else {
    val <- object$lpd[final_iter]
    plot(val, type = "l", ylab="Log Posterior Density", xlab = "Iteration")
    if(object$ncomp > 1) {
      main <- paste("Log Posterior Density traceplot for", object$ncomp, "component", object$model, "mixtures")
    } else {
      main <- paste("Log Posterior Density traceplot fitted (single component)", object$model)
    }
    graphics::title(main = main)
  }

}
