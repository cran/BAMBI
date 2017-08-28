#' @useDynLib BAMBI, .registration = TRUE
#' @import stats
#' @importFrom graphics contour hist plot points persp
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#' @importFrom utils tail txtProgressBar
#' @importFrom parallel detectCores
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom label.switching dataBased ecr.iterative.1



NULL


#' @export
print.angmcmc <- function(x, ...) {

  output <-  paste("Dataset consists of", x$n.data, "observations.")


  if(grepl(x$method, "hmc")) {
    output[2] <- paste(x$ncomp, "cluster", x$model, "mixture fitted via HMC for model parameters.")

    if(x$epsilon.random) {
      output[3] <- paste("epsilon chosen randomly at each iteration with average epsilon =",
                         format(x$epsilon, scientific=TRUE, digits = 2))
    } else {
        output[3] <- paste("epsilon fixed at", x$epsilon)
    }
    if(x$L.random){
      output[4] <- paste("L chosen at each iteration with average L =", x$L)
    } else {
      output[4] <- paste("L fixed at", x$L)
    }
    output[5] <- paste("acceptance rate for model parameters = ",
                       round(100*mean(x$accpt.modelpar), 2), "%.")
  }


   else if(grepl(x$method, "rwmh")) {

    output[2] <- paste(x$ncomp, "cluster", x$model, "mixture fitted via RWMH for model parameters.")

    output[3] <- paste("proposals are independent normal with variances",
                       paste(format(x$propscale.final, scientific=TRUE, digits = 2), sep = "", collapse = ", "),
                       "for", paste(x$par.name[-1], sep = "", collapse=", "))

    output[4] <- paste("acceptance rate for concentration parameters = ",
                       round(100*mean(x$accpt.kappa), 2), "%.")

    output[5] <- paste("acceptance rate for mean parameters = ",
                       round(100*mean(x$accpt.mu), 2), "%.")
  }

  output[6] <- paste("Number of iterations =", x$n.iter)
  cat(output, sep = "\n")
}


.onUnload <- function (libpath) {
  library.dynam.unload("BAMBI", libpath)
}

#' @export
print.stepfit <- function(x, ...)
{
  if(x$check_min) {
    output <- paste("First minimum attained at ncomp =", x$ncomp.best)
    output[2] <- paste("Extract the best fit with the keyword \'fit.best\'")
    cat("\n")
    cat(output, sep = "\n")
  } else {
    warning(paste(toupper(x$crit), "did not attend a first minimum. Probably more clusters are needed."))
    cat("\n")
    cat("Extract the last fit with the keyword \'fit.best\'")
  }
}


#' Angular MCMC (\code{angmcmc}) Object
#' @exportClass angmcmc
#' @description Checking if an R object is angmcmc
#' @param object any R object
#' @return logical. Is the input an angmcmc object?
#' @details
#' \code{angmcmc} objects are classified lists that are created when any of the five mixture model fitting
#' functions, viz., \code{fit_vmmix}, \code{fit_wnormmix}, \code{fit_vmsinmix}, \code{fit_vmcosmix} and
#' \code{fit_wnorm2mix} is used. An \code{angmcmc} object contains a number of elements, including the dataset, the
#' model being fitted on the dataset and dimension of the model (univariate or bivariate), the tuning parameters
#' used, MCMC samples for the mixture model parameters, the (hidden) component or cluster indicators for  data
#' points in each iteration and the (iteration-wise) log likelihood and log posterior density values (both calculated
#' upto some normalizing constants). When printed, an angmcmc object returns a brief summary of the function
#' arguments used to produce the object and the average acceptance rate of the proposals (in HMC and RWMH) used
#' over iterations. An \code{angmcmc} object can be used as an argument for the diagnostic and post-processing
#' functions available in \code{BAMBI} for making further inferences.
#' @examples
#' # illustration only - more iterations needed for convergence
#' fit.vmsin.20 <- fit_vmsinmix(tim8, ncomp = 3, n.iter =  20,
#'                              ncores = 1)
#' is.angmcmc(fit.vmsin.20)
#' @export

is.angmcmc <- function(object) {
  methods::is(object, "angmcmc")
}


