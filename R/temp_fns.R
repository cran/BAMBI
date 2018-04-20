fit_model_2 <- function(data, model = "vmsin", ncomp, start_par = list(), method="hmc", epsilon=0.01, L=10, epsilon.random=TRUE,
                        L.random=FALSE, propscale = rep(0.01, 5), n.iter=500, gam.loc=0, gam.scale=1000, pmix.alpha = 1/2,
                        burnin = 1/3, thin = 1,
                        norm.var=1000, autotune = FALSE,
                        iter.tune = 10, ncores,
                        show.progress = TRUE,
                        accpt.prob.upper = 0.8,
                        accpt.prob.lower = 0.7, tune.incr = 0.05) {

  burnin_iter <- seq_len(floor((n.iter)*burnin))
  burnin_iter_max <- max(burnin_iter)

  if(is.null(dim(data)) | !(mode(data) %in% c("list", "numeric") && ncol(data) == 2)) stop("non-compatible data")

  epsilon.start <- epsilon.orig <- epsilon.final <- rep(epsilon, ncomp)
  L.orig <- L

  kappa_upper <- 150

  curr.model <- model
  data.rad <- rm_NA_rad(data)
  n.data <- nrow(data.rad)

  if(missing(ncores)) {
    ncores <- floor(parallel::detectCores())
  }

  if(ncomp == 1) {

    if(missing(start_par)) {
      starting <- start_clus_kmeans_vmsin(data.rad, ncomp, nstart=5)
      starting$par.mat <- matrix(starting$par.mat, ncol=1)
    } else {
      allpar <- start_par
      if(any(is.null(allpar$kappa1), is.null(allpar$kappa2), is.null(allpar$kappa3),
             is.null(allpar$mu1), is.null(allpar$mu2)) ) {
        stop("too few elements in start_par, with no default")
      }
      allpar1 <- list(allpar$kappa1, allpar$kappa2, allpar$kappa3, allpar$mu1, allpar$mu2)
      allpar_len <- listLen(allpar1)
      if(min(allpar_len) != max(allpar_len)){
        stop("component size mismatch: number of components of the input parameter vectors differ")
      }
      starting <- list("par.mat" = rbind(start_par$kappa1, start_par$kappa2, start_par$kappa3,
                                         start_par$mu1, start_par$mu2), "pi.mix" = 1)
    }
  } else if(ncomp > 1) {
    if(missing(start_par)) {
      starting <- start_clus_kmeans_vmsin(data.rad, ncomp, nstart=5)
    } else {
      allpar <- start_par
      if(any(is.null(allpar$kappa1), is.null(allpar$kappa2), is.null(allpar$kappa3),
             is.null(allpar$mu1), is.null(allpar$mu2), is.null(allpar$pmix)) ) {
        stop("too few elements in start_par, with no default")
      }
      allpar1 <- list(allpar$kappa1, allpar$kappa2, allpar$kappa3, allpar$mu1, allpar$mu2, allpar$pmix)
      allpar_len <- listLen(allpar1)
      allpar_len <- listLen(allpar)
      if(min(allpar_len) != max(allpar_len)){
        stop("component size mismatch: number of components of the input parameter vectors differ")
      }
      starting <- list("par.mat" = rbind(start_par$kappa1, start_par$kappa2, start_par$kappa3,
                                         start_par$mu1, start_par$mu2), "pi.mix" = start_par$pmix)
    }
  }

  starting$par.mat[abs(starting$par.mat) >= kappa_upper/2] <- kappa_upper/2
  starting$l.c.vmsin <- as.numeric(log_const_vmsin_all(starting$par.mat))
  starting$llik <- vapply(1:ncomp,
                          function(j)
                            ifelse(starting$pi.mix[j] == 0, -Inf,
                                   do.call(paste0("llik_", model, "_one_comp"),
                                           list(data = data.rad[starting$clust.ind[[j]], ],
                                                par_vec = starting$par.mat[, comp],
                                                log_c = starting$l.c.vmsin[comp],
                                                ncores = ncores))), 0)
  starting$llik.full <- do.call(paste0("llik_", model, "_full"),
                                list(data = data.rad, par = starting$par.mat,
                                     pi = starting$pi.mix, log_c = starting$l.c.vmsin,
                                     ncores = ncores))

  starting$lprior <- vapply(1:ncomp, function(j)
    sum(ldgamanum(starting$par.mat[1:2, j], gam.loc, gam.scale))
    - 0.5*sum((starting$par.mat[3, j]/norm.var)^2)
    , 0)

  starting$lprior.full <- sum((pmix.alpha-1)*log(starting$pi.mix)) + sum(ldgamanum(starting$par.mat[1:2,], gam.loc, gam.scale)) - 0.5*sum((starting$par.mat[3,]/norm.var)^2)

  starting$lpd <- starting$llik + starting$lprior
  starting$lpd.full <- starting$llik.full + starting$lprior.full

  par.mat.all <- array(0, dim = c(5, ncomp, n.iter+1))
  pi.mix.all <- matrix(1, nrow = ncomp, ncol = n.iter+1)
  llik.full.all <- lprior.full.all <- lpd.full.all <- rep(-Inf, (n.iter+1))
  llik.all <- lprior.all <- lpd.all <- matrix(-Inf, (n.iter+1), ncomp)
  accpt.par.mat.all <- accpt.kappa.all <- accpt.mu.all <- matrix(0, (n.iter+1), ncomp)
  modelpar.names <- c("kappa1", "kappa2", "kappa3", "mu1", "mu2")

  MC <- starting  #simulation results list, 1st entry = method of moments on kmeans output

  par.mat.all[,,1] <- MC$par.mat
  pi.mix.all[,1] <- MC$pi.mix

  llik.all[1, ] <- MC$llik
  lprior.all[1, ] <- MC$lprior
  lpd.all[1, ] <- MC$lpd
  llik.full.all[1 ] <- MC$llik.full
  lprior.full.all[1 ] <- MC$lprior.full
  lpd.full.all[1 ] <- MC$lpd.full

  epsilon_ave <- NULL
  L_ave <- NULL
  propscale_final <- NULL

  clus.ind <- matrix(1, nrow = n.data, ncol = n.iter+2)

  iter <- 2
  ntune <- 0

  if(show.progress) pb <- txtProgressBar(min = 2, max = n.iter+1, style = 3)

  #******************************************************************************************
  # single component model
  #******************************************************************************************

  if(ncomp == 1 && grepl(method, "hmc")) # using hmc
  {
    if(epsilon.random) {
      epsilon_vec <- runif(n.iter, min = 0.9*epsilon, max = 1.1*epsilon)
      # epsilon.0 <- epsilon
    }
    if(L.random) {
      L_vec <- sample(1:L, n.iter, replace = TRUE)
      # L.0 <- L
    }

    while(iter <= (n.iter+1)) {
      broken <- FALSE
      kappa.large <- FALSE

      pi.mix.1 <- 1
      par.mat.old <- MC$par.mat
      l.c.vmsin.old <- MC$l.c.vmsin

      llik_new.pi <- MC$llik

      #----------------------------------------------------------------------------------
      #generating par.mat by HMC
      #----------------------------------------------------------------------------------

      par.mat.1 <- par.mat.old
      lprior.1 <- MC$lprior
      llik.1 <- llik_new.pi
      lpd.1 <- llik.1 + lprior.1
      l.c.vmsin.1 <- MC$l.c.vmsin
      accpt.par.mat <- 0


      current_q <- par.mat.1
      current_p <- matrix(rnorm(5*ncomp,0,1), nrow = 5)  # independent standard normal variates
      p <- current_p
      q <- current_q

      if(L.random) {
        L <- L_vec[iter-1]
        # L <- sample.int(1:L.0, n.iter, replace = TRUE)
      }
      if(epsilon.random) {
        epsilon <- epsilon_vec[iter-1]
        # epsilon <- runif(n.iter, min = 0.9*epsilon.0, max = 1.1*epsilon.0)
      }


      # Do leapfrog with L and epsilon
      {
        # Make a half step for momentum at the beginning

        p <- p - (epsilon/2) * (- grad_vmsin_all_comp(data.rad, q, pi.mix.1, ncores)
                                + matrix(c(1/gam.scale + (1- 1/gam.scale)/q[1:2,], q[3,], rep(0,2)), ncol=1)
        ) # the second term in the bracket arises from prior
        # Alternate full steps for position and momentum

        for (i in 1:L)
        {
          # Make a full step for the position

          q <- q + epsilon * p

          if(all(!is.nan(q)) && any(abs(q[1:3, ]) >= kappa_upper)) {
            kappa.large <- TRUE
            break
          }

          if(any(is.nan(c(q,p))))  {
            broken <- TRUE
            #stop("Algorithm breaks. Try a smaller epsilon.")
            break
          }

          # Make sure the components of q1 are in the proper ranges
          {
            q1 <- q; p1 <- p


            for(j in 1:ncomp) {

              while(q1[1,j] <= 0) {
                q1[1,j] <- -q1[1,j]; p1[1,j] <- -p1[1,j]
              }

              while(q1[2,j] <= 0) {
                q1[2,j] <- -q1[2,j]; p1[2,j] <- -p1[2,j]
              }

              # q3 is unbounded in vmsin

              while(q1[4,j] < 0 || q1[4,j] >= 2*pi) {
                if(q1[4,j] < 0) {
                  q1[4,j] <- - q1[4,j]; p1[4,j] <- -p1[4,j]
                } else {
                  q1[4,j] <- 4*pi - q1[4,j]; p1[4,j] <- -p1[4,j]
                }
              }
              while(q1[5,j] < 0 || q1[5,j] >= 2*pi) {
                if(q1[5,j] < 0) {
                  q1[5,j] <- - q1[5,j]; p1[5,j] <- -p1[5,j]
                } else {
                  q1[5,j] <- 4*pi - q1[5,j]; p1[5,j] <- -p1[5,j]
                }
              }
            }


            p <- p1; q <- q1
          }
          # Make a full step for the momentum, except at end of trajectory

          if(any(is.nan(c(p, q))))  {
            broken <- TRUE
            #stop("Algorithm breaks. Try a smaller epsilon.")
            break
          } else if(i!=L) {
            p <- p - epsilon * (- grad_vmsin_all_comp(data.rad, q, pi.mix.1, ncores)
                                + matrix(c(1/gam.scale + (1- 1/gam.scale)/q[1:2,], q[3,], rep(0,2)), ncol=1 )) # the second term in the bracket arises from prior
          }
        }

        # Make a half step for momentum at the end.
        if(all(!broken, !kappa.large)){
          if(any(is.nan(c(p, q))))  {
            broken <- TRUE
          } else {
            p <- p - (epsilon/2) * (- grad_vmsin_all_comp(data.rad, q, pi.mix.1, ncores)
                                    + matrix(c(1/gam.scale + (1- 1/gam.scale)/q[1:2,], q[3,], rep(0,2)), ncol=1 ) ) # the second term in the bracket arises from prior

          }
        }
        # Negate momentum at end of trajectory to make the proposal symmetric

        if(any(is.nan(c(p, q))))  {
          broken <- TRUE
        } else {
          p <-  -p
        }
      }


      if (iter > 100 && mean(accpt.par.mat.all[1:iter]) < 0.05) {
        broken <- TRUE
      }

      if (broken) {
        print("Acceptance rate too low. Automatically restarting with a smaller \'epsilon\'.")
        iter <- 2
        if(epsilon.random) {
          epsilon_vec <- epsilon_vec/2
        } else {
          epsilon <- epsilon/2
        }


        par.mat.all[,,iter] <- par.mat.1
        pi.mix.all[,iter] <- pi.mix.1
        llik.all[iter] <- llik.1
        lprior.all[iter] <- lprior.1
        lpd.all[iter] <- lpd.1
        accpt.par.mat.all[iter] <- accpt.par.mat


        next
      }



      # Evaluate potential and kinetic energies at start and end of trajectory
      current_U <- -lpd.1
      current_K <- sum(current_p^2) / 2

      if(kappa.large) {
        proposed_U <- proposed_K <- Inf
      } else {

        par.mat.prop <- q
        l.c.vmsin.prop <- log_const_vmsin_all(par.mat.prop)

        lprior.prop <- sum((pmix.alpha-1) * log(pi.mix.1)) + sum(ldgamanum(q[1:2,], gam.loc, gam.scale)) - 0.5*sum((q[3,]/norm.var)^2)

        llik.prop <- llik_vmsin_full(data.rad, q, pi.mix.1, l.c.vmsin.prop, ncores)

        proposed_U <- -(llik.prop + lprior.prop)
        proposed_K <- sum(p^2) / 2
      }

      exp(current_U-proposed_U+current_K-proposed_K)
      # Accept or reject the state at end of trajectory, returning either
      # the position at the end of the trajectory or the initial position

      if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
      {
        # return (q)  # accept
        # accpt = 1
        par.mat.1 <- signif(par.mat.prop, 8)
        lprior.1 <- signif(lprior.prop, 8)
        llik.1 <- signif(llik.prop, 8)
        lpd.1 <- signif(-proposed_U, 8)
        accpt.par.mat <- 1
        l.c.vmsin.1 <- signif(l.c.vmsin.prop, 8)
      }


      MC <- list("par.mat" = par.mat.1, "pi.mix" = pi.mix.1,
                 "l.c.vmsin" = l.c.vmsin.1, "llik" = llik.1, "lprior" = lprior.1, "lpd" = lpd.1,
                 "accpt.par.mat" = accpt.par.mat)

      par.mat.all[,,iter] <- MC$par.mat
      pi.mix.all[,iter] <- MC$pi.mix
      llik.all[iter] <- llik.1
      lprior.all[iter] <- lprior.1
      lpd.all[iter] <- lpd.1
      accpt.par.mat.all[iter] <- accpt.par.mat


      # tuning epsilon with first 20 draws
      if(autotune && iter == iter.tune && mean(accpt.par.mat.all[2:(iter.tune+1)]) < 0.6) {
        iter <- 2
        ntune <- ntune + 1
        if(epsilon.random) {
          epsilon_vec <- epsilon_vec/2
        } else {
          epsilon <- epsilon/2
        }
      }

      if(show.progress && ((iter-1) %% 25 == 0 || iter == n.iter + 1))
        utils::setTxtProgressBar(pb, iter)

      iter <- iter + 1

    }
  }

  if(ncomp == 1 && grepl(method, "rwmh")) # using rwmh
  {
    while(iter <= (n.iter+1)) {
      pi.mix.1 <- 1
      par.mat.old <- MC$par.mat
      l.c.vmsin.old <- MC$l.c.vmsin

      llik_new.pi <- MC$llik
      #----------------------------------------------------------------------------------
      #generating presicion parameters
      #----------------------------------------------------------------------------------

      k1.1.old <- par.mat.old[1, ]
      k2.1.old <- par.mat.old[2, ]
      k3.1.old <- par.mat.old[3, ]

      k1.1.prop <- pmax(k1.1.old + rnorm(ncomp,0,propscale[1]), 1e-6)
      k2.1.prop <- pmax(k2.1.old + rnorm(ncomp,0,propscale[2]), 1e-6)
      k3.1.prop <- k3.1.old + rnorm(ncomp,0,propscale[3])

      prop.mat <- unname(matrix(c(k1.1.prop,k2.1.prop,k3.1.prop, par.mat.old[4:5, ]),ncol=1))
      l.c.vmsin.prop <- as.numeric(log_const_vmsin_all(prop.mat))

      llik_old <- llik_new.pi
      llik_prop <- llik_vmsin_full(data.rad, prop.mat, pi.mix.1, l.c.vmsin.prop, ncores)

      lprior_old <- MC$lprior
      lprior_prop <- sum((pmix.alpha-1) * log(pi.mix.1)) + sum(ldgamanum(prop.mat[1:2,], gam.loc, gam.scale)) - 0.5*sum((prop.mat[3,]/norm.var)^2)

      post.omg_old <- llik_old + lprior_old
      post.omg_prop <- llik_prop + lprior_prop

      if (runif(1) <  exp(post.omg_prop-post.omg_old) ) {
        k1.1 <- k1.1.prop
        k2.1 <- k2.1.prop
        k3.1 <- k3.1.prop
        accpt.kappa <- 1
        l.c.vmsin.1 <- l.c.vmsin.prop
        llik_new.omg <- llik_prop
        lprior.1 <- lprior_prop
        par.mat_new.omg <- prop.mat
      } else {
        k1.1 <- k1.1.old
        k2.1 <- k2.1.old
        k3.1 <- k3.1.old
        accpt.kappa <- 0
        l.c.vmsin.1 <- l.c.vmsin.old
        llik_new.omg <- llik_old
        lprior.1 <- lprior_old
        par.mat_new.omg <- par.mat.old
      }


      #----------------------------------------------------------------------------------
      #generating mu and nu
      #----------------------------------------------------------------------------------
      prop.mu <- prncp_reg(par.mat.old[4, ] + rnorm(ncomp,0,propscale[4]))
      prop.nu <- prncp_reg(par.mat.old[5, ] + rnorm(ncomp,0,propscale[5]))
      #----------------------------------------------------------------------------------
      prop.mat.mean <- matrix(c(par.mat_new.omg[1:3,], prop.mu,prop.nu), ncol=1)

      llik_new.prop <- llik_vmsin_full(data.rad, prop.mat.mean, pi.mix.1, l.c.vmsin.1, ncores)

      if (runif(1) <  exp(llik_new.prop-llik_new.omg) ) {
        par.mat.1 <- prop.mat.mean
        accpt.mu <- 1
        llik.1 <- llik_new.prop
      } else {
        par.mat.1 <- par.mat_new.omg
        accpt.mu <- 0
        llik.1 <- llik_new.omg
      }

      lpd.1 <- llik.1 + lprior.1

      MC <- list("par.mat" = par.mat.1, "pi.mix" = pi.mix.1,
                 "l.c.vmsin" = l.c.vmsin.1, "llik" = llik.1, "lprior" = lprior.1, "lpd" = lpd.1,
                 "accpt.kappa" = accpt.kappa, "accpt.mu" = accpt.mu)

      par.mat.all[,,iter] <- MC$par.mat
      pi.mix.all[,iter] <- MC$pi.mix
      llik.all[iter] <- llik.1
      lprior.all[iter] <- lprior.1
      lpd.all[iter] <- lpd.1
      accpt.kappa.all[iter] <- accpt.kappa
      accpt.mu.all[iter] <- accpt.mu

      # tuning propscale with first 20 draws
      if(autotune && iter == iter.tune && (mean(accpt.kappa.all[2:(iter.tune+1)]) < 0.6 ||
                                           mean(accpt.mu.all[2:(iter.tune+1)]) < 0.6)) {
        iter <- 2
        ntune <- ntune + 1
        propscale <- propscale/2
      }

      if(show.progress && ((iter-1) %% 25 == 0 || iter == n.iter + 1))
        utils::setTxtProgressBar(pb, iter)

      iter <- iter+1

    }
  }
  #******************************************************************************************

  #******************************************************************************************
  # multiple component model
  #******************************************************************************************

  if(ncomp > 1 && grepl(method, "hmc")) # using hmc
  {
    epsilon_vec <- matrix(epsilon, n.iter, ncomp)
    L_vec <- rep(L, n.iter)
    eps.tune <- matrix(NA, n.iter, ncomp)

    while(iter <= (n.iter+1)) {
      #----------------------------------------------------------------------------------
      #generating mixture proportions
      #----------------------------------------------------------------------------------
      pi.mix.old <- MC$pi.mix
      par.mat.old <- MC$par.mat
      l.c.vmsin.old <- MC$l.c.vmsin


      # Gibbs Sampler
      {
        post.wt <- mem_p_sin(data.rad, par.mat.old, pi.mix.old, l.c.vmsin.old, ncores)
        cid.temp <- cID(post.wt, ncomp, runif(n.data))

        # random label switching
        # id.perm <- sample(1:ncomp)
        # cid.curr <- clus.ind[ , iter] <- change_labs(cid.temp, id.perm)
        cid.curr <- clus.ind[ , iter] <- cid.temp

        n.clus <- tabulate(clus.ind[ , iter], nbins = ncomp)
        pi.mix.1 <- as.numeric(rdirichlet(1, (pmix.alpha + n.clus))) #new mixture proportions
        llik_new.pi <- llik_vmsin_full(data.rad, par.mat.old, pi.mix.1, l.c.vmsin.old, ncores)
      }



      #----------------------------------------------------------------------------------
      #generating par.mat by HMC, indep for each component
      #----------------------------------------------------------------------------------
      par.mat.1 <- par.mat.old
      lprior.1 <- MC$lprior
      lpd.1 <- MC$lpd
      llik.1 <- MC$llik


      if(epsilon.random) {
        epsilon_vec[iter-1, ] <- epsilon.final <-
          runif(ncomp, epsilon.orig *.99, epsilon.orig * 1.01)
      } else {
        epsilon.final <- epsilon.orig
      }

      if(L.random) {
        L_vec[iter-1] <- L <- sample(1:L.orig, 1)
      } else {
        L <- L.orig
      }


      for(comp in 1:ncomp) {

        # do hmc for current comp
        clus_obs <- which(cid.curr == comp)
        n_clus_obs <- length(clus_obs)

        lpd_grad_model_1comp_noarg <- function(par, grad) {
          llik.comp <- ifelse(n_clus_obs == 0, -Inf,
                              do.call(paste0("llik_", model, "_one_comp"),
                                      list(data = data.rad[clus_obs, ],
                                           par_vec = par,
                                           log_c = do.call(paste0("log_const_", model,"_all"),
                                                           list(par_mat = as.matrix(par))),
                                           ncores = ncores)))

          if (grad) {
            grad.comp <- ifelse(n_clus_obs == 0, 0,
                                as.numeric(
                                  do.call(paste0("grad_", model, "_one_comp"),
                                          list(data = data.rad[clus_obs, ],
                                               par_vec = par, ncores = ncores))))
            attr(llik.comp, "grad") <- grad.comp
          }
          llik.comp
        }

        if (model == "wnorm2") {
          attr(lpd_grad_model_1comp_noarg, "lower") <- c(0, 0, -Inf, 0, 0)
          attr(lpd_grad_model_1comp_noarg, "upper") <- c(Inf, Inf, Inf, 2*pi, 2*pi)
        } else {
          attr(lpd_grad_model_1comp_noarg, "lower") <- c(0, 0, -Inf, 0, 0)
          attr(lpd_grad_model_1comp_noarg, "upper") <- c(Inf, Inf, Inf, 2*pi, 2*pi)
        }

        hmc.vmsin.comp <- bounded_hmc(
          lpd_grad_model_1comp_noarg,
          step = epsilon.final[comp],
          nsteps = L,
          initial = par.mat.old[, comp],
          return.traj = TRUE)


        par.mat.1[, comp] <- q <- hmc.vmsin.comp$final
        lpd.1[comp] <- hmc.vmsin.comp$lpr[[1]]
        lprior.1[comp] <- sum(ldgamanum(q[1:2],
                                        gam.loc, gam.scale)) -
          0.5*(q[3]/norm.var)^2
        llik.1[comp] <- lpd.1[comp] - lprior.1[comp]
        accpt.par.mat.all[iter, comp] <- hmc.vmsin.comp$acc

        # autotune epsilon in burnin
        if (iter >= 100 & iter %% 5 == 0 & iter <= burnin_iter_max) {
          # mean of acceptance in last 100 for each component
          # if (iter >= 500) browser()
          nterms <- 25
          mean.accpts <- sum(accpt.par.mat.all[(iter-nterms + 1):iter, comp])/nterms
          if (mean.accpts > accpt.prob.upper) {
            epsilon.orig[comp] <- epsilon.orig[comp] * (1 + tune.incr)
            eps.tune[iter-1, comp] <- "incr"
          } else if (mean.accpts < accpt.prob.lower) {
            epsilon.orig[comp] <- epsilon.orig[comp] * (1 - tune.incr)
            eps.tune[iter-1, comp] <- "decr"
          }
        }

      }




      l.c.vmsin.1 <- do.call(paste0("log_const_", model,"_all"),
                             list(par_mat = par.mat.1))
      llik.full.1 <- do.call(paste0("llik_", model, "_full"),
                             list(data = data.rad, par = par.mat.1,
                                  pi = pi.mix.1,
                                  log_c = l.c.vmsin.1,
                                  ncores = ncores))

      lprior.full.1 <- sum((pmix.alpha-1)*log(pi.mix.1)) +
        sum(ldgamanum(par.mat.1[1:2,], gam.loc, gam.scale)) -
        0.5*sum((par.mat.1[3,]/norm.var)^2)

      lpd.full.1 <- llik.full.1 + lprior.full.1



      MC <- list("par.mat" = par.mat.1, "pi.mix" = pi.mix.1,
                 "l.c.vmsin" = l.c.vmsin.1, "llik" = llik.1, "lprior" = lprior.1, "lpd" = lpd.1,
                 "llik.full" = llik.full.1, "lprior.full" = lprior.full.1, "lpd.full" = lpd.full.1,
                 "accpt.par.mat" = accpt.par.mat)

      par.mat.all[,,iter] <- par.mat.1
      pi.mix.all[,iter] <- pi.mix.1
      llik.all[iter, ] <- llik.1
      lprior.all[iter, ] <- lprior.1
      lpd.all[iter, ] <- lpd.1
      llik.full.all[iter] <- llik.full.1
      lprior.full.all[iter] <- lprior.full.1
      lpd.full.all[iter] <- lpd.full.1



      if(show.progress && ((iter-1) %% 25 == 0 || iter == n.iter + 1))
        utils::setTxtProgressBar(pb, iter)

      iter <- iter + 1

    }
  }

  if(ncomp > 1 && grepl(method, "rwmh")) # using rwmh
  {
    while(iter <= (n.iter+1)) {
      #----------------------------------------------------------------------------------
      #generating mixture proportions
      #----------------------------------------------------------------------------------
      pi.mix.old <- MC$pi.mix
      par.mat.old <- MC$par.mat
      l.c.vmsin.old <- MC$l.c.vmsin


      # Gibbs Sampler
      {
        post.wt <- mem_p_sin(data.rad, par.mat.old, pi.mix.old, l.c.vmsin.old, ncores)
        clus.ind[ , iter] <- cID(post.wt, ncomp, runif(n.data))
        n.clus <- tabulate(clus.ind[ , iter], nbins = ncomp)
        pi.mix.1 <- as.numeric(rdirichlet(1, (pmix.alpha + n.clus))) #new mixture proportions
        llik_new.pi <- llik_vmsin_full(data.rad, par.mat.old, pi.mix.1, l.c.vmsin.old, ncores)
      }


      #----------------------------------------------------------------------------------
      #generating presicion parameters
      #----------------------------------------------------------------------------------

      k1.1.old <- par.mat.old[1, ]
      k2.1.old <- par.mat.old[2, ]
      k3.1.old <- par.mat.old[3, ]

      k1.1.prop <- pmax(k1.1.old + rnorm(ncomp,0,propscale[1]), 1e-6)
      k2.1.prop <- pmax(k2.1.old + rnorm(ncomp,0,propscale[2]), 1e-6)
      k3.1.prop <- k3.1.old + rnorm(ncomp,0,propscale[3])

      prop.mat <- unname(rbind(k1.1.prop,k2.1.prop,k3.1.prop, par.mat.old[4:5, ]))
      l.c.vmsin.prop <- as.numeric(log_const_vmsin_all(prop.mat))

      llik_old <- llik_new.pi
      llik_prop <- llik_vmsin_full(data.rad, prop.mat, pi.mix.1, l.c.vmsin.prop, ncores)

      lprior_old <- MC$lprior
      lprior_prop <- sum((pmix.alpha-1) * log(pi.mix.1)) + sum(ldgamanum(prop.mat[1:2,], gam.loc, gam.scale)) - 0.5*sum((prop.mat[3,]/norm.var)^2)

      post.omg_old <- llik_old + lprior_old
      post.omg_prop <- llik_prop + lprior_prop

      if (runif(1) <  exp(post.omg_prop-post.omg_old) ) {
        k1.1 <- k1.1.prop
        k2.1 <- k2.1.prop
        k3.1 <- k3.1.prop
        accpt.kappa <- 1
        l.c.vmsin.1 <- signif(l.c.vmsin.prop, 8)
        llik_new.omg <- signif(llik_prop, 8)
        lprior.1 <- signif(lprior_prop, 8)
        par.mat_new.omg <- signif(prop.mat, 8)
      } else {
        k1.1 <- k1.1.old
        k2.1 <- k2.1.old
        k3.1 <- k3.1.old
        accpt.kappa <- 0
        l.c.vmsin.1 <- l.c.vmsin.old
        llik_new.omg <- llik_old
        lprior.1 <- lprior_old
        par.mat_new.omg <- par.mat.old
      }


      #----------------------------------------------------------------------------------
      #generating mu and nu
      #----------------------------------------------------------------------------------
      prop.mu <- prncp_reg(par.mat.old[4, ] + rnorm(ncomp,0,propscale[4]))
      prop.nu <- prncp_reg(par.mat.old[5, ] + rnorm(ncomp,0,propscale[5]))
      #----------------------------------------------------------------------------------
      prop.mat.mean <- unname(rbind(par.mat_new.omg[1:3,], prop.mu,prop.nu))

      llik_new.prop <- llik_vmsin_full(data.rad, prop.mat.mean, pi.mix.1, l.c.vmsin.1, ncores)

      if (runif(1) <  exp(llik_new.prop-llik_new.omg) ) {
        par.mat.1 <- signif(prop.mat.mean, 8)
        accpt.mu <- 1
        llik.1 <- llik_new.prop
      } else {
        par.mat.1 <- par.mat_new.omg
        accpt.mu <- 0
        llik.1 <- llik_new.omg
      }

      lpd.1 <- llik.1 + lprior.1

      MC <- list("par.mat" = par.mat.1, "pi.mix" = pi.mix.1,
                 "l.c.vmsin" = l.c.vmsin.1, "llik" = llik.1, "lprior" = lprior.1, "lpd" = lpd.1,
                 "accpt.kappa" = accpt.kappa, "accpt.mu" = accpt.mu)

      par.mat.all[,,iter] <- MC$par.mat
      pi.mix.all[,iter] <- MC$pi.mix
      llik.all[iter] <- llik.1
      lprior.all[iter] <- lprior.1
      lpd.all[iter] <- lpd.1
      accpt.kappa.all[iter] <- accpt.kappa
      accpt.mu.all[iter] <- accpt.mu

      # tuning propscale with first 20 draws
      if(autotune && iter == iter.tune && (mean(accpt.kappa.all[2:(iter.tune+1)]) < 0.6 ||
                                           mean(accpt.mu.all[2:(iter.tune+1)]) < 0.6)) {
        iter <- 2
        ntune <- ntune + 1
        propscale <- propscale/2
      }

      if(show.progress && ((iter-1) %% 25 == 0 || iter == n.iter + 1))
        utils::setTxtProgressBar(pb, iter)

      iter <- iter+1

    }
  }
  #******************************************************************************************
  if(grepl(method, "hmc")) {
    if(epsilon.random) {
      epsilon_ave <- colMeans(epsilon_vec)
    } else{
      epsilon_ave <- epsilon.final
    }
    # delete this line later after adjusting printing statements
    epsilon_ave <- mean(epsilon_ave)
    if(L.random) {
      L_ave <- mean(L_vec)
    } else{
      L_ave <- L
    }
  }

  # browser()

  if(grepl(method, "rwmh")) {
    propscale_final <- propscale
  }

  if(show.progress) cat("\n")

  allpar_val <- array(1, dim = c(6, ncomp, n.iter+1))
  allpar_val[1, , ] <- pi.mix.all
  allpar_val[2:6, , ] <- par.mat.all

  allpar_name <- c("pmix", modelpar.names)
  dimnames(allpar_val)[[1]] <- c("pmix", modelpar.names)

  result <- list("par.value" = allpar_val, "par.name" = allpar_name, "llik" = llik.all,
                 "accpt.modelpar" = accpt.par.mat.all,
                 "accpt.kappa" = accpt.kappa.all, "accpt.mu" = accpt.mu.all,
                 "lpd" = lpd.all, "model" = curr.model, "method" = method, "clus.ind" = clus.ind,
                 "epsilon.random" = epsilon.random, "epsilon" = epsilon_ave,
                 "L.random" = L.random, "L" = L_ave,  "type" = "bi",
                 "propscale.final" = propscale_final, "data" = data.rad,
                 "gam.loc" = gam.loc, "gam.scale" = gam.scale, "pmix.alpha" = pmix.alpha, "norm.var" = norm.var,
                 "n.data" = n.data, "ncomp" = ncomp, "n.iter" = n.iter)
  class(result) <- "angmcmc"

  return(result)
}


