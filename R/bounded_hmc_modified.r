# HAMILTONIAN MONTE CARLO UPDATE WITH POSSIBLE LOWER/UPPER BOUNDS.
#
# Radford M. Neal, 2012.
#
# Arguments:
#
#    lpr          Function returning the log probability of the position part 
#                 of the state, plus an arbitrary constant, with gradient
#                 as an attribute if grad=TRUE is passed.  May have attributes
#                 "upper" and/or "lower" giving bounds for variables.
#    initial      The initial position part of the state (a vector).
#    initial.p    The initial momentum part of the state (a vector), default
#                 is to use momentum variables generated as standard normals.
#    lpr.initial  The value of lpr(initial), maybe with grad.  May be omitted.
#    nsteps       Number of steps in trajectory used to propose a new state.
#                 (Default is 1, giving the "Langevin" method.)
#    step         Stepsize or stepsizes.  May be scalar or a vector of length 
#                 equal to the dimensionality of the state.
#    rand.step    Amount of random jitter for step. The step is multiplied
#                 by exp (runif (length(rand.step), -rand.step, rand.step)). 
#                 A scalar or vector of length equal to the dimensionality. 
#                 Default is 0.
#    return.traj  TRUE if values of q and p along the trajectory should be 
#                 returned (default is FALSE).
#
# The value returned is a list containing the following elements:
#
#    final        The new position part of the state
#    final.p      The new momentum part of the state
#    lpr          The value of lpr(final,grad=TRUE)
#    step         Stepsize(s) used (after jittering)
#    acc          1 if the proposal was accepted, 0 if rejected
#    apr          Acceptance probability of the proposal
#    delta        Change in H the acceptance decision was based on
#    reflections  Number of times a reflection off a boundary was done
#
# plus the following elements if return.traj was TRUE:
#
#    traj.q     Matrix of nsteps+1 rows with position values along trajectory
#    traj.p     Matrix of nsteps+1 rows with momentum values along trajectory
#    traj.H     Vector of length nsteps+1 with values of H along trajectory

bounded_hmc <- structure(handles.bounds=TRUE,
               function (lpr, initial, initial.p = rnorm(length(initial)), 
                         lpr.initial = NULL, nsteps = 1, step, rand.step = 0, 
                         return.traj = FALSE)
{
#print(rbind(lower=attr(lpr,"lower"),upper=attr(lpr,"upper")))
  # Check and process the arguments.

  process_nsteps_argument (nsteps)
  step <- process_step_arguments (length(initial), step, rand.step)

  # Get the lower and upper bounds.  Set lower.ix, upper.ix, and both.ix to
  # vectors of indexes of variables with only lower, only upper, or both lower 
  # and upper bounds.

  # browser()
  
  lower <- process_bounds (length(initial), attr(lpr,"lower"))
  upper <- process_bounds (length(initial), attr(lpr,"upper"))

  if (!is.null(lower) && !is.null(upper))
  { lower.present <- lower > -Inf
    upper.present <- upper < Inf
    lower.ix <- which (lower.present & !upper.present)
    upper.ix <- which (upper.present & !lower.present)
    both.ix <- which (lower.present & upper.present)
    if (any(lower[both.ix]>=upper[both.ix]))
    { stop("Lower bound not less than upper bound")
    }
  }
  else if (!is.null(lower))
  { lower.ix <- which (lower > -Inf)
    upper.ix <- both.ix <- NULL
  }
  else if (!is.null(upper))
  { upper.ix <- which (upper < Inf)
    lower.ix <- both.ix <- NULL
  }
  else
  { lower.ix <- upper.ix <- both.ix <- NULL
  }

  # Allocate space for the trajectory, if its return is requested.

  if (return.traj)
  { traj.q <- matrix(NA,nsteps+1,length(initial))
    traj.p <- matrix(NA,nsteps+1,length(initial))
    traj.H <- rep(NA,nsteps+1)
  }

  # Evaluate the log probability and gradient at the initial position, if not 
  # already known.

  if (is.null(lpr.initial) || is.null(attr(lpr.initial,"grad")))
  { lpr.initial <- lpr(initial,grad=TRUE)
  }

  # Compute the kinetic energy at the start of the trajectory.

  kinetic.initial <- sum(initial.p^2) / 2

  # Compute the trajectory by the leapfrog method.

  q <- initial
  p <- initial.p
  gr <- attr(lpr.initial,"grad")

  broken <- FALSE # will be true if breaks
  apr <- NULL #will  be replaced if doesn't break
  delta <- NULL #will be replaced if doesn't break
  
  
  if (return.traj)
  { traj.q[1,] <- q
    traj.p[1,] <- p
    traj.H[1] <- kinetic.initial - lpr.initial
  }

  reflections <- 0
  
  # Make a half step for momentum at the beginning.

  p <- p + (step/2) * gr

  # Alternate full steps for position and momentum.

  for (i in 1:nsteps)
  { 
    # Make a full step for the position.

    q <- q + step * p    
    
    # check if broken
    if (any(is.nan(c(p, q)))) {
      broken <- TRUE
      #stop("Algorithm breaks. Try a smaller epsilon.")
      break
    }

    # Check for bound violations, and adjust position and momentum for a 
    # reflection off the boundary.  For variables with both lower and upper
    # bounds, the bounds must be rechecked until there are no violations.

    for (k in lower.ix)
    { if (q[k]<lower[k])
      { q[k] <- lower[k] + (lower[k] - q[k])
        p[k] <- -p[k]
        reflections <- reflections + 1
      }
    }

    for (k in upper.ix)
    { if (q[k]>upper[k])
      { q[k] <- upper[k] - (q[k] - upper[k])
        p[k] <- -p[k]
        reflections <- reflections + 1
      }
    }

#     for (k in both.ix)
#     { repeat
#       { if (q[k]<lower[k])
#         { q[k] <- lower[k] + (lower[k] - q[k])
#           p[k] <- -p[k]
#           reflections <- reflections + 1
#         }
#         if (q[k]>upper[k])
#         { q[k] <- upper[k] - (q[k] - upper[k])
#           p[k] <- -p[k]
#           reflections <- reflections + 1
# 	}
#         else
#         { break
#         }
#       }
#     }
    
    for(k in both.ix) {
    while(q[k] < lower[k] || q[k] > upper[k]) {
      if (q[k]<lower[k]) {
        q[k] <- lower[k] + (lower[k] - q[k])
        p[k] <- -p[k]
        reflections <- reflections + 1
      } else if (q[k]>upper[k])
      { 
        q[k] <- upper[k] - (q[k] - upper[k])
      p[k] <- -p[k]
      reflections <- reflections + 1
      }
    }
    }  
    
    # for (k in lower.ix)
    # { cond1 <- tryCatch(q[k]<lower[k], error = function(e) "error") 
    # if(any(cond1 == "error")) browser()  
    # else if (cond1)
    # { q[k] <- lower[k] + (lower[k] - q[k])
    # p[k] <- -p[k]
    # reflections <- reflections + 1
    # }
    # }
    # 
    # for (k in upper.ix)
    # { cond2 <- tryCatch(q[k]>upper[k], error = function(e) "error") 
    # if(any(cond2 == "error")) browser()  
    # else if (cond2)
    # { q[k] <- upper[k] - (q[k] - upper[k])
    # p[k] <- -p[k]
    # reflections <- reflections + 1
    # }
    # }
    # 
    # for (k in both.ix)
    # { repeat
    # { cond1 <- tryCatch(q[k]<lower[k], error = function(e) "error") 
    # if(any(cond1 == "error")) browser()  
    # else if (cond1)
    # { q[k] <- lower[k] + (lower[k] - q[k])
    # p[k] <- -p[k]
    # reflections <- reflections + 1
    # }
    # cond2 <- tryCatch(q[k]>upper[k], error = function(e) "error") 
    # if(any(cond2 == "error")) browser()  
    # else if (cond2)
    #   { q[k] <- upper[k] - (q[k] - upper[k])
    #   p[k] <- -p[k]
    #   reflections <- reflections + 1
    #   }
    #   else
    #   { break
    #   }
    # }
    # }
    
    # check if broken
    if (any(is.nan(c(p, q)))) {
      broken <- TRUE
      #stop("Algorithm breaks. Try a smaller epsilon.")
      break
    }
    

    # Evaluate the gradient at the new position.

    lr <- lpr(q,grad=TRUE)
    gr <- attr(lr,"grad")

    # Record trajectory if asked to, with half-step for momentum.

    if (return.traj)
    { traj.q[i+1,] <- q
      traj.p[i+1,] <- p + (step/2) * gr
      traj.H[i+1] <- sum(traj.p[i+1,]^2)/2 - lr
    }
      
    # Make a full step for the momentum, except when we're coming to the end of 
    # the trajectory.  

    if (i!=nsteps)
    { p <- p + step * gr  
    }
  }

  
  # Make a half step for momentum at the end.

  p <- p + (step/2) * gr

  # Negate momentum at end of trajectory to make the proposal symmetric.

  p <- -p

  # Look at log probability and kinetic energy at the end of the trajectory.

  lpr.prop <- lr
  kinetic.prop <- sum(p^2) / 2

  # Accept or reject the state at the end of the trajectory.

  H.initial <- -lpr.initial + kinetic.initial
  H.prop <- -lpr.prop + kinetic.prop
  delta <- H.prop - H.initial
  apr <- min(1,exp(-delta))
  
  
  if (any(is.nan(c(p, q, apr)))) {
    broken <- TRUE
    apr <- NULL
    delta <- NULL
    #stop("Algorithm breaks. Try a smaller epsilon.")
  } 
    
  
  if(broken) # reject
  { 
    final.q <- initial
    final.p <- initial.p
    lpr.final <- lpr.initial
    acc <- 0
  }
  else if (runif(1) > apr) # reject
  { 
    final.q <- initial
    final.p <- initial.p
    lpr.final <- lpr.initial
    acc <- 0
  }
  else  # accept
  { 
    final.q <- q
    final.p <- p
    lpr.final <- lr
    acc <- 1
  }

  # Return new state, its log probability and gradient, plus additional
  # information, including the trajectory, if requested.

  r <- list (final=final.q, final.p=final.p, lpr=lpr.final, step=step, 
             apr=apr, acc=acc, delta=delta, reflections=reflections, broken = broken)

  if (return.traj) 
  { r$traj.q <- traj.q 
    r$traj.p <- traj.p
    r$traj.H <- traj.H 
  }

  r
})
