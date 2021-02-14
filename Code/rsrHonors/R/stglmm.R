stglmm <- function(fixed,
                   data,
                   locations,
                   times,
                   family,
                   covfn,
                   iter,
                   chains,
                   cores,
                   param_start,
                   priors,
                   tuning,
                   nu, # keep 2.5 for now
                   rank,
                   mul = 2){ # rank to reduce spatial matrix to

  start_time <- Sys.time()
  # Figure out how to interpret the formula fixed fixed
  model_frame <- lm(fixed, data = data, method = "model.frame")
  X <- model.matrix(fixed, data = model_frame)
  O <- model.response(model_frame)

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if (family$family == "gaussian") {
    stop("I haven't implemented gaussian families yet")
  } else if (family$family == "poisson") {
    dens_fun_log <- function(x, mean) {dpois(x, lambda = family$linkinv(mean), log = TRUE)}
  } else if (family$family == "binomial") {
    dens_fun_log <- function(x, mean) {dbinom(x, size = 1, prob = family$linkinv(mean), log = TRUE)}
  }

  if(class(locations) == "data.frame") {
    locations = as.matrix(locations)
  }


  # TODO: change this to param_start, then figure out how to get rid of it

  starting <- param_start
  ptm <- proc.time() # starting time
  p = ncol(X) # expected rank of X
  n <- length(O) # number of observations
  rk = rank*mul # rank of approximation matrix
  AP = chol2inv(chol(crossprod(X,X))) %*% t(X) # projection onto column space of X
  n_s <- nrow(locations)
  n_t <- length(time)

  #######################################
  # Prepare MCMC iteration values
  #######################################

  niter <- iter*chains
  message("\nMCMC chain size:", iter, "\n")

  #######################################
  # Define prior parameters
  #######################################
  beta.b <- priors[["beta.normal"]]
  s2.a   <- priors[["s2.IG"]][1]
  s2.b   <- priors[["s2.IG"]][2]
  phi.a  <- priors[["phi.Unif"]][1]
  phi.b  <- priors[["phi.Unif"]][2]


  ###################################
  # Orig: define index for parameters
  # Must mean that parameters are stored in a single vector
  # This vector is called sParams
  ####################################
  nParams <- p + 2;
  beta_index <- 1:p;
  sigma2_index <- p + 1;
  phi_index <- sigma2_index + 1
  delta_index <- (nParams + 1):(nParams + rank)
  w_index <- (nParams + rank + 1):(nParams + rank + n)

  #########################################
  # Orig: initialize some matrix for storage
  # Prepares matrices to hold estimates of parameters
  #########################################
  samples_eta <- matrix(NA,ncol = rank,nrow = niter) # store the r.e samples
  samples_s <- matrix(NA,ncol = nParams,nrow = niter) # store the parameter samples
  samples_w <- matrix(NA,ncol = n,nrow = niter) # store the r.e samples
  sTunings <- sParams <- matrix(NA,nrow = nParams) # update the current parameters for each iteration


  samples <- array(dim = c(iter, chains, nParams + rank + n),
                   dimnames = list("Iteration" = 1:iter,
                                   "Chain" = 1:chains,
                                   "Parameter" = c(paste0("beta_", 0:(p - 1)),
                                                   "sigma2",
                                                   "phi",
                                                   paste0("delta_", 1:rank),
                                                   paste0("w_", 1:n))))


  current_beta <- rep(NA_real_, p)
  current_sigma2 <- NA_real_
  current_phi <- NA_real_
  current_delta <- rep(NA_real_, rank)
  current_w <- rep(NA_real_, rank)

  #########################################
  # Orig: feed in the initial values and tuning params
  # Set starting parameters
  #########################################


  sTunings[beta_index] <- tuning[["beta"]]
  sTunings[sigma2_index] <- tuning[["s2"]]
  sTunings[phi_index] <- tuning[["phi"]]


  wTunings <- rep(tuning[["w"]],rank)
  wTunings <- log(wTunings)
  sTunings <- log(sTunings)



  # Stores info on acceptance rates for MCMC
  # for acceptance rate
  accept_s <- matrix(NA, ncol = nParams, nrow = chains)
  accept_w <- matrix(NA, ncol = rank, nrow = chains)



  ##### start MCMC loop #####
  for (k in 1:chains) {
    sAccepts <- matrix(0,nrow = nParams)
    wAccepts <- matrix(0,nrow = rank)

    # Generate Chain Starting positions following stan recomendation

    current_beta <- runif(p, min = -2, max = 2)
    current_sigma2 <- exp(runif(1, min = -2, max = 2))
    current_phi <- phi.a + (phi.b - phi.a)/(1 + exp(-runif(1, min = -2, max = 2)))



    est_start  <- Sys.time()
    # this is where the first call to the cpp code occurs. It is for the random projections part of the algorithm.
    K = rp(
      current_phi, # single number, phi in starting param list
      locations, #literally x,y locations of obs
      n_s, # number if interations (iter)
      as.integer(rk), # rank times mul (what is mul?)
      nu, #nu as before
      as.integer(cores),# number of cores
      cov_fun = 0
    ) # C++ function for approximating eigenvectors
    est_time  <- Sys.time() - est_start # calculates how long one random projection takes
    message("Estimated time (hrs):",round((chains - k + 1)*iter*2*est_time/3600, 3) ,"\n") # prints out estimate time in hours

    K.rp = list(d = K[[1]],u = K[[2]][,1:rank]) # d is sing values, u is left sing vecs, only take first 1:rank
    d <- (K.rp$d[1:rank])^2 # approximated eigenvalues
    U1 <- U <- u <- K.rp$u[,1:rank] # approximated eigenvectors

    current_delta <- runif(rank, min = -2, max = 2)
    current_w <- U %*% (sqrt(d)*current_delta)


    xbeta <- X %*% current_beta



    for (i in 1:iter) {


      #########################################################
      # Block Update Beta
      #########################################################

      # propose from normal
      beta_proposal <- rnorm(p, current_beta, sd = exp(sTunings[beta_index])) # draw for proposed beta params

      # calculate likelihood ratio
      beta_current_likelihood <- beta_log_full_conditional_st(current_beta, current_w = current_w, X = X,
                                                           O = O,
                                                           beta.b = beta.b,
                                                           p = p,
                                                           dens_fun_log = dens_fun_log,
                                                           n_t = n_t)

      beta_proposal_likelihood <- beta_log_full_conditional_st(beta_proposal, current_w = current_w, X = X,
                                                            O = O,
                                                            beta.b = beta.b,
                                                            p = p,
                                                            dens_fun_log = dens_fun_log,
                                                            n_t = n_t)
      beta_lr <- beta_proposal_likelihood - beta_current_likelihood # log likelihood ratio log(proposed_likelihood/current_likelihood)


      # accept or reject proposal
      if (log(runif(1)) < beta_lr) { # if likelihood ratio is greater than runif(1), accept
        current_beta <- beta_proposal # update params with new guess
        sAccepts[beta_index] <- sAccepts[beta_index] + 1 # mark we accepted
        xbeta <- X %*% current_beta # update xbeta, which will be used in next set of estimation
      }

      #########################################################
      # Block Update Phi
      #########################################################

      # propose from normal
      phi_proposal <-  rnorm(1, current_phi, sd = exp(sTunings[phi_index]))
      phi_proposal_likelihood <- phi_current_likelihood <- phi_sp_log_full_conditional(current_phi, locations = locations, xbeta = xbeta, current_delta = current_delta, U1 = U1, # data and params
                                                                                    O = O, # observations
                                                                                    current_sigma2 = current_sigma2, # priors
                                                                                    nu = nu, n_s = n_s, n_t = n_t, rk = rk, cores = cores, rank = rank, # control params
                                                                                    dens_fun_log = dens_fun_log)



      # checks if guess in bounds of Unif(a, b) prior
      if (phi_proposal < phi.b & phi_proposal > phi.a) {
        phi_proposal_likelihood <- phi_sp_log_full_conditional(phi_proposal, locations = locations, xbeta = xbeta, current_delta = current_delta, U1 = U1, # data and params
                                                            O = O, # observations
                                                            current_sigma2 = current_sigma2, # priors
                                                            nu = nu, n_s = n_s, n_t = n_t, rk = rk, cores = cores, rank = rank, # control params
                                                            dens_fun_log = dens_fun_log)
      } else {
        phi_proposal_likelihood$likelihood <- -Inf
      }
      phi_lr <- phi_proposal_likelihood$likelihood - phi_current_likelihood$likelihood

      if (log(runif(1)) < phi_lr) {
        current_phi <- phi_proposal
        sAccepts[phi_index] <- sAccepts[phi_index] + 1
        phi_current_likelihood <- phi_proposal_likelihood
        d <- phi_current_likelihood$d # approximated eigenvalues^2
        U <- phi_current_likelihood$U # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
        u <- phi_current_likelihood$u # approximated eigenvectors
      }

      twKinvw <- phi_current_likelihood$twKinvw # current delta cross product for some reason

      #########################################################
      # Block Update Sigma 2
      #########################################################

      sigma2_proposal <- rnorm(1, current_sigma2, sd = exp(sTunings[sigma2_index]))
      if (sigma2_proposal > 0) {
        sigma2_current_likelihood <- sigma2_log_full_conditional(sigma2 = current_sigma2, twKinvw = twKinvw,
                                                                 s2.a = s2.a, s2.b = s2.b,
                                                                 rank = rank)

        sigma2_proposal_likelihood <-  sigma2_log_full_conditional(sigma2 = sigma2_proposal, twKinvw = twKinvw,
                                                                   s2.a = s2.a, s2.b = s2.b,
                                                                   rank = rank)
        sigma2_lr <- sigma2_proposal_likelihood - sigma2_current_likelihood
      } else {
        sigma2_lr <- -Inf
      }

      #########################################################
      # Block Update Delta
      #########################################################

      if (log(runif(1)) < sigma2_lr) {
        current_sigma2 <- sigma2_proposal
        sAccepts[sigma2_index] <- sAccepts[sigma2_index] + 1
      }

      # update random effects using multivariate random walk with spherical normal proposal
      delta_proposal <- rnorm(rank, current_delta, sd = exp(wTunings))

      delta_proposal_likelihood <- delta_log_full_conditional_st(delta = delta_proposal, xbeta = xbeta,  U = U, d = d,
                                                              O = O,
                                                              n_t = n_t,
                                                              current_sigma2 = current_sigma2,
                                                              dens_fun_log = dens_fun_log)

      delta_current_likelihood <- delta_log_full_conditional_st(delta = current_delta, xbeta = xbeta,  U = U, d = d,
                                                             O = O,
                                                             n_t = n_t,
                                                             current_sigma2 = current_sigma2,
                                                             dens_fun_log = dens_fun_log)
      delta_lr <- delta_proposal_likelihood$lr - delta_current_likelihood$lr

      if (log(runif(1)) < delta_lr) {
        current_delta <- delta_proposal
        wAccepts <- wAccepts + 1
        delta_current_likelihood <- delta_proposal_likelihood
        current_w <- delta_current_likelihood$w
      }

      samples_s[(k - 1)*iter + i,] <- c(current_beta, current_sigma2, current_phi)
      samples_w[(k - 1)*iter + i,] <- current_w
      samples_eta[(k - 1)*iter + i,] <- current_delta

      samples[i, k, 1:nParams] <- c(current_beta, current_sigma2, current_phi) # stores fixed effects draw and variance and spatial range draw
      samples[i, k, (nParams + 1):(nParams + rank)] <- current_delta # stores draw for synthetic variable of spatial effects
      samples[i, k, (nParams + 1 + rank):(nParams + rank + n)] <- current_w # stores draw of spatial effect (computed from delta draw)
    }

    output_progress(beta_draws = samples[, k, beta_index], sAccepts = sAccepts,
                    iter = iter, chain = k)

    accept_s[k,] <- sAccepts/(iter)
    accept_w[k,] <- wAccepts/(iter)

  }
  accept_s <- apply(accept_s,2,mean)
  accept_w <- apply(accept_w,2,mean)

  stop_time <- Sys.time()

  runtime = stop_time - start_time
  message("runtime", "\n")
  message(round(as.numeric(runtime, units = "mins"), 3), " minutes\n")


  # return(structure(list(coefficients = list(beta = , # median estimates
  #                                           sigma2 = ,
  #                                           phi = ,
  #                                           delta = ,
  #                                           w = ),
  #                       adj_coefficeints = , # adjusted (only beta)
  #                       ses = list(beta = , #mad
  #                                  sigma2 = ,
  #                                  phi = ,
  #                                  delta = ,
  #                                  w = ),
  #                       residuals = ,
  #                       fitted.values = ,
  #                       linear.predictors = ,
  #                       covmat = ,
  #                       family = family,
  #                       formula_fixed = fixed,
  #                       formula_spatial = spatial,
  #                       prior_info = priors,
  #                       tuning_info = tuning,
  #                       samples = samples,
  #                       log_likelihood = ,# log likelihood array for loo?
  #                       model = "rrp_sglmm",
  #                       accept = list(beta = ,
  #                                     sigma2 = ,
  #                                     phi = ,
  #                                     delta = ,
  #                                     w = )),
  #                  class = "rsrHonors_rrp_sglmm"))



  return(list(run.time = runtime,
              p.params = samples[, , 1:nParams],
              model = "rrp/arrp",
              rank = rank,
              eta.params = samples_eta,
              w.params = samples_w,
              accepts = accept_s,
              acceptw = accept_w,
              samples = samples))
}

# --------------------------------------------------------
# functions from source code:
# --------------------------------------------------------

output_progress <- function(beta_draws, sAccepts,
                            iter, chain) {
  message("Chain ",chain,"\n")
  message("-------------------------------------------------------\n")
  message("----------------parameter estimates--------------------\n")
  message(paste0(round(colMeans(beta_draws), 3), collapse = "   "),"\n")
  message("-------------------------------------------------------\n")
  message("----------------acceptance rate------------------------\n")
  message(paste0(round(sAccepts/(iter),3), collapse = "   "),"\n")
  message("-------------------------------------------------------\n")
}

beta_log_full_conditional_st <- function(beta, current_w, X,
                                      O,
                                      beta.b,
                                      p,
                                      dens_fun_log,
                                      n_t){

  z <- X %*% beta + rep(current_w, times = n_t) # if rsr, current_w = L%*%eta
  lf <- sum(dens_fun_log(O, mean = z)) - crossprod(beta)/(p*beta.b)
  return(lf)

}

delta_log_full_conditional_st <- function(delta, xbeta,  U, d,
                                       O,
                                       n_t,
                                       current_sigma2,
                                       dens_fun_log){ # delta is rank-m
  w = U %*% (sqrt(d)*delta)
  z <- xbeta + rep(w, n_t)
  foo2 <- crossprod(delta,delta) # d = D^2 from random projection
  lf <- sum(dens_fun_log(O, mean = z)) - 1/(2*current_sigma2) * foo2
  return(list(lr = lf, twKinvw = foo2, w = w))
}

phi_sp_log_full_conditional <- function(phi, locations, xbeta, current_delta, U1, # data and params
                                     O, # observations
                                     current_sigma2, # priors
                                     nu, n_s, n_t, rk, cores, rank, # control params
                                     dens_fun_log){ # density function
  K1 = rp(phi,
          locations,
          n_s,
          as.integer(rk),
          nu,
          as.integer(cores)) # C++ function for approximating eigenvectors

  K.rp = list(d = K1[[1]],u = K1[[2]][,1:rank])
  d <- (K.rp$d[1:rank])^2 # approximated eigenvalues

  u <- K.rp$u[,1:rank]

  signdiag = sign(diag(t(u) %*% U1)) # find the sign change
  signdiag = as.logical(1 - signdiag)
  u[,signdiag] = -u[,signdiag]

  U <- u
  z <- xbeta + rep(U %*% (sqrt(d)*current_delta), n_t) # U is from random projection
  foo2 <- crossprod(current_delta, current_delta)
  likelihood <- (
    sum(dens_fun_log(O, mean = z)) - 0.5*1/current_sigma2 * foo2 # likelihood
    # + log(phi - phi.a) + log(phi.b - phi) # jacobian
  )
  return(list(likelihood = likelihood, d = d, twKinvw = foo2, U = U, u = u))
}
# (-s2.a - 1 - rank/2)*log(current_sigma2) - (s2.b + 0.5*twKinvw)/current_sigma2
sigma2_sp_log_full_conditional <- function(sigma2, twKinvw,
                                        s2.a, s2.b,
                                        rank) {
  (-s2.a - 1 - rank/2)*log(sigma2) - (s2.b + 0.5*twKinvw)/sigma2
}

alpha_log_full_conditional <- function(delta, xbeta,  U, d,
                                       O,
                                       n_s,
                                       current_sigma2,
                                       dens_fun_log){ # delta is rank-m
  v = U %*% (sqrt(d)*delta)
  z <- xbeta + rep(v, each = n_s)
  foo2 <- crossprod(delta,delta) # d = D^2 from random projection
  lf <- sum(dens_fun_log(O, mean = z)) - 1/(2*current_sigma2) * foo2
  return(list(lr = lf, twKinvw = foo2, v = v))
}

phi_tm_log_full_conditional <- function(phi, times, xbeta, current_delta, U1, # data and params
                                        O, # observations
                                        current_sigma2, # priors
                                        nu, n_s, n_t, rk, cores, rank, # control params
                                        dens_fun_log){ # density function

  K1 = rp(phi, times, n_t, as.integer(rk), nu,as.integer(cores)) # C++ function for approximating eigenvectors
  K.rp = list(d = K1[[1]],u = K1[[2]][,1:rank])
  d <- (K.rp$d[1:rank])^2 # approximated eigenvalues

  u <- K.rp$u[,1:rank]

  signdiag = sign(diag(t(u) %*% U1)) # find the sign change
  signdiag = as.logical(1 - signdiag)
  u[,signdiag] = -u[,signdiag]

  U <- u # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
  z <- xbeta + rep(U %*% (sqrt(d)*current_delta), each = n_s) # U is from random projection
  foo2 <- crossprod(current_delta, current_delta)
  likelihood <- (
    sum(dens_fun_log(O, mean = z)) - 0.5*1/current_sigma2 * foo2 # likelihood
    # + log(phi - phi.a) + log(phi.b - phi) # jacobian
  )
  return(list(likelihood = likelihood, d = d, twKinvw = foo2, U = U, u = u))
}
# (-s2.a - 1 - rank/2)*log(current_sigma2) - (s2.b + 0.5*twKinvw)/current_sigma2
sigma2_tm_log_full_conditional <- function(sigma2, twKinvw,
                                           s2.a, s2.b,
                                           rank) {
  (-s2.a - 1 - rank/2)*log(sigma2) - (s2.b + 0.5*twKinvw)/sigma2
}



