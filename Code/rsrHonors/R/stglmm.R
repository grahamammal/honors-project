stglmm <- function(fixed,
                   data,
                   locations,
                   times,
                   family,
                   covfn,
                   iter,
                   chains,
                   cores,
                   priors,
                   tuning,
                   nu, # keep 2.5 for now
                   rank_sp,
                   rank_tm){ # rank to reduce spatial matrix to

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

  if (class(locations) == "data.frame") {
    locations = as.matrix(locations)
  }


  dist_space <- as.matrix(dist(locations))
  dist_time <- as.matrix(dist(times))


  ptm <- proc.time() # starting time
  p = ncol(X) # expected rank of X
  n <- length(O) # number of observations
  n_s <- nrow(locations)
  n_t <- length(times)

  #######################################
  # Prepare MCMC iteration values
  #######################################

  niter <- iter*chains
  message("\nMCMC chain size:", iter, "\n")

  #######################################
  # Define prior parameters
  #######################################
  beta.b <- priors[["beta.normal"]]
  s2_sp_a   <- priors[["s2_sp_IG"]][1]
  s2_sp_b   <- priors[["s2_sp_IG"]][2]
  phi_sp_a  <- priors[["phi_sp_unif"]][1]
  phi_sp_b  <- priors[["phi_sp_unif"]][2]
  s2_tm_a   <- priors[["s2_tm_IG"]][1]
  s2_tm_b   <- priors[["s2_tm_IG"]][2]
  phi_tm_a  <- priors[["phi_tm_unif"]][1]
  phi_tm_b  <- priors[["phi_tm_unif"]][2]


  ###################################
  # Orig: define index for parameters
  # Must mean that parameters are stored in a single vector
  # This vector is called sParams
  ####################################
  nParams <- p + 4;
  beta_index <- 1:p;
  sigma2_sp_index <- max(beta_index) + 1;
  sigma2_tm_index <- sigma2_sp_index + 1;
  phi_sp_index <- sigma2_tm_index + 1
  phi_tm_index <- phi_sp_index + 1
  delta_index <- (phi_tm_index + 1):(phi_tm_index + rank_sp)
  w_index <- (max(delta_index) + 1):(max(delta_index) + n_s)
  alpha_index <- (max(w_index) + 1):(max(w_index) + rank_tm)
  v_index <- (max(alpha_index) + 1):(max(alpha_index) + n_t)


  #########################################
  # Orig: initialize some matrix for storage
  # Prepares matrices to hold estimates of parameters
  #########################################
  samples_eta <- matrix(NA,ncol = rank_sp, nrow = niter) # store the r.e samples
  samples_w <- matrix(NA,ncol = n,nrow = niter) # store the r.e samples
  sTunings <- sParams <- matrix(NA,nrow = nParams) # update the current parameters for each iteration

  samples <- array(dim = c(iter, chains, max(v_index)),
                   dimnames = list("Iteration" = 1:iter,
                                   "Chain" = 1:chains,
                                   "Parameter" = c(paste0("beta_", 0:(p - 1)),
                                                   "sigma2_sp",
                                                   "phi_sp",
                                                   "sigma2_tm",
                                                   "phi_tm",
                                                   paste0("delta_", 1:length(delta_index)),
                                                   paste0("w_", 1:length(w_index)),
                                                   paste0("alpha_", 1:length(alpha_index)),
                                                   paste0("v_", 1:length(v_index)))))


  current_beta <- rep(NA_real_, p)
  current_sigma2_sp <- NA_real_
  current_phi_sp <- NA_real_
  current_sigma2_tm <- NA_real_
  current_phi_tm <- NA_real_
  current_delta <- rep(NA_real_, rank_sp)
  current_w <- rep(NA_real_, n_s)
  current_alpha <- rep(NA_real_, rank_tm)
  current_v <- rep(NA_real_, n_t)


  #########################################
  # Orig: feed in the initial values and tuning params
  # Set starting parameters
  #########################################


  sTunings[beta_index] <- tuning[["beta"]]
  sTunings[sigma2_sp_index] <- tuning[["s2_sp"]]
  sTunings[phi_sp_index] <- tuning[["phi_sp"]]
  sTunings[sigma2_tm_index] <- tuning[["s2_tm"]]
  sTunings[phi_tm_index] <- tuning[["phi_tm"]]

  deltaTunings <- rep(tuning[["delta"]],rank_sp)
  alphaTunings <- rep(tuning[["alpha"]],rank_tm)




  # Stores info on acceptance rates for MCMC
  # for acceptance rate
  accept_s <- matrix(NA, ncol = nParams, nrow = chains)
  accept_w <- matrix(NA, ncol = rank_sp, nrow = chains)



  ##### start MCMC loop #####
  for (k in 1:chains) {
    sAccepts <- matrix(0,nrow = nParams)
    wAccepts <- matrix(0,nrow = rank_sp)
    vAccepts <- matrix(0, nrow = rank_tm)
    # Generate Chain Starting positions following stan recomendation

    current_beta <- runif(p, min = -2, max = 2)
    current_sigma2_sp <- exp(runif(1, min = -2, max = 2))
    current_phi_sp <- phi_sp_a + (phi_sp_b - phi_sp_a)/(1 + exp(-runif(1, min = -2, max = 2)))

    current_sigma2_tm <- exp(runif(1, min = -2, max = 2))
    current_phi_tm <- phi_tm_a + (phi_tm_b - phi_tm_b) / (1 + exp(-runif(1, min = -2, max = 2)))

    est_start  <- Sys.time()
    # this is where the first call to the cpp code occurs. It is for the random projections part of the algorithm.
    K = rp(
      current_phi_sp, # single number, phi in starting param list
      dist_space, #literally x,y locations of obs
      n_s, # number of points
      rank_sp,
      nu, #nu as before
      as.integer(cores),# number of cores
      cov_fun = 0
    ) # C++ function for approximating eigenvectors

    K.rp = list(d = K[[1]],u = K[[2]][,1:rank_sp]) # d is sing values, u is left sing vecs, only take first 1:rank
    d <- (K.rp$d[1:rank_sp])^2 # approximated eigenvalues
    U1 <- U <- u <- K.rp$u[,1:rank_sp] # approximated eigenvectors

    current_delta <- runif(rank_sp, min = -2, max = 2)
    current_w <- U %*% (sqrt(d)*current_delta)

    J = rp(
      current_phi_tm, # single number, phi in starting param list
      dist_time, #literally x,y locations of obs
      n_t, # number of points
      rank_tm,
      nu, #nu as before
      as.integer(cores),# number of cores
      cov_fun = 1
    ) # C++ function for approximating eigenvectors

    est_time  <- Sys.time() - est_start # calculates how long one random projection takes
    message("Estimated time (hrs):",round((chains - k + 1)*iter*2*est_time/3600, 3) ,"\n") # prints out estimate time in hours


    J.rp = list(l = J[[1]], v = J[[2]][,1:rank_tm]) # d is sing values, u is left sing vecs, only take first 1:rank
    l <- (J.rp$l[1:rank_tm])^2 # approximated eigenvalues
    V1 <- V <- v <- J.rp$v[,1:rank_tm] # approximated eigenvectors

    current_alpha <- runif(rank_tm, min = -2, max = 2)
    current_v <-  V %*% (sqrt(l)*current_alpha)



    xbeta <- X %*% current_beta



    for (i in 1:iter) {


      #########################################################
      # Block Update Beta
      #########################################################

      # propose from normal
      beta_proposal <- rnorm(p, current_beta, sd = sTunings[beta_index]) # draw for proposed beta params

      # calculate likelihood ratio
      beta_current_likelihood <- beta_log_full_conditional_st(current_beta, current_w = current_w, current_v = current_v, X = X,
                                                           O = O,
                                                           beta.b = beta.b,
                                                           p = p,
                                                           dens_fun_log = dens_fun_log,
                                                           n_t = n_t, n_s = n_s)

      beta_proposal_likelihood <- beta_log_full_conditional_st(beta_proposal, current_w = current_w, current_v = current_v, X = X,
                                                            O = O,
                                                            beta.b = beta.b,
                                                            p = p,
                                                            dens_fun_log = dens_fun_log,
                                                            n_t = n_t, n_s = n_s)
      beta_lr <- beta_proposal_likelihood - beta_current_likelihood # log likelihood ratio log(proposed_likelihood/current_likelihood)


      # accept or reject proposal
      if (log(runif(1)) < beta_lr) { # if likelihood ratio is greater than runif(1), accept
        current_beta <- beta_proposal # update params with new guess
        sAccepts[beta_index] <- sAccepts[beta_index] + 1 # mark we accepted
        xbeta <- X %*% current_beta # update xbeta, which will be used in next set of estimation
      }

      #########################################################
      # Block Update Phi Spatial
      #########################################################

      # propose from normal
      phi_sp_proposal <-  rnorm(1, current_phi_sp, sd = sTunings[phi_sp_index])
      phi_sp_proposal_likelihood <- phi_sp_current_likelihood <- phi_sp_log_full_conditional(current_phi_sp, dist_space = dist_space, xbeta = xbeta, current_delta = current_delta, U1 = U1, current_v = current_v,# data and params
                                                                                    O = O, # observations
                                                                                    current_sigma2_sp = current_sigma2_sp, # priors
                                                                                    nu = nu, n_s = n_s, n_t = n_t, rank = rank_sp, cores = cores, # control params
                                                                                    dens_fun_log = dens_fun_log)



      # checks if guess in bounds of Unif(a, b) prior
      if (phi_sp_proposal < phi_sp_b & phi_sp_proposal > phi_sp_a) {
        phi_sp_proposal_likelihood <- phi_sp_log_full_conditional(phi_sp_proposal, dist_space = dist_space, xbeta = xbeta, current_delta = current_delta, U1 = U1, current_v = current_v, # data and params
                                                            O = O, # observations
                                                            current_sigma2_sp = current_sigma2_sp, # priors
                                                            nu = nu, n_s = n_s, n_t = n_t, rank = rank_sp, cores = cores, # control params
                                                            dens_fun_log = dens_fun_log)
      } else {
        phi_sp_proposal_likelihood$likelihood <- -Inf
      }
      phi_sp_lr <- phi_sp_proposal_likelihood$likelihood - phi_sp_current_likelihood$likelihood

      if (log(runif(1)) < phi_sp_lr) {
        current_phi_sp <- phi_sp_proposal
        sAccepts[phi_sp_index] <- sAccepts[phi_sp_index] + 1
        phi_sp_current_likelihood <- phi_sp_proposal_likelihood
        d <- phi_sp_current_likelihood$d # approximated eigenvalues^2
        U <- phi_sp_current_likelihood$U # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
        u <- phi_sp_current_likelihood$u # approximated eigenvectors
      }

      twKinvw <- phi_sp_current_likelihood$twKinvw # current delta cross product for some reason

      #########################################################
      # Block Update Sigma 2 Spatial
      #########################################################

      sigma2_sp_proposal <- rnorm(1, current_sigma2_sp, sd = sTunings[sigma2_sp_index])
      if (sigma2_sp_proposal > 0) {
        sigma2_sp_current_likelihood <- sigma2_log_full_conditional(sigma2 = current_sigma2_sp, twKinvw = twKinvw,
                                                                 s2_a = s2_sp_a, s2_b = s2_sp_b,
                                                                 rank = rank_sp)

        sigma2_sp_proposal_likelihood <-  sigma2_log_full_conditional(sigma2 = sigma2_sp_proposal, twKinvw = twKinvw,
                                                                   s2_a = s2_sp_a, s2_b = s2_sp_b,
                                                                   rank = rank_sp)
        sigma2_sp_lr <- sigma2_sp_proposal_likelihood - sigma2_sp_current_likelihood
      } else {
        sigma2_sp_lr <- -Inf
      }

      if (log(runif(1)) < sigma2_sp_lr) {
        current_sigma2_sp <- sigma2_sp_proposal
        sAccepts[sigma2_sp_index] <- sAccepts[sigma2_sp_index] + 1
      }

      #########################################################
      # Block Update Delta
      #########################################################



      # update random effects using multivariate random walk with spherical normal proposal
      delta_proposal <- rnorm(rank_sp, current_delta, sd = deltaTunings)

      delta_proposal_likelihood <- delta_log_full_conditional_st(delta = delta_proposal, xbeta = xbeta,  U = U, d = d, current_v = current_v,
                                                              O = O,
                                                              n_t = n_t, n_s = n_s,
                                                              current_sigma2_sp = current_sigma2_sp,
                                                              dens_fun_log = dens_fun_log)

      delta_current_likelihood <- delta_log_full_conditional_st(delta = current_delta, xbeta = xbeta,  U = U, d = d, current_v = current_v,
                                                             O = O,
                                                             n_t = n_t, n_s = n_s,
                                                             current_sigma2_sp = current_sigma2_sp,
                                                             dens_fun_log = dens_fun_log)
      delta_lr <- delta_proposal_likelihood$lr - delta_current_likelihood$lr

      if (log(runif(1)) < delta_lr) {
        current_delta <- delta_proposal
        wAccepts <- wAccepts + 1
        delta_current_likelihood <- delta_proposal_likelihood
        current_w <- delta_current_likelihood$w
      }

      #########################################################
      # Block Update Phi Temporal
      #########################################################

      # propose from normal
      phi_tm_proposal <-  rnorm(1, current_phi_tm, sd = sTunings[phi_tm_index])
      phi_tm_proposal_likelihood <- phi_tm_current_likelihood <- phi_tm_log_full_conditional(current_phi_tm, dist_time = dist_time, xbeta = xbeta, current_alpha = current_alpha, V1 = V1, current_w = current_w,# data and params
                                                                                             O = O, # observations
                                                                                             current_sigma2_tm = current_sigma2_tm, # priors
                                                                                             nu = nu, n_s = n_s, n_t = n_t, rank = rank_tm, cores = cores, # control params
                                                                                             dens_fun_log = dens_fun_log)


      # checks if guess in bounds of Unif(a, b) prior
      if (phi_tm_proposal < phi_tm_b & phi_tm_proposal > phi_tm_a) {
        phi_tm_proposal_likelihood <- phi_tm_log_full_conditional(phi_tm_proposal, dist_time = dist_time, xbeta = xbeta, current_alpha = current_alpha, V1 = V1, current_w = current_w, # data and params
                                                                  O = O, # observations
                                                                  current_sigma2_tm = current_sigma2_tm, # priors
                                                                  nu = nu, n_s = n_s, n_t = n_t, rank = rank_tm, cores = cores, # control params
                                                                  dens_fun_log = dens_fun_log)
      } else {
        phi_tm_proposal_likelihood$likelihood <- -Inf
      }
      phi_tm_lr <- phi_tm_proposal_likelihood$likelihood - phi_tm_current_likelihood$likelihood

      if (log(runif(1)) < phi_tm_lr) {
        current_phi_tm <- phi_tm_proposal
        sAccepts[phi_tm_index] <- sAccepts[phi_tm_index] + 1
        phi_tm_current_likelihood <- phi_tm_proposal_likelihood
        l <- phi_tm_current_likelihood$l # approximated eigenvalues^2
        V <- phi_tm_current_likelihood$V # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
        v <- phi_tm_current_likelihood$v # approximated eigenvectors
      }

      twKinvw <- phi_tm_current_likelihood$twKinvw # current delta cross product for some reason



      #########################################################
      # Block Update Sigma 2 Temporal
      #########################################################

      sigma2_tm_proposal <- rnorm(1, current_sigma2_tm, sd = sTunings[sigma2_tm_index])
      if (sigma2_tm_proposal > 0) {
        sigma2_tm_current_likelihood <- sigma2_log_full_conditional(sigma2 = current_sigma2_tm, twKinvw = twKinvw,
                                                                    s2_a = s2_tm_a, s2_b = s2_tm_b,
                                                                    rank = rank_tm)

        sigma2_tm_proposal_likelihood <-  sigma2_log_full_conditional(sigma2 = sigma2_tm_proposal, twKinvw = twKinvw,
                                                                      s2_a = s2_tm_a, s2_b = s2_tm_b,
                                                                      rank = rank_tm)
        sigma2_tm_lr <- sigma2_tm_proposal_likelihood - sigma2_tm_current_likelihood
      } else {
        sigma2_tm_lr <- -Inf
      }

      if (log(runif(1)) < sigma2_tm_lr) {
        current_sigma2_tm <- sigma2_tm_proposal
        sAccepts[sigma2_tm_index] <- sAccepts[sigma2_tm_index] + 1
      }

      #########################################################
      # Block Update Alpha
      #########################################################

      # update random effects using multivariate random walk with spherical normal proposal
      alpha_proposal <- rnorm(rank_tm, current_alpha, sd = alphaTunings)

      alpha_proposal_likelihood <- alpha_log_full_conditional(alpha = alpha_proposal, xbeta = xbeta,  V = V, l = l,
                                                                 O = O,
                                                                 n_s = n_s,
                                                                 current_sigma2_tm = current_sigma2_tm,
                                                                 dens_fun_log = dens_fun_log)

      alpha_current_likelihood <- alpha_log_full_conditional(alpha = current_alpha, xbeta = xbeta,  V = V, l = l,
                                                                O = O,
                                                                n_s = n_s,
                                                                current_sigma2_tm = current_sigma2_tm,
                                                                dens_fun_log = dens_fun_log)
      alpha_lr <- alpha_proposal_likelihood$lr - alpha_current_likelihood$lr

      if (log(runif(1)) < alpha_lr) {
        current_alpha <- alpha_proposal
        vAccepts <- vAccepts + 1
        alpha_current_likelihood <- alpha_proposal_likelihood
        current_v <- alpha_current_likelihood$v
      }


      samples_w[(k - 1)*iter + i,] <- current_w
      samples_eta[(k - 1)*iter + i,] <- current_delta

      samples[i, k, c(beta_index, sigma2_sp_index, phi_sp_index, sigma2_tm_index, phi_tm_index)] <- c(current_beta, current_sigma2_sp, current_phi_sp, current_sigma2_tm, current_phi_tm) # stores fixed effects draw and variance and spatial range draw
      samples[i, k, delta_index] <- current_delta # stores draw for synthetic variable of spatial effects
      samples[i, k, w_index] <- current_w # stores draw of spatial effect (computed from delta draw)
      samples[i, k, alpha_index] <- current_alpha
      samples[i, k, v_index] <- current_v
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
  #                                           sigma2_sp = ,
  #                                           phi = ,
  #                                           delta = ,
  #                                           w = ),
  #                       adj_coefficeints = , # adjusted (only beta)
  #                       ses = list(beta = , #mad
  #                                  sigma2_sp = ,
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
  #                                     sigma2_sp = ,
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

beta_log_full_conditional_st <- function(beta, current_w, current_v, X,
                                      O,
                                      beta.b,
                                      p,
                                      dens_fun_log,
                                      n_t, n_s){

  z <- X %*% beta + rep(current_w, times = n_t) + rep(current_v, times = n_s) # if rsr, current_w = L%*%eta
  lf <- sum(dens_fun_log(O, mean = z)) - crossprod(beta)/(p*beta.b)
  return(lf)

}

delta_log_full_conditional_st <- function(delta, xbeta,  U, d,
                                       O,
                                       current_v,
                                       n_t, n_s,
                                       current_sigma2_sp,
                                       dens_fun_log){ # delta is rank-m
  w = U %*% (sqrt(d)*delta)
  z <- xbeta + rep(w, n_t) + rep(current_v, each = n_s)
  foo2 <- crossprod(delta,delta) # d = D^2 from random projection
  lf <- sum(dens_fun_log(O, mean = z)) - 1/(2*current_sigma2_sp) * foo2
  return(list(lr = lf, twKinvw = foo2, w = w))
}

phi_sp_log_full_conditional <- function(phi, dist_space, xbeta, current_delta, U1, # data and params
                                     O, # observations
                                     current_v,
                                     current_sigma2_sp,
                                     nu, n_s, n_t, cores, rank, # control params
                                     dens_fun_log){ # density function
  K1 = rp(phi,
          dist_space,
          n_s,
          rank,
          nu,
          as.integer(cores),
          cov_fun = 0) # C++ function for approximating eigenvectors

  K.rp = list(d = K1[[1]],u = K1[[2]][,1:rank])
  d <- (K.rp$d[1:rank])^2 # approximated eigenvalues

  u <- K.rp$u[,1:rank]

  signdiag = sign(diag(t(u) %*% U1)) # find the sign change
  signdiag = as.logical(1 - signdiag)
  u[,signdiag] = -u[,signdiag]

  U <- u
  z <- xbeta + rep(U %*% (sqrt(d)*current_delta), n_t) + rep(current_v, each = n_s) # U is from random projection
  foo2 <- crossprod(current_delta, current_delta)
  likelihood <- (
    sum(dens_fun_log(O, mean = z)) - 0.5*1/current_sigma2_sp * foo2 # likelihood
    # + log(phi - phi_sp_a) + log(phi_sp_b - phi) # jacobian
  )
  return(list(likelihood = likelihood, d = d, twKinvw = foo2, U = U, u = u))
}
# (-s2_a - 1 - rank/2)*log(current_sigma2_sp) - (s2_b + 0.5*twKinvw)/current_sigma2_sp
sigma2_log_full_conditional <- function(sigma2, twKinvw,
                                        s2_a, s2_b,
                                        rank) {
  (-s2_a - 1 - rank/2)*log(sigma2) - (s2_b + 0.5*twKinvw)/sigma2
}

alpha_log_full_conditional <- function(alpha, xbeta,  V, l,
                                       O,
                                       current_w,
                                       n_s, n_t,
                                       current_sigma2_tm,
                                       dens_fun_log){ # delta is rank-m


  v = V %*% (sqrt(l) * alpha)

  z <- xbeta + rep(current_w, times = n_t) + rep(v, each = n_s)
  foo2 <- crossprod(alpha, alpha) # d = D^2 from random projection
  lf <- sum(dens_fun_log(O, mean = z)) - 1/(2*current_sigma2_tm) * foo2
  return(list(lr = lf, twKinvw = foo2, v = v))
}

phi_tm_log_full_conditional <- function(phi, dist_time, xbeta, current_alpha, V1, # data and params
                                        O, # observations
                                        current_w,
                                        current_sigma2_tm, # priors
                                        nu, n_s, n_t, cores, rank, # control params
                                        dens_fun_log){ # density function

  J1 = rp(phi, dist_time, n_t, rank, nu, as.integer(cores), cov_fun = 1) # C++ function for approximating eigenvectors
  J.rp = list(l = J1[[1]], v = J1[[2]][,1:rank])
  l <- (J.rp$l[1:rank])^2 # approximated eigenvalues

  v <- J.rp$v[,1:rank]

  signdiag = sign(diag(t(v) %*% V1)) # find the sign change
  signdiag = as.logical(1 - signdiag)
  v[,signdiag] = -v[,signdiag]

  V <- v # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
  z <- xbeta + rep(current_w, times = n_t) + rep(V %*% (sqrt(l)*current_alpha), each = n_s) # U is from random projection
  foo2 <- crossprod(current_alpha, current_alpha)
  likelihood <- (
    sum(dens_fun_log(O, mean = z)) - 0.5*1/current_sigma2_tm * foo2 # likelihood
    # + log(phi - phi_tm_a) + log(phi_tm_b - phi) # jacobian
  )
  return(list(likelihood = likelihood, l = l, twKinvw = foo2, V = V, v = v))
}




