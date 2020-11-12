rrp_glm <- function(fixed,
                    spatial,
                    data,
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


  # Spatial effect
  coords <- as.matrix(lm(spatial, data = data, method = "model.frame"))

  # TODO: change this to param_start, then figure out how to get rid of it

  starting <- param_start
  ptm <- proc.time() # starting time
  p = ncol(X) # expected rank of X
  n <- length(O) # number of observations
  rk = rank*mul # rank of approximation matrix
  AP = chol2inv(chol(crossprod(X,X))) %*% t(X) # projection onto column space of X
  PPERP <- diag(n) - X %*% AP # projection onto space orthogonal to column space of X


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
  samples_arrp <- matrix(NA,ncol = p,nrow = niter)
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
    sParams[phi_index] <- phi.a + (phi.b - phi.a)/(1 + exp(-runif(1, min = -2, max = 2)))



    est_start  <- Sys.time()
    # this is where the first call to the cpp code occurs. It is for the random projections part of the algorithm.
    K = rp(
              sParams[phi_index], # single number, phi in starting param list
              coords, #literally x,y coords of obs
              as.integer(n), # number if interations (iter)
              as.integer(rk), # rank times mul (what is mul?)
              nu, #nu as before
              as.integer(cores)# number of cores
              ) # C++ function for approximating eigenvectors
    est_time  <- Sys.time() - est_start # calculates how long one random projection takes
    message("Estimated time (hrs):",round((chains - k + 1)*iter*2*est_time/3600, 3) ,"\n") # prints out estimate time in hours

    K.rp = list(d = K[[1]],u = K[[2]][,1:rank]) # d is sing values, u is left sing vecs, only take first 1:rank
    d <- (K.rp$d[1:rank])^2 # approximated eigenvalues
    U1 <- U <- u <- K.rp$u[,1:rank] # approximated eigenvectors
    U <- mmC(PPERP,U,as.integer(n),as.integer(rank),as.integer(cores)) # compute PPERP%*%u restrict random effect to be orthogonal to fix effect

    etaParams <- runif(rank, min = -2, max = 2)
    wParams <- U %*% (sqrt(d)*etaParams)


    xbeta <- X %*% current_beta



    for (i in 1:iter) {

      # block update beta
      betastar <- rnorm(p, current_beta, sd = exp(sTunings[beta_index])) # draw for proposed beta params

      beta.lfcur <- beta_log_full_conditional(current_beta, wParams = wParams, X = X,
                                              O = O,
                                              beta.b = beta.b,
                                              p = p,
                                              dens_fun_log = dens_fun_log)

      beta.lfcand <- beta_log_full_conditional(betastar, wParams = wParams, X = X,
                                              O = O,
                                              beta.b = beta.b,
                                              p = p,
                                              dens_fun_log = dens_fun_log)
      lr <- beta.lfcand - beta.lfcur # log likelihood ratio log(proposed_likelihood/current_likelihood)


      # this is metropolis hastings step
      if (log(runif(1)) < lr) { # if likelihood ratio is greater than runif(1), accept
        current_beta <- betastar # update params with new guess
        sAccepts[beta_index] <- sAccepts[beta_index] + 1 # mark we accepted
        xbeta <- X %*% current_beta # update xbeta, which will be used in next set of estimation
      }


      # project beta guess to be orthoganal to
      AParams = current_beta - AP %*% u %*% (sqrt(d)*etaParams) # adjust the random effects to get ARRP


      # guess phi params
      phistar <-  rnorm(1, sParams[phi_index], sd = exp(sTunings[phi_index]))
      phi.lfcand <- phi.lfcur <- phi_log_full_conditional(sParams[phi_index], coords = coords, xbeta = xbeta, etaParams = etaParams, U1 = U1, PPERP = PPERP, # data and params
                                                          O = O, # observations
                                                          current_sigma2 = current_sigma2, # priors
                                                          nu = nu, n = n, rk = rk, cores = cores, rank = rank, # control params
                                                          dens_fun_log = dens_fun_log)



      # checks if guess in bounds of Unif(a, b) prior
      if (phistar < phi.b & phistar > phi.a) {
        phi.lfcand <- phi_log_full_conditional(phistar, coords = coords, xbeta = xbeta, etaParams = etaParams, U1 = U1, PPERP = PPERP, # data and params
                                               O = O, # observations
                                               current_sigma2 = current_sigma2, # priors
                                               nu = nu, n = n, rk = rk, cores = cores, rank = rank, # control params
                                               dens_fun_log = dens_fun_log)
      } else {
        phi.lfcand$lr <- -Inf
      }
      lr <- phi.lfcand$lr - phi.lfcur$lr

      if (log(runif(1)) < lr) {
        sParams[phi_index] <- phistar
        sAccepts[phi_index] <- sAccepts[phi_index] + 1
        phi.lfcur <- phi.lfcand
        d <- phi.lfcur$d # approximated eigenvalues^2
        U <- phi.lfcur$U # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
        u <- phi.lfcur$u # approximated eigenvectors
      }

      twKinvw <- phi.lfcur$twKinvw # etaParams cross product for some reason

      # update s2
      s2star <- rnorm(1, current_sigma2, sd = exp(sTunings[sigma2_index]))
      if (s2star > 0) {
        s2.lfcur <- sigma2_log_full_conditional(sigma2 = current_sigma2, twKinvw = twKinvw,
                                                s2.a = s2.a, s2.b = s2.b,
                                                rank = rank)

        s2.lfcand <-  sigma2_log_full_conditional(sigma2 = s2star, twKinvw = twKinvw,
                                                  s2.a = s2.a, s2.b = s2.b,
                                                  rank = rank)
        lr <- s2.lfcand - s2.lfcur
      } else {
        lr <- -Inf
      }


      if (log(runif(1)) < lr) {
        current_sigma2 <- s2star
        sAccepts[sigma2_index] <- sAccepts[sigma2_index] + 1
      }

      # update random effects using multivariate random walk with spherical normal proposal
      deltastar <- rnorm(rank, etaParams, sd = exp(wTunings))

      delta.lfcand <- delta_log_full_conditional(delta = deltastar, xbeta = xbeta,  U = U, d = d,
                                                 O = O,
                                                 current_sigma2 = current_sigma2,
                                                 dens_fun_log = dens_fun_log)

      delta.lfcur <- delta_log_full_conditional(delta = etaParams, xbeta = xbeta,  U = U, d = d,
                                                O = O,
                                                current_sigma2 = current_sigma2,
                                                dens_fun_log = dens_fun_log)
      lr <- delta.lfcand$lr - delta.lfcur$lr

      if (log(runif(1)) < lr) {
        etaParams <- deltastar
        wAccepts <- wAccepts + 1
        delta.lfcur <- delta.lfcand
        wParams <- delta.lfcur$w
      }

      samples_arrp[(k - 1)*iter + i,] <- AParams
      samples_s[(k - 1)*iter + i,] <- c(current_beta, current_sigma2, sParams[phi_index])
      samples_w[(k - 1)*iter + i,] <- wParams
      samples_eta[(k - 1)*iter + i,] <- etaParams

      samples[i, k, 1:nParams] <- c(current_beta, current_sigma2, sParams[phi_index])
      samples[i, k, (nParams + 1):(nParams + rank)] <- etaParams
      samples[i, k, (nParams + 1 + rank):(nParams + rank + n)] <- wParams
    }

    output_progress(samples_s = samples_s, sAccepts = sAccepts,
                    iter = iter, k = k)

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
  #                       log_likelihood = ,# log likelihood array for loo
  #                       model = "rrp/arrp",
  #                       accept = list(beta = ,
  #                                     sigma2 = ,
  #                                     phi = ,
  #                                     delta = ,
  #                                     w = )),
  #                  class = "rsrHonors_rrp_sglm"))



  return(list(run.time = runtime,
              p.params = samples_s ,
              arrp.params = samples_arrp,
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

output_progress <- function(samples_s, sAccepts,
                            iter, k) {
  message("Batch ",k,"\n")
  message("-------------------------------------------------------\n")
  message("----------------parameter estimates--------------------\n")
  message(paste0(round(bmmat(samples_s[((k - 1)*iter + 1):(k*iter),])[,1], 3), collapse = "   "),"\n")
  message("-------------------------------------------------------\n")
  message("----------------acceptance rate------------------------\n")
  message(paste0(round(sAccepts/(iter),3), collapse = "   "),"\n")
  message("-------------------------------------------------------\n")
}

# covfndef is in util.r; include all util functions
# define covariance function
covfndef <- function(nu){
  # exponential
  if (nu == 1/2) covfn <- function(dist,phi) {
    K = exp(-1/phi*dist)
    return(K)
  }
  # matern 1.5
  if (nu == 1.5) covfn <- function(dist,phi) {
    K = (1 + sqrt(3)/phi*dist)*exp(-sqrt(3)/phi*dist)
    return(K)
  }
  # matern 2.5
  if (nu == 2.5 ) covfn <- function(dist,phi) {
    K = (1 + sqrt(5)/phi*dist + 5/(3*phi^2)*dist^2)*exp(-sqrt(5)/phi*dist)
    return(K)
  }
  # square exponential
  if (nu == 10) covfn <- function(dist,phi) {
    K = exp(-1/(2*phi^2)*dist^2)
    return(K)
  }
  return(covfn)
}

# generate multivariate normal
rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if (any(is.na(match(dim(V), p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

# compute batchmeans, function from R packages batchmeans
bm <- function(x, size = "sqroot", warn = FALSE)
{
  n = length(x)
  if (n < 1000) {
    if (warn)
      warning("too few samples (less than 1,000)")
    if (n < 10)
      return(NA)
  }
  if (size == "sqroot") {
    b = floor(sqrt(n))
    a = floor(n/b)
  }
  else if (size == "cuberoot") {
    b = floor(n^(1/3))
    a = floor(n/b)
  }
  else {
    if (!is.numeric(size) || size <= 1 || size == Inf)
      stop("'size' must be a finite numeric quantity larger than 1.")
    b = floor(size)
    a = floor(n/b)
  }
  y = sapply(1:a, function(k) return(mean(x[((k - 1) * b +
                                               1):(k * b)])))
  mu.hat = mean(y)
  var.hat = b * sum((y - mu.hat)^2)/(a - 1)
  se = sqrt(var.hat/n)
  list(est = mu.hat, se = se)
}

# compute batchmeans, function from R packages batchmeans
bmmat <- function(x)
{
  if (!is.matrix(x) && !is.data.frame(x))
    stop("'x' must be a matrix or data frame.")
  num = ncol(x)
  bmvals = matrix(NA, num, 2)
  rownames(bmvals) = colnames(x)
  bmres = apply(x, 2, bm)
  for (i in 1:num) bmvals[i, ] = c(bmres[[i]]$est, bmres[[i]]$se)
  bmvals
}

beta_log_full_conditional <- function(beta, wParams, X,
                                      O,
                                      beta.b,
                                      p,
                                      dens_fun_log){

    z <- X %*% beta + wParams # if rsr, wParams = L%*%eta
    lf <- sum(dens_fun_log(O, mean = z)) - crossprod(beta)/(p*beta.b)
    return(lf)

}

delta_log_full_conditional <- function(delta, xbeta,  U, d,
                                       O,
                                       current_sigma2,
                                       dens_fun_log){ # delta is rank-m
    w = U %*% (sqrt(d)*delta)
    z <- xbeta + w
    foo2 <- crossprod(delta,delta) # d = D^2 from random projection
    lf <- sum(dens_fun_log(O, mean = z)) - 1/(2*current_sigma2) * foo2
    return(list(lr = lf, twKinvw = foo2, w = w))
}

phi_log_full_conditional <- function(phi, coords, xbeta, etaParams, U1, PPERP, # data and params
                                     O, # observations
                                     current_sigma2, # priors
                                     nu, n, rk, cores, rank, # control params
                                     dens_fun_log){ # density function
  K1 = rp(phi,coords,as.integer(n) ,as.integer(rk),nu,as.integer(cores)) # C++ function for approximating eigenvectors
  K.rp = list(d = K1[[1]],u = K1[[2]][,1:rank])
  d <- (K.rp$d[1:rank])^2 # approximated eigenvalues

  u <- K.rp$u[,1:rank]

  signdiag = sign(diag(t(u) %*% U1)) # find the sign change
  signdiag = as.logical(1 - signdiag)
  u[,signdiag] = -u[,signdiag]

  U <- mmC(PPERP,u,as.integer(n),as.integer(rank),as.integer(cores)) # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
  z <- xbeta + U %*% (sqrt(d)*etaParams) # U is from random projection
  foo2 <- crossprod(etaParams,etaParams)
  lr <- (
    sum(dens_fun_log(O, mean = z)) - 0.5*1/current_sigma2 * foo2 # likelihood
    # + log(phi - phi.a) + log(phi.b - phi) # jacobian
  )
  return(list(lr = lr,d = d,twKinvw = foo2, U = U, u = u))
}
# (-s2.a - 1 - rank/2)*log(current_sigma2) - (s2.b + 0.5*twKinvw)/current_sigma2
sigma2_log_full_conditional <- function(sigma2, twKinvw,
                                        s2.a, s2.b,
                                        rank) {
  (-s2.a - 1 - rank/2)*log(sigma2) - (s2.b + 0.5*twKinvw)/sigma2
}

# R wrappers to call cpp function for random projection
rpcpp <- function(r, coords, phi, nu){ # r is demension selected, K,n is loaded from global environment
  n = nrow(coords)
  rp(phi, coords, as.integer(n), as.integer(r), nu, as.integer(1))
}
