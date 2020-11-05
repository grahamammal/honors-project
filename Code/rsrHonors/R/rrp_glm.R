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


  call <- match.call()


  # Figure out how to interpret the formula fixed fixed
  model_frame <- lm(fixed, data = data, method = "model.frame")
  X <- model.matrix(fixed, data = model_frame)
  y <- model.response(model_frame)

  # family things
  if (family != "poisson" || family$family != "poisson") {
    stop("I haven't implmented this family lol")
  }

  # Spatial effect

  spatial_frame <- lm(spatial, data = data, method = "model.frame")

  ptm <- proc.time() # starting time
  p = ncol(X) # expected rank of X
  n <- length(O) # number of observations
  rk = rank*mul # rank of approximation matrix
  AP = chol2inv(chol(crossprod(X,X))) %*% t(X) # projection onto column space of X
  PPERP <- diag(n) - X %*% AP # projection onto space orthogonal to column space of X


  #######################################
  # Prepare MCMC iteration values
  #######################################
  batchlength = adapt[["batchlength"]]
  n.batch = adapt[["n.batch"]]
  niter <- batchlength*n.batch
  cat("MCMC chain size:", batchlength,"\n")

  #######################################
  # Define prior parameters
  #######################################
  beta.b <- priors[["beta.normal"]]
  s2.a   <- priors[["s2.IG"]][1]
  s2.b   <- priors[["s2.IG"]][2]
  phi.a  <- priors[["phi.Unif"]][1]
  phi.b  <- priors[["phi.Unif"]][2]



  # might be doing gibbs sampler.
  # Full conditional sounds like gibbs sampler
  # Might be doing component wise metropolis hastings.
  #

  # log full conditional of beta
  # NOTE: should it be (p*beta.b) instead of (2*beta.b)???
  #
  # Estimates fixed effects parameters
  # lf: look at notes
  beta.lf <- function(beta){
    z <- X %*% beta + wParams # if rsr, wParams = L%*%eta
    lf <- sum(dpois(O, exp(z), log = TRUE)) - crossprod(beta)/(p*beta.b)
    return(lf)
  }
  # log full conditional of delta

  delta.lf <- function(delta){ # delta is rank-m
    w = U %*% (sqrt(d)*delta)
    z <- xbeta + w
    foo2 <- crossprod(delta,delta) # d = D^2 from random projection
    lf <- sum(dpois(O, exp(z), log = TRUE)) - 1/(2*sParams[s2indx]) * foo2
    return(list(lr = lf, twKinvw = foo2, w = w))
  }

  phi.lf <- function(phi){
    K1 = rp(phi,coords,as.integer(n) ,as.integer(rk),nu,as.integer(core)) # C++ function for approximating eigenvectors
    K.rp = list(d = K1[[1]],u = K1[[2]][,1:rank])
    d <- (K.rp$d[1:rank])^2 # approximated eigenvalues

    u <- K.rp$u[,1:rank]

    signdiag = sign(diag(t(u) %*% U1)) # find the sign change
    signdiag = as.logical(1 - signdiag)
    u[,signdiag] = -u[,signdiag]

    U <- mmC(PPERP,u,as.integer(n),as.integer(rank),as.integer(core)) # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
    z <- xbeta + U %*% (sqrt(d)*etaParams) # U is from random projection
    foo2 <- crossprod(etaParams,etaParams)
    lr <- (
      sum(dpois(O,exp(z),log = TRUE)) - 0.5*1/sParams[s2indx] * foo2 # likelihood
      # + log(phi - phi.a) + log(phi.b - phi) # jacobian
    )
    return(list(lr = lr,d = d,twKinvw = foo2, U = U, u = u))
  }


  ###################################
  # Orig: define index for parameters
  # Must mean that parameters are stored in a single vector
  # This vector is called sParams
  ####################################
  nParams <- p + 2;
  betaindx <- 1:p;
  s2indx <- p + 1;
  phiindx <- s2indx + 1



  #########################################
  # Orig: initialize some matrix for storage
  # Prepares matrices to hold estimates of parameters
  #########################################
  samples_eta <- matrix(NA,ncol = rank,nrow = niter) # store the r.e samples
  samples_s <- matrix(NA,ncol = nParams,nrow = niter) # store the parameter samples
  samples_w <- matrix(NA,ncol = n,nrow = niter) # store the r.e samples
  samples_arrp <- matrix(NA,ncol = p,nrow = niter)
  sTunings <- sParams <- matrix(NA,nrow = nParams) # update the current parameters for each iteration

  #########################################
  # Orig: feed in the initial values and tuning params
  # Set starting parameters
  #########################################
  sParams[betaindx] <- starting[["beta"]];  sParams[s2indx] <- starting[["s2"]];  sParams[phiindx] <- starting[["phi"]]
  sTunings[betaindx] <- tuning[["beta"]];  sTunings[s2indx] <- tuning[["s2"]];  sTunings[phiindx] <- tuning[["phi"]]
  wTunings <- rep(tuning[["w"]],rank)
  wTunings <- log(wTunings)
  sTunings <- log(sTunings)


  # Stores info on acceptance rates for MCMC
  # for acceptance rate
  accept_s <- matrix(NA, ncol = nParams, nrow = n.batch)
  accept_w <- matrix(NA, ncol = rank, nrow = n.batch)


  est.time  <- proc.time()
  # this is where the first call to the cpp code occurs. It is for the random projections part of the algorithm.
  K = rp(
            sParams[phiindx], # single number, phi in starting param list
            coords, #literally x,y coords of obs
            as.integer(n), # number if interations (batchlength)
            as.integer(rk), # rank times mul (what is mul?)
            nu, #nu as before
            as.integer(core)# number of cores
            ) # C++ function for approximating eigenvectors
  est.time  <- proc.time() - est.time # calculates how long one random projection takes
  cat("Estimated time (hrs):",niter*2*est.time[3]/3600 ,"\n") # prints out estimate time in hours

  K.rp = list(d = K[[1]],u = K[[2]][,1:rank]) # d is sing values, u is left sing vecs, only take first 1:rank
  d <- (K.rp$d[1:rank])^2 # approximated eigenvalues
  U1 <- U <- u <- K.rp$u[,1:rank] # approximated eigenvectors
  U <- mmC(PPERP,U,as.integer(n),as.integer(rank),as.integer(core)) # compute PPERP%*%u restrict random effect to be orthogonal to fix effect

  if (is.null(starting[["w"]])) {
    etaParams <- rep(0,rank)
    wParams <- U %*% (sqrt(d)*etaParams)
  }else{
    wParams <- starting[["w"]]
    etaParams <- 1/sqrt(d)*(t(U) %*% wParams)
  }

  xbeta <- X %*% sParams[betaindx]

  ##### start MCMC loop #####
  for (k in 1:n.batch) {
    sAccepts <- matrix(0,nrow = nParams)
    wAccepts <- matrix(0,nrow = rank)


    for (i in 1:batchlength) {

      # block update beta
      betastar <- rnorm(p, sParams[betaindx], sd = exp(sTunings[betaindx])) # draw for proposed beta params
      beta.lfcur <- beta.lf(sParams[betaindx]) # log likelihood of current beta params
      beta.lfcand <- beta.lf(c(betastar)) # log likelihood of proposed beta params
      lr <- beta.lfcand - beta.lfcur # log likelihood ratio log(proposed_likelihood/current_likelihood)


      # this is metropolis hastings step
      if (log(runif(1)) < lr) { # if likelihood ratio is greater than runif(1), accept
        sParams[betaindx] <- betastar # update params with new guess
        sAccepts[betaindx] <- sAccepts[betaindx] + 1 # mark we accepted
        xbeta <- X %*% sParams[betaindx] # update xbeta, which will be used in next set of estimation
      }

      # project beta guess to be orthoganal to
      AParams = sParams[betaindx] - AP %*% u %*% (sqrt(d)*etaParams) # adjust the random effects to get ARRP

      # update phi
      # phistar <- logit.inv(rnorm(1,logit(sParams[phiindx],phi.a,phi.b),
      # sd=exp(sTunings[phiindx])),phi.a,phi.b)

      # phi.lfcur <- phi.lf(sParams[phiindx])
      # phi.lfcand <- phi.lf(phistar)

      phistar <-  rnorm(1, sParams[phiindx], sd = exp(sTunings[phiindx]))
      phi.lfcand <- phi.lfcur <- phi.lf(sParams[phiindx])

      if (phistar < phi.b & phistar > phi.a) {
        phi.lfcand <- phi.lf(phistar)
      } else {
        phi.lfcand$lr <- -Inf
      }
      lr <- phi.lfcand$lr - phi.lfcur$lr

      if (log(runif(1)) < lr) {
        sParams[phiindx] <- phistar
        sAccepts[phiindx] <- sAccepts[phiindx] + 1
        phi.lfcur <- phi.lfcand
        d <- phi.lfcur$d # approximated eigenvalues^2
        U <- phi.lfcur$U # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
        u <- phi.lfcur$u # approximated eigenvectors
      }

      twKinvw <- phi.lfcur$twKinvw

      # update s2
      s2star <- rnorm(1, sParams[s2indx], sd = exp(sTunings[s2indx]))
      if (s2star > 0) {
        s2.lfcur <-  (-s2.a - 1 - rank/2)*log(sParams[s2indx]) - (s2.b + 0.5*twKinvw)/sParams[s2indx]
        s2.lfcand <-  (-s2.a - 1 - rank/2)*log(s2star) - (s2.b + 0.5*twKinvw)/s2star
        lr <- s2.lfcand - s2.lfcur
      } else {
        lr <- -Inf
      }


      if (log(runif(1)) < lr) {
        sParams[s2indx] <- s2star
        sAccepts[s2indx] <- sAccepts[s2indx] + 1
      }

      #       # update re one-by-one # either this or block update will work
      #
      #       deltastar<-etaParams
      #       delta.lfcur<- delta.lf(etaParams)
      #       for(j in 1:rank){
      #         deltastar[j] <- rnorm(1,etaParams[j],exp(wTunings[j]))
      #         delta.lfcand <- delta.lf(deltastar)
      #         lr <- delta.lfcand$lr - delta.lfcur$lr
      #
      #         if(log(runif(1)) < lr) {
      #           etaParams[j] <- deltastar[j]
      #           wAccepts[j] <- wAccepts[j] +1
      #           delta.lfcur <- delta.lfcand
      #         } else{
      #           deltastar[j] <- etaParams[j]
      #         }
      #       }.Call("rp", sParams[phiindx]
      #       wParams <- delta.lfcur$w

      # update random effects using multivariate random walk with spherical normal proposal
      deltastar <- rnorm(rank, etaParams, sd = exp(wTunings))
      delta.lfcand <- delta.lf(deltastar)
      delta.lfcur <- delta.lf(etaParams)
      lr <- delta.lfcand$lr - delta.lfcur$lr

      if (log(runif(1)) < lr) {
        etaParams <- deltastar
        wAccepts <- wAccepts + 1
        delta.lfcur <- delta.lfcand
        wParams <- delta.lfcur$w
      }

      samples_arrp[(k - 1)*batchlength + i,] <- AParams
      samples_s[(k - 1)*batchlength + i,] <- sParams
      samples_w[(k - 1)*batchlength + i,] <- wParams
      samples_eta[(k - 1)*batchlength + i,] <- etaParams
    }

    cat(" Batch ",k,"\n")
    cat("-------------------------------------------------------\n")
    cat("----------------est------------------------------------\n")
    cat(bmmat(samples_s[((k - 1)*batchlength + 1):(k*batchlength),])[,1],"\n")
    cat("-------------------------------------------------------\n")

    cat("acceptance rate:", round(sAccepts/(batchlength),4),"\n")
    cat("-------------------------------------------------------\n")

    accept_s[k,] <- sAccepts/(batchlength)
    accept_w[k,] <- wAccepts/(batchlength)

  }
  accept_s <- apply(accept_s,2,mean)
  accept_w <- apply(accept_w,2,mean)
  runtime = proc.time() - ptm
  cat("runtime", "\n")
  print(runtime)

  return(list(run.time = runtime,
              p.params = samples_s ,
              arrp.params = samples_arrp,
              model = "rrp/arrp",
              rank = rank,
              eta.params = samples_eta,
              w.params = samples_w,
              accepts = accept_s,
              acceptw = accept_w))
}

# --------------------------------------------------------
# functions from source code:
# --------------------------------------------------------

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

#logit transform
logit <- function(theta, a, b){log((theta - a)/(b - theta))}
logit.inv <- function(z, a, b){b - (b - a)/(1 + exp(z))}

# generate inverse gamma random variable
rinvgamma <- function(n, shape, scale = 1)
{
  return(1/rgamma(n = n, shape = shape, rate = scale))
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
  colnames(bmvals) = c("est", "se")
  rownames(bmvals) = colnames(x)
  bmres = apply(x, 2, bm)
  for (i in 1:num) bmvals[i, ] = c(bmres[[i]]$est, bmres[[i]]$se)
  bmvals
}

beta_log_full_conditional <- function(beta){
    z <- X %*% beta + wParams # if rsr, wParams = L%*%eta
    lf <- sum(dpois(O, exp(z), log = TRUE)) - crossprod(beta)/(p*beta.b)
    return(lf)
}

delta_log_full_conditional <- function(delta){ # delta is rank-m
    w = U %*% (sqrt(d)*delta)
    z <- xbeta + w
    foo2 <- crossprod(delta,delta) # d = D^2 from random projection
    lf <- sum(dpois(O, exp(z), log = TRUE)) - 1/(2*sParams[s2indx]) * foo2
    return(list(lr = lf, twKinvw = foo2, w = w))
}

phi_log_full_conditional <- function(phi){
  K1 = rp(phi,coords,as.integer(n) ,as.integer(rk),nu,as.integer(core)) # C++ function for approximating eigenvectors
  K.rp = list(d = K1[[1]],u = K1[[2]][,1:rank])
  d <- (K.rp$d[1:rank])^2 # approximated eigenvalues

  u <- K.rp$u[,1:rank]

  signdiag = sign(diag(t(u) %*% U1)) # find the sign change
  signdiag = as.logical(1 - signdiag)
  u[,signdiag] = -u[,signdiag]

  U <- mmC(PPERP,u,as.integer(n),as.integer(rank),as.integer(core)) # gives PPERP%*%u restrict random effect to be orthogonal to fix effect
  z <- xbeta + U %*% (sqrt(d)*etaParams) # U is from random projection
  foo2 <- crossprod(etaParams,etaParams)
  lr <- (
    sum(dpois(O,exp(z),log = TRUE)) - 0.5*1/sParams[s2indx] * foo2 # likelihood
    # + log(phi - phi.a) + log(phi.b - phi) # jacobian
  )
  return(list(lr = lr,d = d,twKinvw = foo2, U = U, u = u))
}


# R wrappers to call cpp function for random projection
rpcpp <- function(r, coords, phi, nu){ # r is demension selected, K,n is loaded from global environment
  n = nrow(coords)
  rp(phi, coords, as.integer(n), as.integer(r), nu, as.integer(1))
}
