library(withr)


# simulate data
set.seed(451)
n <- 100
nu <- 2.5
sigma2 <- 1
phi <- 0.2
tau2 <- 0.1

x0 <- rep(1, n)
x1 <- runif(n)
x2 <- runif(n)
x3 <- rnorm(n)

X <- cbind(x0, x1, x2, x3)

beta <- c(0,1,1,-2)


h <- as.matrix(dist(data.frame(x1 = x1, x2 = x2)))
matern_covariance <- sigma2 * (1 + sqrt(5)*h/phi + 5*h^2 / (3*phi^2)) * exp(-1*sqrt(5)*h/phi)

W <- rmvn(n = 1, mu = rep(0, n), V = matern_covariance)

y_mean <- as.numeric(X %*% beta + W)

y_pois <- rpois(n, exp(y_mean))
y_norm <- rnorm(n, y_mean, sqrt(tau2))
y_binom <- rbinom(n, 1, exp(y_mean)/(1 + exp(y_mean)))

sim_data <- data.frame(x1 = x1, x2 = x2, x3 = x3, W = W,
                       y_mean = y_mean,
                       y_pois = y_pois,
                       y_norm = y_norm,
                       y_binom = y_binom)

test_that("test poisson model runs", {
  # run poisson model

  poisson_model <- rrp_glm(
        fixed = y_pois ~ x1 + x2 + x3,
        spatial =  ~ x1 + x2,
        data = sim_data,
        family = poisson(),
        covfn = covfndef(nu),
        iter = 1000,
        chains = 2,
        cores = 1,
        param_start = list("beta" = rnorm(4, 5, sd = 0.5),
                           "s2" = 2,
                           "phi" = 0.5),
        priors = list("beta.normal" = 100, # variance of beta prior
                      "s2.IG" = c(2,2), # inverse gamma params
                      "phi.Unif" = c(0.01, 1.5)), # uniform prior on phi
        tuning = list("beta" = runif(4, 0.005, 0.01),
                      "s2" = 0.1,
                      "phi" = 0.01,
                      "w" = 0.1),
        nu = nu,
        rank = 10,
        mul = 2)

  expect_type(poisson_model, "list")
})

test_that("test linear model runs", {
  linear_model <- rrp_glm(
        fixed = y_norm ~ x1 + x2 + x3,
        spatial = ~ x1 + x2,
        data = sim_data,
        family = gaussian(),
        covfn = covfndef(nu),
        iter = 1000,
        chains = 2,
        cores = 1,
        param_start = list("beta" = rnorm(4, 5, sd = 0.5),
                           "s2" = 2,
                           "phi" = 0.5),
        priors = list("beta.normal" = 100, # variance of beta prior
                      "s2.IG" = c(2,2), # inverse gamma params
                      "phi.Unif" = c(0.01, 1)), # uniform prior on phi
        tuning = list("beta" = rnorm(4, 1, 0.2),
                      "s2" = 0.1,
                      "phi" = 0.01,
                      "w" = 0.1),
        nu = nu,
        rank = 10,
        mul = 2)

  expect_type(linear_model, "list")
})

test_that("test logistic model runs", {
  expect_type(rrp_glm(
        fixed = y_binom ~ x1 + x2 + x3,
        spatial = ~ x1 + x2,
        data = sim_data,
        family = binomial(),
        covfn = covfndef(nu),
        iter = 1000,
        chains = 2,
        cores = 1,
        param_start = list("beta" = rnorm(4, 5, sd = 0.5),
                           "s2" = 2,
                           "phi" = 0.5),
        priors = list("beta.normal" = c(100), # variance of beta prior
                      "s2.IG" = c(2,2), # inverse gamma params
                      "phi.Unif" = c(0.01, 1.5)), # uniform prior on phi
        tuning = list("beta" = runif(4, 0.05, 0.015),
                      "s2" = 0.1,
                      "phi" = 0.01,
                      "w" = 0.1),
        nu = nu,
        rank = 10,
        mul = 2), "list")
})

test_that("poisson fit matches spatialPoisson fit", {
  starting <- list("beta" = rnorm(4, 5, sd = 0.5),
                   "s2" = 2,
                   "phi" = 0.5)

  priors <- list("beta.normal" = 100, # variance of beta prior
                 "s2.IG" = c(2,2), # inverse gamma params
                 "phi.Unif" = c(0.01, 1))

  tuning <- list("beta" = rnorm(4, 1, 0.2),
                 "s2" = 0.1,
                 "phi" = 0.01,
                 "w" = 0.1)

  rank = 10

  set.seed(451)
  my_fit <- rrp_glm(
      fixed = y_pois ~ x1 + x2 + x3,
      spatial =  ~ x1 + x2,
      data = sim_data,
      family = poisson(),
      covfn = covfndef(nu),
      iter = 1000,
      chains = 2,
      cores = 1,
      param_start = starting,
      priors = priors,
      tuning = tuning,
      nu = nu,
      rank = rank,
      mul = 2)

  design_matrix <- model.matrix(~ x1 + x2 + x3,
                              data = sim_data)

  set.seed(451)
  orig_fit <- poi_gp_arrpfit(
        sim_data$y_pois,
        coords = cbind(x1, x2),
        X = design_matrix,
        mul = 2,
        covfn = covfndef(nu),
        adapt = c("batchlength" = 1000,
                   "n.batch" = 2),
        nu = nu,
        core = 1,
        starting = starting,
        tuning = tuning,
        priors = priors,
        rank = 10)

  expect_equal(my_fit$p.params, orig_fit$p.params)
})

test_that("beta full conditional correct", {

  dens_fun_log <- function(x, mean) {dpois(x, lambda = poisson()$linkinv(mean), log = TRUE)}

  expect_snapshot_value(beta_log_full_conditional(beta = beta,
                            X = X,
                            wParams = W,
                            O = y_pois,
                            p = 4,
                            beta.b = 100,
                            dens_fun_log = dens_fun_log), style = "serialize")

})


test_that("delta full conditional correct", {
dens_fun_log <- function(x, mean) {dpois(x, lambda = poisson()$linkinv(mean), log = TRUE)}

  cov_eigen <- eigen(matern_covariance)

  U <- cov_eigen$vectors[,1:10]
  d <- cov_eigen$values[1:10] ^ 2

  expect_snapshot_value(delta_log_full_conditional(delta = -4:5/10,
                            O = y_pois,
                            U = U,
                            d = d,
                            xbeta = X %*% beta,
                            sParams = 1,
                            s2indx = 1,
                            dens_fun_log = dens_fun_log), style = "serialize")
})


test_that("phi full conditional correct", {
  dens_fun_log <- function(x, mean) {dpois(x, lambda = poisson()$linkinv(mean), log = TRUE)}

  cov_eigen <- eigen(matern_covariance)

  U <- cov_eigen$vectors[,1:10]
  d <- cov_eigen$values[1:10] ^ 2

  AP = chol2inv(chol(crossprod(X,X))) %*% t(X) # projection onto column space of X
  PPERP <- diag(n) - X %*% AP # projection onto space orthogonal to column space of X


  expect_snapshot_value(phi_log_full_conditional(phi = 0.2,
                            coords = cbind(x1, x2),
                            O = y_pois,
                            n = n,
                            rk = 10,
                            nu = nu,
                            cores = 1,
                            dens_fun_log = dens_fun_log,
                            U1 = U,
                            PPERP = PPERP,
                            rank = 10,
                            xbeta = X %*% beta,
                            etaParams = -4:5/10,
                            sParams = 1,
                            s2indx = 1), style = "serialize")})
