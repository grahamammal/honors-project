library(withr)


# simulate data
set.seed(451)
n_s <- 100
n_t <- 10


nu <- 2.5
sigma2_sp <- 1
sigma2_tm <- 1
phi_sp <- 0.5
phi_tm <- 1
tau2 <- 0.1

x0 <- rep(1, n_s * n_t)
x1 <- rep(runif(n_s), times = n_t)
x2 <- rep(runif(n_s), times = n_t)
x3 <- rnorm(n_s*n_t)

time <- 1:n_t / 10

locations <- data.frame(x1 = x1[1:n_s], x2 = x2[1:n_s])


X <- cbind(x0, x1, x2, x3)
beta <- c(0,1,1,-2)


d <- as.matrix(dist(data.frame(x1 = x1[1:n_s],
                               x2 = x2[1:n_s])))
matern_covariance <- sigma2_sp * (1 + sqrt(5)*d/phi_sp + 5*d^2 / (3*phi_sp^2)) * exp(-1*sqrt(5)*d/phi_sp)

h <- as.matrix(dist(time))
exp_cov <- sigma2_tm * exp(-h/phi_tm)


W <- rmvn(n = 1, mu = rep(0, n_s), V = matern_covariance)
V <- rmvn(n = 1, mu = rep(0, n_t), V = exp_cov)


y_mean <- as.numeric(X %*% beta + rep(W, times = n_t) + rep(V, each = n_s))

y_pois <- rpois(n_s * n_t, exp(y_mean))
y_norm <- rnorm(n_s * n_t, y_mean, sqrt(tau2))
y_binom <- rbinom(n_s * n_t, 1, exp(y_mean)/(1 + exp(y_mean)))

sim_data <- data.frame(x1 = x1, x2 = x2, x3 = x3,
                       W = rep(W, times = n_t),
                       V = rep(V, each = n_s),
                       y_mean = y_mean,
                       y_pois = y_pois,
                       y_norm = y_norm,
                       y_binom = y_binom)

test_that("test poisson stglmm runs", {
  # run poisson model

  poisson_model <- stglmm(
    fixed = y_pois ~ x1 + x2 + x3,
    data = sim_data,
    locations = locations,
    times = time,
    family = poisson(),
    covfn = covfndef(nu),
    iter = 1000,
    chains = 2,
    cores = 1,
    param_start = list("beta" = rnorm(4, 5, sd = 0.5),
                       "s2" = 2,
                       "phi" = 0.5),
    priors = list("beta.normal" = 100, # variance of beta prior
                  "s2_sp_IG" = c(2,2), # inverse gamma params
                  "s2_tm_IG" = c(2,2),
                  "phi_sp_unif" = c(0.01, 1.5),
                  "phi_tm_unif" = c(0.01, 1.5)), # uniform prior on phi
    tuning = list("beta" = rnorm(4, 1, 0.2),
                  "s2_sp" = 0.1,
                  "phi_sp" = 0.01,
                  "s2_tm" = 0.1,
                  "phi_tm" = 0.01,
                  "delta" = 0.1,
                  "alpha" = 0.1),
    nu = nu,
    rank_sp = 10,
    rank_tm = 10)

  expect_type(poisson_model, "list")
  expect_snapshot(poisson_model$p.params)
})

test_that("test linear stglmm runs", {
  linear_model <- stglmm(
    fixed = y_norm ~ x1 + x2 + x3,
    locations = locations,
    times = time,
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
                  "s2_sp_IG" = c(2,2), # inverse gamma params
                  "s2_tm_IG" = c(2,2),
                  "phi_sp_unif" = c(0.01, 1.5),
                  "phi_tm_unif" = c(0.01, 1.5)),
    tuning = list("beta" = rnorm(4, 1, 0.2),
                  "s2_sp" = 0.1,
                  "phi_sp" = 0.01,
                  "s2_tm" = 0.1,
                  "phi_tm" = 0.01,
                  "delta" = 0.1,
                  "alpha" = 0.1),
    nu = nu,
    rank_sp = 10,
    rank_tm = 10)

  expect_type(linear_model, "list")
})

test_that("test logistic stglmm runs", {
  expect_type(stglmm(
    fixed = y_binom ~ x1 + x2 + x3,
    locations = locations,
    times = time,
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
                  "s2_sp_IG" = c(2,2), # inverse gamma params
                  "s2_tm_IG" = c(2,2),
                  "phi_sp_unif" = c(0.01, 1.5),
                  "phi_tm_unif" = c(0.01, 1.5)),
    tuning = list("beta" = rnorm(4, 1, 0.2),
                  "s2_sp" = 0.1,
                  "phi_sp" = 0.01,
                  "s2_tm" = 0.1,
                  "phi_tm" = 0.01,
                  "delta" = 0.1,
                  "alpha" = 0.1),
    nu = nu,
    rank_sp = 10,
    rank_tm = 10), "list")
})

test_that("stglmm beta full conditional correct", {
  set.seed(451)
  dens_fun_log <- function(x, mean) {dpois(x, lambda = poisson()$linkinv(mean), log = TRUE)}

  expect_snapshot_value(beta_log_full_conditional_st(beta = beta,
                                                  X = X,
                                                  current_w = W,
                                                  O = y_pois,
                                                  p = 4,
                                                  beta.b = 100,
                                                  dens_fun_log = dens_fun_log,
                                                  n_t = n_t), style = "serialize")

})


test_that("stglmm delta full conditional correct", {
  set.seed(451)
  dens_fun_log <- function(x, mean) {dpois(x, lambda = poisson()$linkinv(mean), log = TRUE)}

  cov_eigen <- eigen(matern_covariance)

  U <- cov_eigen$vectors[,1:10]
  d <- cov_eigen$values[1:10] ^ 2

  expect_snapshot_value(delta_log_full_conditional_st(delta = -4:5/10,
                                                   O = y_pois,
                                                   U = U,
                                                   d = d,
                                                   n_t = n_t,
                                                   xbeta = X %*% beta,
                                                   current_sigma2_sp = 1,
                                                   dens_fun_log = dens_fun_log), style = "serialize")
})


test_that("stglmm phi_sp full conditional correct", {
  set.seed(451)
  dens_fun_log <- function(x, mean) {dpois(x, lambda = poisson()$linkinv(mean), log = TRUE)}

  cov_eigen <- eigen(matern_covariance)

  U <- cov_eigen$vectors[,1:10]
  d <- cov_eigen$values[1:10] ^ 2


  expect_snapshot_value(phi_sp_log_full_conditional(phi = 0.2,
                                                    dist_space = as.matrix(dist(locations)),
                                                 O = y_pois,
                                                 n_s = n_s,
                                                 n_t = n_t,
                                                 nu = nu,
                                                 cores = 1,
                                                 dens_fun_log = dens_fun_log,
                                                 U1 = U,
                                                 rank = 10,
                                                 xbeta = X %*% beta,
                                                 current_delta = -4:5/10,
                                                 current_sigma2_sp = 1), style = "serialize")})


test_that("stglmm alpha full conditional correct", {
  set.seed(451)
  dens_fun_log <- function(x, mean) {dpois(x, lambda = poisson()$linkinv(mean), log = TRUE)}

  cov_eigen <- eigen(exp_cov)

  U <- cov_eigen$vectors[,1:10]
  d <- cov_eigen$values[1:10] ^ 2

  expect_snapshot_value(alpha_log_full_conditional(alpha = -4:5/10,
                                                   O = y_pois,
                                                   U = U,
                                                   d = d,
                                                   n_s = 100,
                                                   xbeta = X %*% beta,
                                                   current_sigma2_tm = 1,
                                                   dens_fun_log = dens_fun_log), style = "serialize")
})


test_that("stglmm phi_tm full conditional correct", {
  set.seed(451)
  dens_fun_log <- function(x, mean) {dpois(x, lambda = poisson()$linkinv(mean), log = TRUE)}

  cov_eigen <- eigen(exp_cov)

  U <- cov_eigen$vectors[,1:10]
  d <- cov_eigen$values[1:10] ^ 2

  expect_snapshot_value(phi_tm_log_full_conditional(phi = 0.2,
                                                    dist_time = as.matrix(dist(time)),
                                                    O = y_pois,
                                                    n_s = n_s,
                                                    n_t = n_t,
                                                    rank = 10,
                                                    nu = nu,
                                                    cores = 1,
                                                    dens_fun_log = dens_fun_log,
                                                    U1 = U,
                                                    xbeta = X %*% beta,
                                                    current_alpha = -4:5/10,
                                                    current_sigma2_tm = 1), style = "serialize")})
