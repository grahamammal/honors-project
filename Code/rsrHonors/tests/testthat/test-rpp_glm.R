context("rrp_glm")

# simulate data
set.seed(451)
n <- 100
nu <- 2.5
sigma2 <- 1
phi <- 0.2
tau2 <- 0.1


x1 <- runif(n)
x2 <- runif(n)
x3 <- rnorm(n)

h <- as.matrix(dist(data.frame(x1 = x1, x2 = x2)))
matern_covariance <- sigma2 * (1 + sqrt(5)*h/phi + 5*h^2 / (3*phi^2)) * exp(-1*sqrt(5)*h/phi)

W <- rmvn(n = 1, mu = rep(0, n), V = matern_covariance)

y_mean <- 5*x1 + 5*x2 + 5*x3 + W

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
          spatial = y_pois ~ x1 + x2,
          data = sim_data,
          family = poisson(),
          covfn = covfndef(nu),
          iter = 2000,
          chains = 2,
          cores = 1,
          param_start = list("beta" = rnorm(3, 5, sd = 0.5),
                             "s2" = 2,
                             "phi" = 0.5),
          priors = list("beta.normal" = 100, # variance of beta prior
                        "s2.IG" = c(2,2), # inverse gamma params
                        "phi.unif" = c(0.01, 1)), # uniform prior on phi
          tuning = list("beta" = rnorm(3, 1, 0.2),
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
        spatial = y_norm ~ x1 + x2,
        data = sim_data,
        family = poisson(),
        covfn = covfndef(nu),
        iter = 2000,
        chains = 2,
        cores = 1,
        param_start = list("beta" = rnorm(3, 5, sd = 0.5),
                           "s2" = 2,
                           "phi" = 0.5),
        priors = list("beta.normal" = 100, # variance of beta prior
                      "s2.IG" = c(2,2), # inverse gamma params
                      "phi.unif" = c(0.01, 1)), # uniform prior on phi
        tuning = list("beta" = rnorm(3, 1, 0.2),
                      "s2" = 0.1,
                      "phi" = 0.01,
                      "w" = 0.1),
        nu = nu,
        rank = 10,
        mul = 2)

  expect_type(linear_model, "list")
})

test_that("test logistic model runs", {
  logistic_model <- rrp_glm(
        fixed = y_binom ~ x1 + x2 + x3,
        spatial = y_binom ~ x1 + x2,
        data = sim_data,
        family = poisson(),
        covfn = covfndef(nu),
        iter = 2000,
        chains = 2,
        cores = 1,
        param_start = list("beta" = rnorm(3, 5, sd = 0.5),
                           "s2" = 2,
                           "phi" = 0.5),
        priors = list("beta.normal" = 100, # variance of beta prior
                      "s2.IG" = c(2,2), # inverse gamma params
                      "phi.unif" = c(0.01, 1)), # uniform prior on phi
        tuning = list("beta" = rnorm(3, 1, 0.2),
                      "s2" = 0.1,
                      "phi" = 0.01,
                      "w" = 0.1),
        nu = nu,
        rank = 10,
        mul = 2)

  expect_type(logistic_model, "list")
})

test_that("poisson fit matches spatialPoisson fit", {
  starting <- list("beta" = rnorm(3, 5, sd = 0.5),
                   "s2" = 2,
                   "phi" = 0.5)

  priors <- list("beta.normal" = 100, # variance of beta prior
                 "s2.IG" = c(2,2), # inverse gamma params
                 "phi.unif" = c(0.01, 1))

  tuning <- list("beta" = rnorm(3, 1, 0.2),
                 "s2" = 0.1,
                 "phi" = 0.01,
                 "w" = 0.1)

  rank = 10

  set.seed(451)
  my_fit <- rrp_glm(
      fixed = y_pois ~ x1 + x2 + x3,
      spatial = y_pois ~ x1 + x2,
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
        core = core,
        starting = starting,
        tuning = tuning,
        priors = priors,
        nu = nu,
        rank = 10)

  expect_equal(my_fit$p.params, orig_fit$p.params)
})

test_that("beta full conditional correct", {

})


test_that("delta full conditional correct", {
})

test_that("phi full conditional correct", {

})
