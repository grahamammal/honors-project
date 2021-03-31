library(purrr)
library(tidyr)
library(glue)
library(rstanarm)
library(rsrHonors)
library(MASS)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)


# simulate data
set.seed(451)

numCores <- detectCores()

registerDoParallel(numCores-1)
run_sims <- function(phi) {
start_time <- Sys.time()
sim_runs <- foreach (i = 1:100, .verbose = T) %dopar% {
  ##############################################
  # Simulate Data
  ##############################################
  n_s <- 100
  n_t <- 50


  nu <- 2.5
  sigma2_sp <- 1
  sigma2_tm <- 1
  phi_sp <- phi
  phi_tm <- phi
  tau2 <- 0.1

  x1 <- rep(runif(n_s), times = n_t)
  x2 <- rep(runif(n_s), times = n_t)
  x3 <- rnorm(n_s*n_t, sd = 0.05)
  x4 <- rnorm(n_s*n_t, sd = 0.05)
  x5 <- rnorm(n_s*n_t, sd = 0.05)

  time <- 1:n_t / n_t

  locations <- data.frame(x1 = x1[1:n_s], x2 = x2[1:n_s])


  X <- cbind(x3, x4, x5)
  beta <- c(1,1,-1)


  d <- as.matrix(dist(data.frame(x1 = x1[1:n_s],
                                 x2 = x2[1:n_s])))
  matern_covariance <- sigma2_sp * (1 + sqrt(5)*d/phi_sp + 5*d^2 / (3*phi_sp^2)) * exp(-1*sqrt(5)*d/phi_sp)

  h <- as.matrix(dist(time))
  exp_cov <- sigma2_tm * exp(-(h/phi_tm)^2)

  W <- mvrnorm(1, rep(0, n_s), matern_covariance)
  V <- mvrnorm(1, rep(0, n_t), exp_cov)


  y_mean <- 1 + as.numeric(X %*% beta + rep(W, times = n_t) + rep(V, each = n_s))

  y_pois <- rpois(n_s * n_t, exp(y_mean))
  y_norm <- rnorm(n_s * n_t, y_mean, sqrt(tau2))
  y_binom <- rbinom(n_s * n_t, 1, exp(y_mean)/(1 + exp(y_mean)))

  sim_data <- data.frame(x1 = x1, x2 = x2, x3 = x3,
                         x4 = x4, x5 = x5,
                         W = rep(W, times = n_t),
                         V = rep(V, each = n_s),
                         y_mean = y_mean,
                         y_pois = y_pois,
                         y_norm = y_norm,
                         y_binom = y_binom)



  ################################################################
  # Fit GLM (negative binomial link function)
  ################################################################

  my_start_time_stan <- Sys.time()
  simulated_glm_negbin <- stan_glm(y_pois ~ x3 + x4 + x5,
                                   data = sim_data,
                                   family = neg_binomial_2(), iter = 5000, chains = 2,
                                   refresh = 0)
  my_end_time_stan <- Sys.time()
  simulated_glm_negbin_draws <- as.matrix(simulated_glm_negbin)

  glm_params_quantiles <- map_dfr(c("(Intercept)", "x3", "x4", "x5"), ~quantile(simulated_glm_negbin_draws[,.x], probs = c(0.05, 0.25, 0.5, 0.75, 0.95))) %>%
    mutate(variable = c("beta_0", "beta_1", "beta_2", "beta_3"))




  ################################################################
  # Fit STGLMM
  ################################################################
  linmod <- glm(y_pois ~ x3 + x4 + x5, data = sim_data, family=poisson(link="log")) # find starting values

  my_time_start <- Sys.time()
  simulated_stglmm <- stglmm(
    fixed = y_pois ~ x3 + x4 + x5,
    data = sim_data,
    locations = locations,
    times = time,
    family = poisson(),
    covfn = 2.5,
    iter = 30000,
    chains = 2,
    cores = 1,
    param_start = list("beta" = linmod$coefficients,
                       "s2_sp" = 3,
                       "phi_sp" = 0.5,
                       "s2_tm" = 1,
                       "phi_tm" = 0.5),
    priors = list("beta.normal" = 100, # variance of beta prior
                  "s2_sp_IG" = c(2,2), # inverse gamma params
                  "s2_tm_IG" = c(2,2),
                  "phi_sp_unif" = c(0.01, 1.5),
                  "phi_tm_unif" = c(0.01, 1.5)), # uniform prior on phi
    tuning = list("beta" = c(sqrt(diag(vcov(linmod)))[2]/4, sqrt(diag(vcov(linmod)))[2:4]),
                  "s2_sp" = 1.5,
                  "phi_sp" = 0.02,
                  "s2_tm" = 1.5,
                  "phi_tm" = 0.01,
                  "delta" = 0.03,
                  "alpha" = 0.02),
    nu = nu,
    rank_sp = 50,
    rank_tm = 10)
  my_time_end <- Sys.time()

  non_burnout_stglmm_samples <- simulated_stglmm$samples[(nrow(simulated_stglmm$samples)/2):nrow(simulated_stglmm$samples),,]

  params_of_interest <- non_burnout_stglmm_samples[,,c("beta_0", "beta_1", "beta_2", "beta_3",
                                 "sigma2_sp", "sigma2_tm",
                                 "phi_sp", "phi_tm")]

  stglmm_params_quantiles <- map_dfr(c("beta_0", "beta_1", "beta_2", "beta_3",
        "sigma2_sp", "sigma2_tm",
        "phi_sp", "phi_tm"), ~quantile(params_of_interest[,,.x], probs = c(0.05, 0.25, 0.5, 0.75, 0.95))) %>%
    mutate(variable = c("beta_0", "beta_1", "beta_2", "beta_3",
                        "sigma2_sp", "sigma2_tm",
                        "phi_sp", "phi_tm"))
  

  list(glm_params_quantiles = glm_params_quantiles,
       stglmm_params_quantiles = stglmm_params_quantiles)
  
}
end_time <- Sys.time()

end_time - start_time

save(sim_runs, file = glue("simulations/simulations_phi_{phi}.RData"))
}

run_sims(0.2)
run_sims(0.5)
run_sims(1)