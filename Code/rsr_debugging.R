library(dplyr)
library(ggplot2)
library(here)
library(ngspatial)
library(MASS)
library(broom)
library(spaMM)
library(fields)
library(Rcpp)
library(glue)
library(bayesplot)


source("../RSRcode/spatialPoisson.R")




n <- 500


set.seed(454)

x1 <- runif(n)
x2 <- runif(n)

nu <- 2.5
sigma2 <- 1
phi <- 0.2

tau2 <- 0.1

h <- as.matrix(dist(data.frame(x1 = x1, x2 = x2)))

matern_covariance <- sigma2 * (1 + sqrt(5)*h/phi + 5*h^2 / (3*phi^2)) * exp(-1*sqrt(5)*h/phi)


W <- mvrnorm(n = 1, mu = rep(0, n), Sigma = matern_covariance)



y_mean <- exp(x1 + x2 + W)
y <- rpois(n, y_mean)


sim_data <- tibble(x1 = x1, x2 = x2, W = W, y = y)


sim_data %>%
  ggplot() +
  geom_point(aes(x = x1, y = x2, color = W)) +
  theme_minimal() +
  scale_color_viridis_c(option = "C")

linear_model <- glm(y ~ x1 + x2,
                    data = sim_data,
                    family = "poisson")
linear_model %>%
  tidy(conf.int = TRUE)

sglmm <- fitme(y ~ x1 + x2 + Matern(1|x1 + x2), data = sim_data, family = poisson)
summary(sglmm)


design_matrix <- model.matrix(~ x1 + x2,
                              data = sim_data)

covfn <- covfndef(nu)
adapt <- c("batchlength" = 1e4,
           "n.batch" = 100) # set Monte Carlo sample size = batchlength*n.batch

starting <- list("beta" = coef(linear_model),
                 "s2" = 2,
                 "phi" = 0.5) # set starting values

tuning   <- list("beta" = c(sqrt(diag(vcov(linear_model)))),
                 "s2" = 0.2,
                 "phi" = 0.05,
                 "w" = 0.1) # set tuning parameters

priors   <- list("beta.normal" = c(100),
                 "s2.IG" = c(2,2),
                 "phi.Unif" = c(0.01, 1.5)) # set priors

core = 1
mul = 2
rank = 10


rsr_fit <- poi_gp_arrpfit(O = y,
                          X = design_matrix,
                          coords = cbind(x1, x2),
                          covfn = covfn,
                          adapt = adapt,
                          nu = nu,
                          starting = starting,
                          tuning = tuning,
                          priors = priors,
                          core = core,
                          mul = mul,
                          rank = rank)




model_frame <- lm(~ x1 + x2, data = sim_data, method = "model.frame")

model_matrix <- model.matrix(model_frame)


rrp_output <- rrp_glm(
          fixed = y ~ x1 + x2,
          spatial = ~ x1 + x2,
          data = sim_data,
          family = poisson(),
          covfn = covfndef(nu),
          iter = 2000,
          chains = 2,
          cores = 1,
          param_start = starting,
          priors = priors, # uniform prior on phi
          tuning = tuning,
          nu = nu,
          rank = 10,
          mul = 2)


mcmc_trace(rrp_output$samples,
           pars = glue("delta_3"))

