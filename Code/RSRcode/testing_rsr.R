library(dplyr)
library(ggplot2)
library(here)
library(ngspatial)


data(infant)

infant_small <- infant[1:500,]

infant_small$low_weight = infant_small$low_weight / infant_small$births
Z = infant_small$deaths
X = cbind(1, infant_small$low_weight, infant_small$black, infant_small$hispanic, infant_small$gini, infant_small$affluence, infant_small$stability)
data(A)

A_small <- A[1:500, 1:500]

# takes 23 minutes to run be careful

set.seed(123456)
fit = sparse.sglmm(Z ~ X - 1 + offset(log(infant_small$births)), family = poisson, A = A_small, method = "RSR",
                   tune = list(sigma.s = 0.02), verbose = TRUE)
summary(fit)


fit_lm <- glm(Z ~ X - 1 + offset(log(infant_small$births)), family = poisson)
summary(fit_lm)
confint(fit_lm)

fit_sparse.sglmm <- fit

# these are completley different results from what the Khan paper says we should see!?!?!