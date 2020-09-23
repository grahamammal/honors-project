context("Random Projections")


test_that("rp returns list", {
  expect_type(rp(phi_r = 0.2,
                 coords = cbind(rnorm(100), rnorm(100)),
                 n_r = 100,
                 r_r = 50,
                 nu_r = 2.5,
                 num_r = 1), "list")
})
