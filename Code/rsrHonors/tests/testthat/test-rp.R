

test_that("rp returns list", {
  expect_type(rp(phi_r = 0.2,
                 coords = cbind(rnorm(100), rnorm(100)),
                 n_r = 100,
                 r_r = 50,
                 nu_r = 2.5,
                 num_r = 1), "list")
})

test_that("rp close to actual values", {
  n = 1000
  coords = cbind(runif(n), runif(n))
  phi = 0.187
  r = n/10
  nu = 2.5

  rp_out = rp(phi_r = phi,
              coords = coords,
              n_r = n,
              r_r = r,
              nu_r = nu,
              num_r = 1)

  h <- as.matrix(dist(coords))

  K <- 1 * (1 + sqrt(5)*h/phi + 5*h^2 / (3*phi^2)) * exp(-1*sqrt(5)*h/phi)

  K_eigen <- eigen(K)

  expect_true(mean(rp_out[[1]]^2 - K_eigen$values[1:100]) < 1)

})
