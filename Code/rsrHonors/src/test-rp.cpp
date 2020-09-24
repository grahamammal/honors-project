/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>
#include "random_projection.h"

// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.
// int twoPlusTwo() {
//   return 2 + 2;
// }

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("Testing RP C++") {

  // The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
  // test_that("MakeCov2.5 is symmetric") {
  //   expect_true(makeCov1(arma::mat(100, 2, arma::fill::randn),
  //                        0.2,
  //                        100,
  //                        1).is_symmetric());
  // }
//
//   test_that("Test pls") {
//     expect_true(2 + 2 == 4);
//   }
//
//   test_that("Test pls2") {
//     expect_true(2 + 2 != 4);
//   }
//
//   test_that("Test vector norm closeness") {
//     int n = 100;
//
//     arma::mat coords = arma::mat(n, 2);
//
//     arma::mat::iterator coords_begin = coords.begin();
//     arma::mat::iterator coords_end = coords.end();
//
//     // for(; coords_begin != coords_end; ++coords_begin) {
//     //   (*coords_begin) = R::runif(0,1);
//     // }
//     // arma::mat K = arma::mat(n,n);
//     // K = makeCov1(coords, 0.189, n, 1);
//     //
//     // arma::vec real_eigen_val = arma::vec(n);
//     // arma::mat real_eigen_vec = arma::mat(n,n);
//     // arma::mat real_eigen_vec2 = arma::mat(n,n);
//     //
//     // //arma::svd_econ(real_eigen_vec, real_eigen_val, real_eigen_vec2, K, "left");
//     // //List rp_out = rp(0.189, coords, n, n / 10, 2.5, 1);
//     //
//     // //NumericVector est_eigen_val = rp_out[0];
//     //
//     // real_eigen_val.resize(n / 10);
//     // real_eigen_vec.resize(n, n / 10);
//     //
//     // // double sum = 0;
//     // // for(int i = 0; i < 100; i++) {
//     // //   sum += estimated_eigen_val[i];
//     // // }
//     //
//     //
//     // double tmp = coords(0,1);
//     expect_true(coords(0,1) == 1);
//   }

}
