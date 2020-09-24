#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat makeCov1(arma::mat coords, double phi, int n, int num);

List rp(double phi_r, // phi_starting
        arma::mat coords, //x,y coords of obs
        int n_r, // batchlength
        int r_r, // rank*mul, not sure what mul is, rank is number of svd vectors.
        double nu_r, // nu
        int num_r // number of cores
          );

arma::mat mmC(arma::mat A,  arma::mat B, int n, int b, int num);
