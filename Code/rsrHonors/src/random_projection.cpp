#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#define pi 3.14159265358979323846


// calcualte matern 2.5
arma::mat makeCov1(arma::mat coords, double phi, int n, int num){
  double oneoverphi = sqrt(5)/phi;
  int i,j;
  double dist;
  arma::mat C = arma::mat(n, n, arma::fill::none);


  for(i = 0; i < n; i++){
    for(j = i; j < n; j++){
      dist = sqrt(pow(coords(i,0)-coords(j,0),2) + pow(coords(i,1)-coords(j,1),2));
      C(i,j) = C(j,i) = (1.0+ oneoverphi*dist + 1.0/3.0*pow(oneoverphi,2)*pow(dist,2))*exp(-1*oneoverphi*dist);
    }
  }

  return(C);
}

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
List rp(double phi_r, // phi_starting
        arma::mat coords, //x,y coords of obs
        int n_r, // batchlength
        int r_r, // rank*mul, not sure what mul is, rank is number of svd vectors.
        double nu_r, // nu
        int num_r // number of cores
          ){

  int i,j,k;
  int n = n_r; // Cast into int
  int r = r_r; // Cast into int
  double sum;
  int num = num_r; // Cast into int
  double nu = nu_r; // cast into double
  double sd0 = 1/sqrt(r); //
  double phi = phi_r; // cast to double

  arma::mat (*covfn) (arma::mat coords, double phi, int n, int num); // allocate create pointer to covariance matrix creating function
  arma::mat K = arma::mat(n, n); // pointer to covariance matrix
  arma::mat Omega = arma::mat(n, r); // Make pointer to double matrix of size n*r
  arma::mat KOmega= arma::mat(n, r); // Make pointer to double matrix of size n*r
  arma::mat K11 = arma::mat(r, r);
  arma::mat KPhit = arma::mat(n, r);
  double ALPHA; // constant, not sure what it does
  ALPHA= 1.0;


  // Make function to make K
  if(nu==2.5) covfn = makeCov1;// matern


  // Make K (Covariance matrix of Coords)
  K = covfn(coords ,  phi , n , num); // K is the covariance matrix

  // Build Omega (random projection matrix
  for(i = 0; i < n; i++){
    for(j = 0; j < r; j++){
      Omega(i,j) = as<double>(rnorm(1, 0, sd0)); // simulate omega used to project into lower dim
    }
  }

  // Called Phi in paper, last part of step 1 in algo 1
  // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
  //   n, r, n, 1.0 , K, n, omega, n, 0.0 ,Komega, n);	 // n*r mm multiplication
  //
  KOmega = K * Omega;

  arma::vec d = arma::vec(r); //store eigenvalues
  arma::mat u= arma::mat(r, r); // not referenced
  arma::mat VT= arma::mat(r, r); // not used


  // F77_NAME(dgesdd)("O", &n, &r, Komega, &n, d, u, &n, VT, &r, &wkopt, &lwork, iwork, &info);
  // Komega contains u now after above computation
  svd_econ(u, d, VT, KOmega, "left");

  // F77_NAME(dgesdd)("O", &n, &r, Komega, &n, d, u, &n, VT, &r, work, &lwork, iwork, &info);
  // Komega is actually u from before. Now, performing svd on this, and again overwriting Komega with the new u
  svd_econ(u, d, VT, u, "left");


  // Komega is overwritten by the basis phit
  // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
  //   n, r, n, 1.0 , K, n, Komega, n, 0.0 ,Kphit, n);	 // n*r mm multiplication
  KPhit = K * u;


  // Komega transpose * Kphit = rxn * nxr -> phi*K*phit = VT
  // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
  // r, r, n, 1.0 , Komega, n, Kphit, n, 0.0 ,VT, r);	 // rxn * nxr mm multiplication
  // recall, Komega above is actually left sing values of previous svd
  VT = u.t() * KPhit;

  // chol(VT) upper VT is replaced by cholesky upper tri-matrix
  // F77_NAME(dpotrf)("U", &r, VT, &r, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
  VT = chol(VT, "upper");

  Rcout << VT

  // matrix solve
  // cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans,
  //   CblasNonUnit, n, r, ALPHA, VT, r, Kphit, n);
  KPhit = solve(VT, KPhit);

  List result;

  arma::vec d_r = arma::vec(r, 1, arma::fill::none);
  arma::mat u_r = arma::mat(n, r, arma::fill::none);

  // svd Kphit
  // F77_NAME(dgesdd)("S", &n, &r, Kphit, &n, REAL(d_r),  REAL(u_r), &n, VT, &r, work, &lwork, iwork, &info);
  // int nResultListObjs = 2;
  svd_econ(u_r, d_r, VT, KPhit, "left");

  int nResultListObjs = 2;



// comment this part for testing
  result = List::create(d_r, u_r);
  return(result);

// uncomment this part for testing
  /* // // Kphit becomes the final matrix
   SEXP out;
   out = PROTECT(allocMatrix(REALSXP, n, r));



  #pragma omp parallel shared(out,Kphit,r,n) private(i, j)
  {

   #pragma omp for schedule(static)
    for(i = 0; i < n; i++){
      for(j = 0; j < r; j++){
        REAL(out)[n*j+i] = Kphit[n*j+i] ;
      }
    }
  }

  UNPROTECT(1);
  return(out); */

}

//' Matrix Multiplication
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
arma::mat mmC(arma::mat A,  arma::mat B, int n, int b, int num){
  int i,j;
  double alpha, beta;
  double sum;
  alpha = 1.0; beta = 0.0;

  //double *C = (double *)malloc( n*b*sizeof( double ), 64 );
  arma::mat C = arma::mat(n, b, arma::fill::none);
  //omp_set_num_threads(num);

  C = A * B;

  return(C);
}


