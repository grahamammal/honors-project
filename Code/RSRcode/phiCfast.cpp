#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

#define pi 3.14159265358979323846



/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  Function prototypes                                           :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
// double deg2rad(double);
// double rad2deg(double);
// 
// double distance(double lat1, double lon1, double lat2, double lon2){
//   double dist;
//   double R =  6378.388;
//   
//   dist = cos(deg2rad(lat1)) * cos(deg2rad(lon1)) * cos(deg2rad(lat2)) * cos(deg2rad(lon2))+sin(deg2rad(lon1))* sin(deg2rad(lon2))*cos(deg2rad(lat1))*cos(deg2rad(lat2)) + sin(deg2rad(lat2))*sin(deg2rad(lat1));
//   dist = acos(fabs(dist) > 1.0 ? 1.0 * (dist > 0.0 ? 1.0 : -1.0) : dist);
//   dist = dist * R; 
// }
// 
// /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
// /*::  This function converts decimal degrees to radians             :*/
// /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
// double deg2rad(double deg){
//   return(deg * pi / 180.0);
// }
// 
// /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
// /*::  This function converts radians to decimal degrees             :*/
// /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
// double rad2deg(double rad){
//   return(rad * 180.0 / pi);
// }
// 
// 
// // calculate the earth distance induced covariance matrix for the bird data
// //calculate the exponential cov fn
// double* makeCov000(double *coords, double phi,int n,int num){
//   // double oneoverphi = -1/phi;
//   double *C = (double *) R_alloc(n*n, sizeof(double));
//   //omp_set_num_threads(num);	
//   int i,j;
//   double dist = 0.0;
//   
//   #pragma omp parallel shared(coords) private(i,j) 
//   {
//     #pragma omp for schedule(static)        	
//     for(i = 0; i < n; i++){
//       for(j = i; j < n; j++){
//         dist = distance(coords[i], coords[n+i],coords[j], coords[n+j]);
//         // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
//         C[n*j+i] = C[n*i+j] = exp(-dist/phi);
//       }	
//     }
//   }
//   return(C);
// }
// 
// // calculate the covariance matrix for the digit data
// double* makeCov00(double *coords, double phi,int n,int num){
//   double oneoverphi = -1/(2*pow(phi,2));
//   double *C = (double *) R_alloc(n*n, sizeof(double));
//   //omp_set_num_threads(num);	
//   int i,j,k;
//   double dist = 0.0;
//   
//   #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
//   {
//     #pragma omp for schedule(static)        	
//     for(i = 0; i < n; i++){
//       for(j = i; j < n; j++){
//         dist = 0.0;
//         for(k = 0; k < 256; k++){
//           dist += pow(coords[k*n+i]-coords[k*n+j],2);
//         }
//         // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*dist);
//         C[n*j+i] = C[n*i+j] = exp(oneoverphi*dist);
//       }	
//     }
//   }
//   return(C);
// }
// 
// // calculate the covariance matrix for the gene data matern 1.5
// double* makeCovGene(double *coords, double phi,int n,int num){
//   double oneoverphi = sqrt(3)/phi;//-1/(2*pow(phi,2));
//   double *C = (double *) R_alloc(n*n, sizeof(double));
//   //omp_set_num_threads(num);	
//   int i,j,k;
//   double dist = 0.0;
//   
//   #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
//   {
//     #pragma omp for schedule(static)        	
//     for(i = 0; i < n; i++){
//       for(j = i; j < n; j++){
//         dist = 0.0;
//         for(k = 0; k < 19; k++){
//           dist += pow(coords[k*n+i]-coords[k*n+j],2);
//         }
//         dist = sqrt(dist);
//         // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*dist);
//         C[n*j+i] = C[n*i+j] = (1.0+ oneoverphi*dist)*exp(-1*oneoverphi*dist);//exp(oneoverphi*dist);
//       }	
//     }
//   }
//   return(C);
// }
// // calculate the sq exp covariance matrix for 2-d 
// double* makeCov0(double *coords, double phi,int n,int num){
//   double oneoverphi = -1/(2*pow(phi,2));
//   double *C = (double *) R_alloc(n*n, sizeof(double));
//   //omp_set_num_threads(num);	
//   int i,j;
//   
//   #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
//   {
//     #pragma omp for schedule(static)        	
//     for(i = 0; i < n; i++){
//       for(j = i; j < n; j++){
//         // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
//         C[n*j+i] = C[n*i+j] = exp(oneoverphi*(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
//       }	
//     }
//   }
//   return(C);
// }	
// //calculate the exponential cov fn
// double* makeCov(double *coords, double phi,int n,int num){
//   double oneoverphi = -1/phi;
//   double *C = (double *) R_alloc(n*n, sizeof(double));
//   //omp_set_num_threads(num);	
//   int i,j;
//   
//   #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
//   {
//     #pragma omp for schedule(static)        	
//     for(i = 0; i < n; i++){
//       for(j = i; j < n; j++){
//         // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
//         C[n*j+i] = C[n*i+j] = exp(oneoverphi*sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
//       }	
//     }
//   }
//   return(C);
// }

// calcualte matern 2.5
mat makeCov1(mat coords, double phi, int n, int num){
  double oneoverphi = sqrt(5)/phi;
  int i,j;
  double dist;
  mat C = mat(n, n, fill::none);
  
      	
  for(i = 0; i < n; i++){
    for(j = i; j < n; j++){
      dist = sqrt(pow(coords(i,1)-coords(j,1),2) + pow(coords(i,2)-coords(j,2),2));
      C(i,j) = C(j,i) = (1.0+ oneoverphi*dist + 1.0/3.0*pow(oneoverphi,2)*pow(dist,2))*exp(-1*oneoverphi*dist);
    }	
  }

  return(C);
}
// // calcualte matern 1.5
// double* makeCov15(double *coords, double phi,int n,int num){
//   double oneoverphi = sqrt(3)/phi;
//   int i,j;
//   double dist;
//   double *C = (double *) R_alloc(n*n, sizeof(double));
//   //omp_set_num_threads(num);	
//   
//   #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
//   {
//     #pragma omp for schedule(static)        	
//     for(i = 0; i < n; i++){
//       for(j = i; j < n; j++){
//         dist = sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2));
//         C[n*j+i] = C[n*i+j] = (1.0+ oneoverphi*dist)*exp(-1*oneoverphi*dist);
//       }	
//     }
//   }
//   return(C);
// }

// dyn.load(paste0(sourceDir,"phiC.so"))
// .Call("rp",0.2,coords,as.integer(n),as.integer(r),2.5,as.integer(num))


// [[Rcpp::export]]
List rp(double phi_r, // phi_starting
        mat coords, //x,y coords of obs
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

  mat (*covfn) (mat coords, double phi, int n, int num); // allocate create pointer to covariance matrix creating function
  mat K = mat(n, n, fill::none); // pointer to covariance matrix
  mat omega = mat(n, r, fill::none); // Make pointer to double matrix of size n*r
  mat Komega= mat(n, r, fill::none); // Make pointer to double matrix of size n*r
  mat Kphit = mat(n, r, fill::none); // Make pointer to double matrix of size n*r
  double ALPHA; // constant, not sure what it does
  ALPHA= 1.0;
  
  // if(nu==200) covfn = makeCov000; // only for the bird count
  // if(nu==100) covfn = makeCov00; // only for the binary digit classification
  // if(nu==19) covfn = makeCovGene; // only for the binary gene classification
  // if(nu==10) covfn = makeCov0; // for square exponential
  // if(nu==0.5) covfn = makeCov; // matern 
  if(nu==2.5) covfn = makeCov1;// matern 
  // if(nu==1.5) covfn = makeCov15;// matern 
  K = covfn(coords ,  phi , n , num ); // K is the covariance matrix
  
  for(i = 0; i < n; i++){
    for(j = 0; j < r; j++){
      omega(i,j) = as<double>(rnorm(1, 0, sd0)); // simulate omega used to project into lower dim
    }
  }
  
  Komega = K * omega;
  
  vec d = vec(r, fill::none); //store eigenvalues      
  mat u= mat(n, n, fill::none); // not referenced
  mat VT= mat(r, r, fill::none); // not used
  int lwork;
  double wkopt;
  double *work;
  int *iwork = (int*) R_alloc(8*r, sizeof(double));
  int info;
  
  lwork = -1;
  // 
  //
  //
  //
  //
  //
  //
  //
  svd_econ(u, d, VT, Komega);
  

  //Komega is overwritten by the basis phit
  
  Kphit = K * Komega;
  
  // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
  // n, r, n, 1.0 , K, n, Komega, n, 0.0 ,Kphit, n);	 // n*r mm multiplication
  
  VT = Komega.t() * Kphit;
  
// Komega transpose * Kphit = rxn * nxr -> phi*K*phit = VT
  // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
  // r, r, n, 1.0 , Komega, n, Kphit, n, 0.0 ,VT, r);	 // rxn * nxr mm multiplication
   
// chol(VT) upper VT is replaced by cholesky upper tri-matrix    
  
  VT = chol(VT, "upper");
  
  // LAPACKE_dpotrf(CblasColMajor,
  //                'U', 
  //                r, 
  //                VT, 
  //                r);
  
  
  Kphit = solve(VT, Kphit);
  // matrix solve
  //dtrsm("R", "U", "N", "N", &n, &r, &ALPHA, VT, &r, Kphit, &n);
  // cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans,
  //             CblasNonUnit, n, r, ALPHA, VT, r, Kphit, n);
  
  List result;
  
  vec d_r = vec(r, 1, fill::none);
  mat u_r = mat(n, r, fill::none);

  
  svd_econ(u_r, d_r, VT, Kphit);
// svd Kphit
 // LAPACKE_dgesdd(CblasColMajor,
 //               'S', 
 //               n, 
 //               r, 
 //               Kphit, 
 //               n, 
 //               REAL(d_r),  
 //               REAL(u_r), 
 //               n, 
 //               VT, 
 //               r);
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

// // matrix multiplication nxn *nxr
// NumericMatrix mmC(SEXP A_r,  SEXP B_r, SEXP n_r, SEXP b_r,SEXP num_r){
//   int i,j;
//   double alpha, beta;
//   int n = INTEGER(n_r)[0];
//   int b = INTEGER(b_r)[0];	
//   int num = INTEGER(num_r)[0];
//   double *A = REAL(A_r);
//   double *B = REAL(B_r);
//   double sum;
//   alpha = 1.0; beta = 0.0;
//   
//   SEXP Cfoo;
//   Cfoo = NumericMatrix(n, b);
//   
//   //double *C = (double *)malloc( n*b*sizeof( double ), 64 );
//   double *C = (double *)malloc( n*b*sizeof( double ) );
//   //omp_set_num_threads(num);
//   
//   cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
//   n, b, n, alpha, A, n, B, n, beta,C, n);	
//   
//   #pragma omp parallel shared(C,Cfoo,n,b) private(i, j)
//   {
//     
//     #pragma omp for schedule(static)        	
//     for(i = 0; i < n; i++){
//       for(j = 0; j < b; j++){
//         REAL(Cfoo)[n*j+i] = C[n*j+i] ;
//       }
//     }
//   }
//   
//   free(C);
//   return(Cfoo);
// }

