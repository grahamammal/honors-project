#include <algorithm>
#include <string>
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
//#include <malloc.h>
//#include <cuda_runtime.h>
//#include "cublas_v2.h"
#include <Rmath.h>
#include <Rinternals.h>
// #include <R_ext/Linpack.h>
// #include <R_ext/Lapack.h>
// #include <R_ext/BLAS.h>
//#include <omp.h>
//#include <mkl.h>
//#include <mkl_lapacke.h>
//#include <mkl_blas.h>
//#include <Accelerate/Accelerate.h>
#include <vecLib/clapack.h>
#include <vecLib/cblas.h>

#define pi 3.14159265358979323846

extern "C" {
  
  /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
  /*::  Function prototypes                                           :*/
  /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
  double deg2rad(double);
  double rad2deg(double);
  
  double distance(double lat1, double lon1, double lat2, double lon2){
    double dist;
    double R =  6378.388;
    
    dist = cos(deg2rad(lat1)) * cos(deg2rad(lon1)) * cos(deg2rad(lat2)) * cos(deg2rad(lon2))+sin(deg2rad(lon1))* sin(deg2rad(lon2))*cos(deg2rad(lat1))*cos(deg2rad(lat2)) + sin(deg2rad(lat2))*sin(deg2rad(lat1));
    dist = acos(fabs(dist) > 1.0 ? 1.0 * (dist > 0.0 ? 1.0 : -1.0) : dist);
    dist = dist * R; 
    return(dist);
  }
  
  /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
  /*::  This function converts decimal degrees to radians             :*/
  /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
  double deg2rad(double deg){
    return(deg * pi / 180.0);
  }
  
  /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
  /*::  This function converts radians to decimal degrees             :*/
  /*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
  double rad2deg(double rad){
    return(rad * 180.0 / pi);
  }
  
  
  // calculate the earth distance induced covariance matrix for the bird data
  //calculate the exponential cov fn
  double* makeCov000(double *coords, double phi,int n,int num){
    // double oneoverphi = -1/phi;
    double *C = (double *) R_alloc(n*n, sizeof(double));
    //omp_set_num_threads(num);	
    int i,j;
    double dist = 0.0;
    
    #pragma omp parallel shared(coords) private(i,j) 
    {
      #pragma omp for schedule(static)        	
      for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
          dist = distance(coords[i], coords[n+i],coords[j], coords[n+j]);
          // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
          C[n*j+i] = C[n*i+j] = exp(-dist/phi);
        }	
      }
    }
    return(C);
  }
  
  // calculate the covariance matrix for the digit data
  double* makeCov00(double *coords, double phi,int n,int num){
    double oneoverphi = -1/(2*pow(phi,2));
    double *C = (double *) R_alloc(n*n, sizeof(double));
    //omp_set_num_threads(num);	
    int i,j,k;
    double dist = 0.0;
    
    #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
    {
      #pragma omp for schedule(static)        	
      for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
          dist = 0.0;
          for(k = 0; k < 256; k++){
            dist += pow(coords[k*n+i]-coords[k*n+j],2);
          }
          // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*dist);
          C[n*j+i] = C[n*i+j] = exp(oneoverphi*dist);
        }	
      }
    }
    return(C);
  }
  
  // calculate the covariance matrix for the gene data matern 1.5
  double* makeCovGene(double *coords, double phi,int n,int num){
    double oneoverphi = sqrt(3)/phi;//-1/(2*pow(phi,2));
    double *C = (double *) R_alloc(n*n, sizeof(double));
    //omp_set_num_threads(num);	
    int i,j,k;
    double dist = 0.0;
    
    #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
    {
      #pragma omp for schedule(static)        	
      for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
          dist = 0.0;
          for(k = 0; k < 19; k++){
            dist += pow(coords[k*n+i]-coords[k*n+j],2);
          }
          dist = sqrt(dist);
          // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*dist);
          C[n*j+i] = C[n*i+j] = (1.0+ oneoverphi*dist)*exp(-1*oneoverphi*dist);//exp(oneoverphi*dist);
        }	
      }
    }
    return(C);
  }
  // calculate the sq exp covariance matrix for 2-d 
  double* makeCov0(double *coords, double phi,int n,int num){
    double oneoverphi = -1/(2*pow(phi,2));
    double *C = (double *) R_alloc(n*n, sizeof(double));
    //omp_set_num_threads(num);	
    int i,j;
    
    #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
    {
      #pragma omp for schedule(static)        	
      for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
          // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
          C[n*j+i] = C[n*i+j] = exp(oneoverphi*(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
        }	
      }
    }
    return(C);
  }	
  //calculate the exponential cov fn
  double* makeCov(double *coords, double phi,int n,int num){
    double oneoverphi = -1/phi;
    double *C = (double *) R_alloc(n*n, sizeof(double));
    //omp_set_num_threads(num);	
    int i,j;
    
    #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
    {
      #pragma omp for schedule(static)        	
      for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
          // REAL(C)[n*j+i] = REAL(C)[n*i+j] = exp(oneoverphi*sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
          C[n*j+i] = C[n*i+j] = exp(oneoverphi*sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2)));
        }	
      }
    }
    return(C);
  }
  
  // calcualte matern 2.5
  double* makeCov1(double *coords, double phi,int n,int num){
    double oneoverphi = sqrt(5)/phi;
    int i,j;
    double dist;
    double *C = (double *) R_alloc(n*n, sizeof(double));
    //omp_set_num_threads(num);	
    
    #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
    {
      #pragma omp for schedule(static)        	
      for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
          dist = sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2));
          C[n*j+i] = C[n*i+j] = (1.0+ oneoverphi*dist + 1.0/3.0*pow(oneoverphi,2)*pow(dist,2))*exp(-1*oneoverphi*dist);
        }	
      }
    }
    return(C);
  }
  // calcualte matern 1.5
  double* makeCov15(double *coords, double phi,int n,int num){
    double oneoverphi = sqrt(3)/phi;
    int i,j;
    double dist;
    double *C = (double *) R_alloc(n*n, sizeof(double));
    //omp_set_num_threads(num);	
    
    #pragma omp parallel shared(oneoverphi,coords) private(i,j) 
    {
      #pragma omp for schedule(static)        	
      for(i = 0; i < n; i++){
        for(j = i; j < n; j++){
          dist = sqrt(pow(coords[i]-coords[j],2) + pow(coords[n+i]-coords[n+j],2));
          C[n*j+i] = C[n*i+j] = (1.0+ oneoverphi*dist)*exp(-1*oneoverphi*dist);
        }	
      }
    }
    return(C);
  }
  
  // dyn.load(paste0(sourceDir,"phiC.so"))
  // .Call("rp",0.2,coords,as.integer(n),as.integer(r),2.5,as.integer(num))
  SEXP rp(SEXP phi_r, SEXP coords_r, SEXP n_r, SEXP r_r, SEXP nu_r, SEXP num_r){
    
    int i,j,k;
    int n = INTEGER(n_r)[0];
    int r = INTEGER(r_r)[0];
    double sum;
    int num = INTEGER(num_r)[0];
    double nu = REAL(nu_r)[0];
    double sd0 = 1/sqrt(r);
    double phi = REAL(phi_r)[0];
    double *coords = REAL(coords_r);
    double *(*covfn) (double *coords,double phi ,int n,int num);
    double *K = (double *) R_alloc(n*n, sizeof(double));
    double *omega = (double *) R_alloc(n*r, sizeof(double));
    double *Komega=(double*) R_alloc(n*r, sizeof(double));
    double *Kphit=(double*) R_alloc(n*r, sizeof(double));
    double ALPHA;
    ALPHA= 1.0;
    
    if(nu==200) covfn = makeCov000; // only for the bird count
    if(nu==100) covfn = makeCov00; // only for the binary digit classification
    if(nu==19) covfn = makeCovGene; // only for the binary gene classification
    if(nu==10) covfn = makeCov0; // for square exponential
    if(nu==0.5) covfn = makeCov; // matern 
    if(nu==2.5) covfn = makeCov1;// matern 
    if(nu==1.5) covfn = makeCov15;// matern 
    K = covfn(coords ,  phi , n , num );
    
    for(i = 0; i < n; i++){
      for(j = 0; j < r; j++){
        omega[n*j+i] = rnorm(0,sd0);
      }
    }
    
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
    n, r, n, 1.0 , K, n, omega, n, 0.0 ,Komega, n);	 // n*r mm multiplication
    
    
    double *d=(double*) R_alloc(r,sizeof(double)); //store eigenvalues      
    double *u=(double*) R_alloc(n*n,sizeof(double)); // not referenced
    double *VT=(double*) R_alloc(r*r,sizeof(double)); // not used
    int lwork;
    double wkopt;
    double *work;
    int *iwork = (int*) R_alloc(8*r, sizeof(double));
    int info;
    
    lwork = -1;
    F77_NAME(dgesdd)("O", &n, &r, Komega, &n, d, u, &n, VT, &r, &wkopt, &lwork, iwork, &info);
    lwork = (int) wkopt;
    work = (double*) R_alloc(lwork, sizeof(double));
    F77_NAME(dgesdd)("O", &n, &r, Komega, &n, d, u, &n, VT, &r, work, &lwork, iwork, &info);
    
    if (info > 0) {
      printf("SVN failed to converge");
    }
	
    //Komega is overwritten by the basis phit
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
    n, r, n, 1.0 , K, n, Komega, n, 0.0 ,Kphit, n);	 // n*r mm multiplication
    
	// Komega transpose * Kphit = rxn * nxr -> phi*K*phit = VT
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
    r, r, n, 1.0 , Komega, n, Kphit, n, 0.0 ,VT, r);	 // rxn * nxr mm multiplication
    
	// chol(VT) upper VT is replaced by cholesky upper tri-matrix    
    F77_NAME(dpotrf)("U", &r, VT, &r, &info); if(info != 0){error("c++ error: dpotrf failed\n");} 
    // matrix solve
    //dtrsm("R", "U", "N", "N", &n, &r, &ALPHA, VT, &r, Kphit, &n);
    cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans,
                CblasNonUnit, n, r, ALPHA, VT, r, Kphit, n);
    
	SEXP d_r,u_r,result;
     d_r = PROTECT(allocMatrix(REALSXP, r, 1));
     u_r = PROTECT(allocMatrix(REALSXP, n,r));
	
	// svd Kphit
	 F77_NAME(dgesdd)("S", &n, &r, Kphit, &n, REAL(d_r),  REAL(u_r), &n, VT, &r, work, &lwork, iwork, &info);
     int nResultListObjs = 2;
    
	// comment this part for testing
    PROTECT(result = allocVector(VECSXP, nResultListObjs));
    SET_VECTOR_ELT(result, 0, d_r);
    SET_VECTOR_ELT(result, 1, u_r);
	UNPROTECT(3);
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
  
  // matrix multiplication nxn *nxr
  SEXP mmC(SEXP A_r,  SEXP B_r, SEXP n_r, SEXP b_r,SEXP num_r){
    int i,j;
    double alpha, beta;
    int n = INTEGER(n_r)[0];
    int b = INTEGER(b_r)[0];	
    int num = INTEGER(num_r)[0];
    double *A = REAL(A_r);
    double *B = REAL(B_r);
    double sum;
    alpha = 1.0; beta = 0.0;
    
    SEXP Cfoo;
    Cfoo = PROTECT(allocMatrix(REALSXP, n, b));
    
    //double *C = (double *)malloc( n*b*sizeof( double ), 64 );
    double *C = (double *)malloc( n*b*sizeof( double ) );
    //omp_set_num_threads(num);
    
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
    n, b, n, alpha, A, n, B, n, beta,C, n);	
    
    #pragma omp parallel shared(C,Cfoo,n,b) private(i, j)
    {
      
      #pragma omp for schedule(static)        	
      for(i = 0; i < n; i++){
        for(j = 0; j < b; j++){
          REAL(Cfoo)[n*j+i] = C[n*j+i] ;
        }
      }
    }
    
    free(C);
    UNPROTECT(1);
    return(Cfoo);
  }
}
