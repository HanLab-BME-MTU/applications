#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs,mxArray **plhs,int nrhs,const mxArray **prhs){

  if( nrhs != 1 )
    mexErrMsgTxt("poissonRV takes exactly 2D array as input\n");

  size_t nDim=mxGetNumberOfDimensions(prhs[0]);

  if( nDim != 2 )
    mexErrMsgTxt("poissonRV operates on 2D arrays only\n");

  size_t nx,ny;
  size_t i,j;

  double lambda;
  long int xi;


  /* get columns and rows of input array */
  ny=mxGetN(prhs[0]);
  nx=mxGetM(prhs[0]);

  gsl_rng *RNG;
  gsl_matrix *A;
  gsl_matrix_long *B;

  RNG=gsl_rng_alloc(gsl_rng_mt19937);
  A=gsl_matrix_alloc(nx,ny);
  B=gsl_matrix_long__alloc(nx,ny);

  /* copy data to C array */
  memcpy(A->data,mxGetPr(prhs[0]),sizeof(double)*nx*ny);

  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      lambda=gsl_matrix_get(A,i,j);
      xi=gsl_ran_poisson(RNG,lambda);
      gsl_matrix_set(B,i,j,xi);
    }
  }

  /* output */
  if( nlhs > 1 )
    mexErrMsgTxt("poissonRV must have one putput argument\n");

  plhs[0]=mxCreateDoubleMatrix(ny,nx,mxREAL);
  memcpy(mxGetPr(plhs[0]),B->data,sizeof(long int)*nx*ny);

  gsl_rng_free(RNG);
  gsl_matrix_free(A);
}
