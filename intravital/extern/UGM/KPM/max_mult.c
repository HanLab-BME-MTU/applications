/* C mex version of max_mult.m in BPMRF2 directory  */
/* gcc -Wall -I/mit/matlab_v6.5/distrib/bin/glnx86 -c max_mult.c */

// MAX_MULT Like matrix multiplication, but sum gets replaced by max
// function y=max_mult(A,x) y(j) = max_i A(j,i) x(i)
// MATLAB version:
//     X=x*ones(1,size(A,1)); % X(i,j) = x(i)
//     y=max(A'.*X)';


#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int rows,cols,common,m,n,p;
  double y1, y2;
  double *arr1, *arr2, *arr3;


  if (nrhs!=2 || nlhs>1)
    mexErrMsgTxt("max_mult requires two inputs and one output");
  if (mxIsChar(prhs[0]) || mxIsClass(prhs[0], "sparse") || mxIsComplex(prhs[0])
      || mxIsChar(prhs[1]) || mxIsClass(prhs[1], "sparse") || mxIsComplex(prhs[1]))
    mexErrMsgTxt("Inputs must be real, full, and nonstring");
  if (mxGetN(prhs[0])!=mxGetM(prhs[1]))
    mexErrMsgTxt("The number of columns of A must be the same as the number of rows of x");
  

  arr1=mxGetPr(prhs[0]);
  arr2=mxGetPr(prhs[1]);
  p=mxGetN(prhs[0]);
  m=mxGetM(prhs[0]);
  n=mxGetN(prhs[1]);
    // correction 
	//plhs[0]=mxCreateDoubleMatrix(m, n, mxREAL);
  //arr3=mxMalloc(m*n*sizeof(double));
  // uday's code
  plhs[0] = mxCreateNumericMatrix( m*n, 1, mxDOUBLE_CLASS, mxREAL);
  arr3 = (double*)mxGetData(plhs[0]);
	//
  for (rows=0; rows<m ; rows++)
    for (cols=0; cols<n ; cols++)
    {
      y1=arr1[rows]*arr2[cols*p];  
      for (common=1; common<p; common++)
      {
	y2=arr1[rows+common*m]*arr2[common+cols*p];
        if (y2>y1)
          y1=y2;
      }
      arr3[rows+cols*m]=y1;
    }

  //mxSetPr(plhs[0], arr3);
	 
}
