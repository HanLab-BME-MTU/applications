/* ------------------------------------------------------- */
/*                                                         */
/* comp2Vec             [MATLAB C-MEX]                     */
/*                                                         */
/* ------------------------------------------------------- */
/*                                                         */
/* Files:                                                  */
/*                                                         */
/*     comp2Vec.c - MEX interface                          */
/*                                                         */
/*                                                         */
/* First version: Achim Besser - 12/02/2010                */
/* ------------------------------------------------------- */

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    /* Initialize pointers for 2 inputs and 1 output */
    double *M, *N, *T;
	//double *D;

	/* Initialize int variables to store matrix dimensions */
	int Mrows,Mcols;
	int Nrows,Ncols;
    
    /* Initialize int variables for the 2 for loops: */
    int idx,jdx;
    
    /* Check that the number of input and output parameters is valid */
	if(nrhs != 2)
		mexErrMsgTxt("Two input parameters required.");
	if(nlhs > 1)
		mexErrMsgTxt("One output parameter required.");	
	
	/* Read input parameter dimensions */
	Mrows=mxGetM(prhs[0]);
	Mcols=mxGetN(prhs[0]);
	Nrows=mxGetM(prhs[1]);
	Ncols=mxGetN(prhs[1]);
	
	/* Check input parameter dimension */
	if ((Mcols>1) || (Ncols>1))
		mexErrMsgTxt("Input has to be column vectors.");
	
	/* Create matrix for the return argument D */
	plhs[0]=mxCreateDoubleMatrix(Mrows,1, mxREAL);
    
    /* Assign pointers to each input and output */
	M=mxGetPr(prhs[0]);
	N=mxGetPr(prhs[1]);
	T=mxGetPr(plhs[0]);
    
    int D[Mrows];
    
    for(idx=0; idx < Mrows; idx++) {
        D[idx]=0;
        for (jdx=0; jdx < Nrows; jdx++) {
            if(N[jdx]<M[idx]){
                D[idx]= D[idx] + 1;
            }
        }
        D[idx]= D[idx] + 1;
        //printf("D[%d] = %d\n", idx, D[idx]);
        
        T[idx]=D[idx];
    }    
}
