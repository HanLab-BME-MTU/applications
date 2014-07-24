
To solve graph labeling problem

Original available from:

http://www.cs.ubc.ca/~schmidtm/Software/UGM.html

CBL modifications:

- Fixed a major memory bug in \UGM\KPM\max_mult.c - by Uday Kurkure, which is also available separately from: http://www.cs.ubc.ca/~schmidtm/Software/UGM/user/max_mult.zip

- Added a MATLAB wrapper to solve graph labeling problem using LBP or TRBP - by Deepak Chittajallu:
UGM\MRFEnergyMin_BeliefPropagation.m
Usage: [labels ] = MRFEnergyMin_BeliefPropagation( FirstOrderCliquePotential, SecondOrderCliques, SecondOrderCliquePotential, 'max_iters' , 100,'BP_TYPE' , 'LBP' )

- Changed mexAll.m for better folder structure management - by Uday Kurkure

-------

Explanation on the memory bug in max_mult.c

There is one function in the UGM library, “max_mult”, which is used during message passing in LBP optimization. So, more the number of iterations, more number of times this function is called. It is coded in C with matlab mex interface.

The author creates a mxMatrix, say m1, with allocated memory and
then creates another matrix, say m2, and allocates memory to it
then fills m2 with values and
later assigns the pointer from m2 to m1 without freeing the memory of m1 which was anyway redundantly allocated.
MATLAB by itself doesn’t clear the memory (though I think it claims to). Thus, with increasing number of optimization cycles, the amount of RAM used by MATLAB increases without any way to free it even after the script has finished. The only way to free the memory is to close MATLAB.

Original code:

plhs[0]=mxCreateDoubleMatrix(m, n, mxREAL); 
arr3=mxMalloc(m*n*sizeof(double));

… do something important…

mxSetPr(plhs[0], arr3); 

Correction:

plhs[0] = mxCreateNumericMatrix( m*n, 1, mxDOUBLE_CLASS, mxREAL); 
arr3 = (double*)mxGetData(plhs[0]);

… do something important…