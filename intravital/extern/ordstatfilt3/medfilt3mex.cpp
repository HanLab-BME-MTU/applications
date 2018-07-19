#include "mex.h"
#include <vector>
#include <algorithm>

using namespace std;

void exch(std::vector<double> &a, int i, int j)
{
    double temp;
    
    temp = a[i];
    a[i] = a[j];
    a[j] = temp;
}

int partition(std::vector<double> &a, int lbound, int rbound)
{
  int i = lbound;
  int j = rbound + 1;

  while(true)
  {
    while( a[++i] < a[lbound] )
    {
      if(i >= rbound)
        break;
    }

    while( a[--j] > a[lbound] )
    {
      if(j <= lbound)
          break;
    }

    if(i >= j)
      break;

    exch(a, i, j);
  }
  
  exch(a, lbound, j);

  return j;
}

void quickselect(std::vector<double> &a, int lbound, int rbound, int k)
{
  if( rbound <= lbound )
    return;

  int i = partition(a, lbound, rbound);

  if( i < k )
    quickselect(a, i, rbound, k);
  
  if( i > k )
    quickselect(a, lbound, i, k);
}

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, mxArray *prhs[])
{
    
    /* check for proper number of input and output arguments */
    if( nrhs != 2 )
    {
        mexErrMsgIdAndTxt( "MATLAB:medfilt3mex:invalidNumInputs",
                "Two input arguments are required.");
    }
    
    if(nlhs > 1)
    {
        mexErrMsgIdAndTxt( "MATLAB:medfilt3mex:invalidNumOutputs",
                "Too many output arguments.");
    }

    /* check the validity of the input arguments */
    if( mxGetNumberOfDimensions( prhs[0] ) != 3 )
    {
        mexErrMsgIdAndTxt( "MATLAB:medfilt3mex:invalidInputImageDimension",
                "This function works only on 3D images.");
    }
        
    if( mxGetNumberOfElements(prhs[1]) != 3 )
    {
        mexErrMsgIdAndTxt( "MATLAB:medfilt3mex:invalidNeighRadius",
                "The neighborhood radius argument must contain 3 elements.");
    }
    
    /* Get the data */
    const mwSize *pImageSize;
    double *pInputImage, *pFilteredImage;
    double *pNeighRad;
            
    pImageSize = mxGetDimensions( prhs[0] );
    pInputImage = mxGetPr( prhs[0] );    
    pNeighRad = mxGetPr( prhs[1] );
    
    plhs[0] = mxCreateNumericArray(3, pImageSize, mxDOUBLE_CLASS, mxREAL );
    pFilteredImage = mxGetPr( plhs[0] );
    
    /* perform median filtering */
    int numNeighElements = 1;
    
    for(int i = 0; i < 3; i++ )
    {
      numNeighElements = numNeighElements * (2*pNeighRad[i]+1);
    }
    
    std::vector<double> neighElementVec(numNeighElements, 0);

    int sz_r, sz_c, sz_rc;
            
    sz_r = pImageSize[0];
    sz_c = pImageSize[1];
    sz_rc = sz_r * sz_c;

    for(int r = pNeighRad[0]; r < pImageSize[0] - pNeighRad[0]; r++)
    {
        for(int c = pNeighRad[1]; c < pImageSize[1] - pNeighRad[1]; c++)
        {
            for(int z = pNeighRad[2]; z < pImageSize[2] - pNeighRad[2]; z++)
            {
                // get all the elements in the neighborhood
                int i = 0;
                for(int nr = -pNeighRad[0]; nr <= pNeighRad[0]; nr++)
                {
                    for(int nc = -pNeighRad[1]; nc <= pNeighRad[1]; nc++)
                    {
                        for(int nz = -pNeighRad[2]; nz <= pNeighRad[2]; nz++)
                        {
                            int cur_rind, cur_cind, cur_zind;
                            double cur_val;
                            
                            cur_rind = r + nr;
                            cur_cind = c + nc;
                            cur_zind = z + nz;
                            
                            int cur_neighlinind = sz_rc * cur_zind + sz_c * cur_cind + cur_rind;
                            
                            neighElementVec[i++] = pInputImage[ cur_neighlinind ];
                            
                        }
                    }
                }
                
                // get median 
                std::random_shuffle(neighElementVec.begin(), neighElementVec.end());
                
                int k = (numNeighElements-1)/2;
                quickselect(neighElementVec, 0, numNeighElements-1, k);

                int cur_outlinind = sz_rc * z + sz_c * c + r;                
                pFilteredImage[ cur_outlinind ] = neighElementVec[k];
                
            }
        }
    }
    
}