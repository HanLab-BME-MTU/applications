#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mex.h"
#include "matrix.h"

long int getNeighborPB(long int site,long int step,long int border){
    
    long int erg;
    
    erg=site+step;
    
    if( erg < 1 )
        erg=erg+border;
    else if( erg > border )
        erg=erg-border;
    
    return erg;
}

long int getNeighborNPB(long int site,long int step,long int border){
    
    long int erg;
    
    erg=site+step;
    
    if( erg < 1 || erg > border )
        erg=-1;
    
    return erg;
}

void mexFunction(int nlhs,mxArray **plhs,int nrhs,const mxArray **prhs){
    
    double *dim;
    mwSize nx,ny;
    mxLogical pbc;
    
    mwSize i,j,k,m;
    mwSize xx,yy,zz;
    
    mxArray *clique;
    clique=mxCreateDoubleMatrix(8,2,mxREAL);
    
    double *nb=mxGetPr(clique);
    
    if( nrhs != 2 )
        mexErrMsgTxt("constructNeighborList takes exactly two arguments\n");
    
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
        mexErrMsgTxt("first argument must be real double array\n");
    
    if( !mxIsLogical( prhs[1] ) )
        mexErrMsgTxt("second argument must be logical\n");
    
    /* get size of grid: nx -> rows, ny -> cols */
    dim=mxGetPr(prhs[0]);
    nx=(mwSize)dim[0];
    ny=(mwSize)dim[1];
    /* periodic boundary conditions? */
    pbc=mxGetScalar(prhs[1]);
    
    if( nlhs != 1 )
        mexErrMsgTxt("function generates exactly one output argument\n");
    plhs[0]=mxCreateCellMatrix(nx,ny);
    
    if( pbc ){
        
        for(i=1;i<=nx;i++){
            for(j=1;j<=ny;j++){
                
                zz=(j-1)*nx+(i-1);
                
                m=0;
                /* row above current bin: xx -> i-1 */
                xx=getNeighborPB(i,-1,nx);
                for(k=-1;k<=1;k++){
                    yy=getNeighborPB(j,k,ny);
                    nb[m]=(double)xx;
                    nb[m+8]=(double)yy;
                    m++;
                }
                /* row below current bin: xx -> i+1 */
                xx=getNeighborPB(i,+1,nx);
                for(k=-1;k<=1;k++){
                    yy=getNeighborPB(j,k,ny);
                    nb[m]=(double)xx;
                    nb[m+8]=(double)yy;
                    m++;
                }
                /* bins left and right from current bin */
                xx=i;
                yy=getNeighborPB(j,-1,ny);
                nb[m]=(double)xx;
                nb[m+8]=(double)yy;
                m++;
                
                yy=getNeighborPB(j,+1,ny);
                nb[m]=(double)xx;
                nb[m+8]=(double)yy;
                m++;
                
                mxSetCell(plhs[0],zz,mxDuplicateArray(clique));
            }
        }
    }
    else{
        
        for(i=1;i<=nx;i++){
            for(j=1;j<=ny;j++){
                
                zz=(j-1)*nx+(i-1);
                
                m=0;
                
                xx=getNeighborNPB(i,-1,nx);
                for(k=-1;k<=1;k++){
                    yy=getNeighborNPB(j,k,ny);
                    nb[m]=xx;
                    nb[m+8]=yy;
                    m++;
                }
                
                xx=getNeighborNPB(i,+1,nx);
                for(k=-1;k<=1;k++){
                    yy=getNeighborNPB(j,k,ny);
                    nb[m]=xx;
                    nb[m+8]=yy;
                    m++;
                }
                
                xx=i;
                yy=getNeighborNPB(j,-1,ny);
                nb[m]=xx;
                nb[m+8]=yy;
                m++;
                
                yy=getNeighborNPB(j,+1,ny);
                nb[m]=xx;
                nb[m+8]=yy;
                m++;
                
                mxSetCell(plhs[0],zz,mxDuplicateArray(clique));
            }
        }
    }
}
