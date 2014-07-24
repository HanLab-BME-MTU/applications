#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "mex.h"
#include "matrix.h"

#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )

struct STORMexp{

  double lambda;
  double kon,koff,kb,kth;

  size_t imSize;
  double pxSize;
  
  double *pos;
  size_t nEmi;

  double tAct,tExp;

  size_t nCycles,nRep;
  size_t ngamma;

  int gain;
  double bg;

} STORMexp;

void GaussSpot(double x,double y,double amp,double sigma,
    size_t imSize,double *frame){

  size_t xi,yi;
  size_t ig,jg,i,j;
  size_t wg;

  double xo,yo;
  double xx,yy;
  double r2;

  double *area;

  /* coordinates of nearest pixel */
  xi=(size_t)floor(x);
  yi=(size_t)floor(y);
  /* offset from nearest pixel */
  xo=x-(double)xi;
  yo=y-(double)yi;
  /* size of Gaussian spot grid: should be odd */
  wg=(size_t)MAX(floor(sigma),1);
  wg=8*wg+1;

  area=calloc(wg*wg,sizeof(double));

  for(i=0;i<wg;i++){
    xx=xo-(double)(i*1.0-(double)(wg/2));
    xx=xx/sigma;
    for(j=0;j<wg;j++){
      yy=yo-(double)(j*1.0-(double)(wg/2));
      yy=yy/sigma;

      r2=xx*xx+yy*yy;
      area[j*wg+i]=amp*exp(-r2/2.0);
    }
  }

  /* map spot onto current frame */
  for(i=0;i<wg;i++){
    ig=xi+(i-wg/2)-1;
    if( ig < 0 || ig >= imSize )
      continue;

    for(j=0;j<wg;j++){
      jg=yi+(j-wg/2)-1;
      if( jg < 0 || jg >= imSize )
	continue;
      frame[jg*imSize+ig]+=area[j*wg+i];
    }
  }

  free(area);
}

void addNoise(double *frame,size_t imSize,gsl_rng *RNG){

  size_t i,j;

  double val;

  for(i=0;i<imSize;i++){
    for(j=0;j<imSize;j++){
      val=round(frame[j*imSize+i]);
      frame[j*imSize+i]=(double)gsl_ran_poisson(RNG,val);
    }
  }
}

int MLalgo(struct STORMexp *E,double *stack,int *emStateM){

  /* variables obtained through experiment structure */
  double lambda;
  double kon,koff,kb,kth;

  size_t imSize;
  double pxSize;
  
  double *pos;
  size_t nEmi;

  double tAct,tExp;

  size_t nCycles,nRep;
  size_t ngamma;

  int gain;
  double bg;

  /* further variables */
  double sigmaPSF,sigmaPX;
  double aveAmp;

  size_t w;

  int *emState;

  double *frame;
  double sqrtA;

  double eta,zeta;
  size_t nFrames;
  double x,y;
  double amp;
  size_t xi,yi;

  size_t e,i,j,k,l,n;

  gsl_rng *RNG;

  /* get parameters from experiment structure */
  lambda=E->lambda;
  kon=E->kon;
  koff=E->koff;
  kb=E->kb;
  kth=E->kth;

  imSize=E->imSize;
  pxSize=E->pxSize;

  pos=E->pos;
  nEmi=E->nEmi;

  tAct=E->tAct;
  tExp=E->tExp;

  nCycles=E->nCycles;
  nRep=E->nRep;
  ngamma=E->ngamma;

  gain=E->gain;
  bg=E->bg;

  /* derived variables */
  sigmaPSF=(lambda/2.0)/2.35482;
  sigmaPX=sigmaPSF/pxSize;
  w=(size_t)floor(sigmaPX);

  RNG=gsl_rng_alloc(gsl_rng_mt19937);

  emState=calloc(nEmi,sizeof(int));
  frame=calloc(imSize*imSize,sizeof(double));

  aveAmp=gain*ngamma/(2*M_PI*sigmaPX*sigmaPX);
  sqrtA=sqrt(aveAmp);

  for(e=1;e<nEmi;e++){
    eta=gsl_rng_uniform(RNG);
    if( eta < 0.01 )
      emState[e]=1;
  }

  nFrames=0;
  for(n=0;n<nCycles;n++){

    /* activation pulse */
    for(e=0;e<nEmi;e++){
      if( emState[e] == 0 ){
	eta=gsl_rng_uniform(RNG);
	if( eta < tAct*kon )
	  emState[e]=1;
      }
    }

    for(i=0;i<nRep;i++){

      /* background */
      for(xi=0;xi<imSize;xi++){
	for(yi=0;yi<imSize;yi++)
	  frame[yi*imSize+xi]=bg;
      }

      for(e=0;e<nEmi;e++){

        /* save current fluorescent state in return variable */
	emStateM[nFrames*nEmi+e]=emState[e];

	switch( emState[e] ){
	  case 1:
	    /* generate Gaussian spot at (x,y) with given amplitude */
	    x=E->pos[e];
	    y=E->pos[e+nEmi];
	    amp=aveAmp+sqrtA*gsl_ran_ugaussian(RNG);
            GaussSpot(x,y,amp,sigmaPX,imSize,frame);

	    /* does emitter turn off after imaging? */
	    eta=gsl_rng_uniform(RNG);
	    if( eta < tExp*(koff+kb) ){
	      zeta=gsl_rng_uniform(RNG);
	      if( zeta < koff/(koff+kb) )
		emState[e]=0;
	      else
		emState[e]=2;
	    }
	    break;
	  case 0:
	    eta=gsl_rng_uniform(RNG);
	    if( eta < kth*tExp )
	      emState[e]=1;
	    break;
	  default:
	    break;
	}
      }
      /* add noise and store in stack[] */
      addNoise(frame,imSize,RNG);
      for(xi=0;xi<imSize;xi++){
	for(yi=0;yi<imSize;yi++)
	  stack[nFrames*imSize*imSize+yi*imSize+xi]=frame[yi*imSize+xi];
      }
      nFrames+=1;
    }
  }

  gsl_rng_free(RNG);
  free(emState);
  free(frame);

  return 0;
}

void mexFunction(int nlhs,mxArray **plhs,int nrhs,const mxArray **prhs){

  if( nrhs != 1 || !mxIsStruct(prhs[0]) )
    mexErrMsgTxt("Please provide input as single structure\n");

  if( nlhs != 2 )
    mexErrMsgTxt("Please provide exactly 2 output arguments\n");

  struct STORMexp E;

  int success;
  int dims[3];
  const int *cdims;

  double *stack;
  int *ems;

  mxArray *movie;
  mxArray *fstate;

  E.lambda=mxGetScalar(mxGetField(prhs[0],0,"lambda"));
  E.kon=mxGetScalar(mxGetField(prhs[0],0,"kon"));
  E.koff=mxGetScalar(mxGetField(prhs[0],0,"koff"));
  E.kb=mxGetScalar(mxGetField(prhs[0],0,"kb"));
  E.kth=mxGetScalar(mxGetField(prhs[0],0,"kth"));
  E.imSize=mxGetScalar(mxGetField(prhs[0],0,"imSize"));
  E.pxSize=mxGetScalar(mxGetField(prhs[0],0,"pxSize"));

  E.nEmi=mxGetScalar(mxGetField(prhs[0],0,"nEmi"));
  E.tAct=mxGetScalar(mxGetField(prhs[0],0,"tAct"));
  E.tExp=mxGetScalar(mxGetField(prhs[0],0,"tExp"));
  E.nCycles=mxGetScalar(mxGetField(prhs[0],0,"nCycles"));
  E.nRep=mxGetScalar(mxGetField(prhs[0],0,"nRep"));
  E.ngamma=mxGetScalar(mxGetField(prhs[0],0,"ngamma"));
  E.gain=mxGetScalar(mxGetField(prhs[0],0,"gain"));
  E.bg=mxGetScalar(mxGetField(prhs[0],0,"bg"));

  E.pos=mxGetPr(mxGetField(prhs[0],0,"pos"));

  dims[0]=E.imSize;
  dims[1]=E.imSize;
  dims[2]=E.nCycles*E.nRep;

  cdims=dims;

  movie=mxCreateNumericArray(3,cdims,mxDOUBLE_CLASS,mxREAL);
  fstate=mxCreateNumericMatrix(E.nEmi,E.nCycles*E.nRep,mxINT32_CLASS,mxREAL);

  stack=mxGetPr(movie);
  ems=(int*)mxGetData(fstate);

  /*
  ems[0]=234;

  mexPrintf("lambda= %lf\n",E.lambda);
  mexPrintf("kon= %lf\n",E.kon);
  mexPrintf("koff= %lf\n",E.koff);
  mexPrintf("kb= %lf\n",E.kb);

  mexPrintf("nEmi= %ld\n",E.nEmi);
  mexPrintf("ngamma= %ld\n",E.ngamma);
  mexPrintf("tAct= %lf\n",E.tAct);
  mexPrintf("bg= %lf\n",E.bg);

  mexPrintf("nCycles= %ld\n",E.nCycles);
  mexPrintf("nRep= %ld\n",E.nRep);

  for(i=0;i<10;i++)
    mexPrintf("%10.6lf  %10.6lf\n",E.pos[i],E.pos[i+E.nEmi]);

  for(i=0;i<E.imSize;i++){
    for(j=0;j<E.imSize;j++){
      stack[i*E.imSize+j]=i*1.0+j*1.0;
    }
  }
  */

  /* call subroutine */
  success=MLalgo(&E,stack,ems);

  /* store output in plhs[] */
  plhs[0]=mxDuplicateArray(movie);
  plhs[1]=mxDuplicateArray(fstate);
}

/*
int main(int argc,char **argv){

  struct STORMexp E;

  double *stack;
  double pos[]={64.0,64.0,128.0,64.0,196.0,64.0};

  int *ems;

  size_t i,j;

  E.lambda=672.0;
  E.kon=80.0,
  E.koff=2.0;
  E.kb=0.01;
  E.imSize=256;
  E.pxSize=64.0;

  E.nEmi=3;
  E.tAct=0.01;
  E.tExp=0.1;
  E.nCycles=1;
  E.nRep=10;
  E.ngamma=6000;
  E.gain=4;
  E.bg=100.0;

  E.pos=pos;

  ems=calloc(E.nEmi*E.nCycles*E.nRep,sizeof(int));
  stack=calloc(E.imSize*E.imSize*E.nCycles*E.nRep,sizeof(double));
  
  MLalgo(&E,stack,ems);

  for(i=0;i<E.imSize;i++){
    for(j=0;j<E.imSize;j++)
      printf("%ld  ",(size_t)floor(stack[i*E.imSize+j]));
    printf("\n");
  }

  free(ems);
  free(stack);

  return 0;
}
*/
