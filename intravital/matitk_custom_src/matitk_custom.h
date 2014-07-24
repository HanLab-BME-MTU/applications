/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: matitk.h,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef MATITK_H
#define MATITK_H

#include "itkImage.h"
#include "matrix.h"
#include "itkImportImageFilter.h"
#include <math.h>
#include "mex.h"
#include <streambuf>
#include <ostream>

#define DIMENSION 3
#define LRESULTINDEX 0
#define ROPINDEX 0
#define RSEEDINDEX 4
#define RINPUTAINDEX 2
#define RINPUTBINDEX 3
#define RPARAMINDEX 1
#define IMPORTFILTERA 0
#define IMPORTFILTERB 1
#define RSPACINGINDEX 5

//#define DOUBLEONLY /*Enable this define allow faster compilation, but only support double image type*/

typedef const mxArray MATPARAMTYPE;
typedef double MATSEEDTYPE;

// Macros for mex conversion


class FunctionCall{
	typedef void (*pt2Function)();
public:
	char* OpCode;
	char* OpName;
	pt2Function ptrFcn;
	const char* OpDesc;
	FunctionCall(char* pstrzCode, char* pstrzName, pt2Function fcn){
		OpCode=pstrzCode;
		OpName=pstrzName;
		ptrFcn=fcn;
		OpDesc="Trivial/Unavailable.";
	}
	FunctionCall(char* pstrzCode, char* pstrzName, pt2Function fcn, const char* pstrzDesc):
	OpDesc(pstrzDesc)
	{
		OpCode=pstrzCode;
		OpName=pstrzName;
		ptrFcn=fcn;		
	}
};

template <class ITKPIXELTYPE>
class ITKPIXELTYPEArray{
	#include "typedefs.inl"
public:
	ITKPIXELTYPE* theArray;
	unsigned long numElem;
	bool needTranspose;
	ITKPIXELTYPEArray(ITKPIXELTYPE* ptrDoubleArray, unsigned long lnumberOfElements, bool bneedTranspose=true){
		theArray=ptrDoubleArray;
		numElem=lnumberOfElements;
		needTranspose = bneedTranspose;
	}
};

class mstream : public std::streambuf 
{
	protected:  
	
		virtual std::streamsize xsputn(const char *s, std::streamsize n)  
		{
			mexPrintf("%.*s",n,s);
			mexEvalString("pause(.001);"); 
			return n;
		}
		
		virtual int overflow(int c = EOF)
		{
		    if (c != EOF) 
			{
		      mexPrintf("%.1s",&c);
			  mexEvalString("pause(.001);"); 
			}

			return 1;
		}
};

mstream mout;

// Time and memory probes -- ITK
#include "itkTimeProbesCollectorBase.h"

#define itkProbesCreate() itk::TimeProbesCollectorBase chronometer
#define itkProbesStart( text ) chronometer.Start( text )
#define itkProbesStop( text )  chronometer.Stop( text )
#define itkProbesReport( stream )  chronometer.Report( stream )

#endif