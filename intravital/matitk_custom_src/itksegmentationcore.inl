/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: itksegmentationcore.inl,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef ITKSEGMENTATIONCORE_TXX
#define ITKSEGMENTATIONCORE_TXX

#include "itksegmentationcore.h"
#include "vnl/vnl_math.h"

#include <sstream>

template <class ITKPIXELTYPE>
void ITKSegClass<ITKPIXELTYPE>::ITKSegmentationEntry(MATITKTemplatedVariables<ITKPIXELTYPE>& GTV)
{
	////////////////////////////////start need to append if additing a new function/////////////////////////
	pixelContainer = GTV.pixelContainer;
	seedsIndex = GTV.seedsIndex;
	pixelContainers = GTV.pixelContainers;
	importFilter[0] = GTV.importFilter[0];
	importFilter[1] = GTV.importFilter[1];

	const FunctionCall operations[]={
			FunctionCall("SOTSU","OtsuThresholdImageFilter",&ITKSegClass<ITKPIXELTYPE>::segmentOtsuThreshold,ITKOTSUTHRESHOLDIMAGEFILTERDESC)
	};


	////////////////////////////////end need to append if additing a new function/////////////////////////
	const int nFcn = sizeof(operations)/sizeof(*operations);
	bool listFunctions=false;

	//dispatch the correct one
	bool found=false;
	for (int i=0;i<nFcn;i++){
    if (!STRICMP(operations[i].OpCode,pstrzOp)) {
			mexPrintf("\n%s is being executed...\n",pstrzOp);
			if (bHelpMode) {
				mexPrintf("\n***********Begin description of %s(%s)***********\n\n", operations[i].OpName,pstrzOp);
				mexPrintf(operations[i].OpDesc);
				mexPrintf("\n***************************End description***************************\n\n");
			}
			mexEvalString("drawnow");
			operations[i].ptrFcn();				
			mexPrintf("%s has completed.\n",pstrzOp);
			found=true;
			break;
		}
	}
	if (!found) 	
	{
		mexPrintf("\nThe following segmentation functions are implemented:\n"); 	
		for (int i=0;i<nFcn;i++){
			mexPrintf("%s: %s\n",operations[i].OpCode,operations[i].OpName);
		}
		//if (pstrzOp[0]=='S') mexErrMsgTxt("Unknown Opcode");
	}

	GTV.pixelContainer = pixelContainer;
	GTV.pixelContainers = pixelContainers;
	GTV.importFilter[0] = importFilter[0];
	GTV.importFilter[1] = importFilter[1];
	return;
}

/**************************************************************************

                 itk::OtsuThresholdImageFilter

**************************************************************************/

template <class ITKPIXELTYPE>
void ITKSegClass<ITKPIXELTYPE>::segmentOtsuThreshold()
{	
	const char* PARAM[]={"NumberOfHistogramBins"};
	const char* SUGGESTVALUE[]={"255"};

	const int nMinParam = 1;
	const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);

	/////////////////////Begin Core Filter Code////////////////////////////////////
	int numHistogramBins=(int) mxGetScalar( paramIterator.getCurrentParam(0) );
	typedef itk::OtsuThresholdImageFilter< InternalImageType, InternalImageType > FilterType;
	typename FilterType::Pointer  filter = FilterType::New();
	filter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());
	filter->SetNumberOfHistogramBins( numHistogramBins );
	filter->Update();
	pixelContainer = filter->GetOutput()->GetPixelContainer(); 
	/////////////////////End Core Filter Code////////////////////////////////////
}

#endif