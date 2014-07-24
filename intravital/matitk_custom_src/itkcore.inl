/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: itkcore.inl,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#include "itkcore.h"

#define DEBUG_FILTERING
#define DEBUG_SEGMENTATION
#define DEBUG_REGISTRATION

#ifdef DEBUG_FILTERING
#include "itkfiltercore.h"
#endif
#ifdef DEBUG_SEGMENTATION
#include "itksegmentationcore.h"
#endif
#ifdef DEBUG_REGISTRATION
#include "itkregistrationcore.h"
#endif
#include "MATITKTemplatedVariables.h"


template <class ITKPIXELTYPE>
class ITKCore{
public:
	static void ITKEntry(MATITKTemplatedVariables<ITKPIXELTYPE>& GTV)
	{
		if (pstrzOp[0]=='F' || pstrzOp[0]=='f') {
#ifdef DEBUG_FILTERING
			ITKFilterClass<ITKPIXELTYPE>::ITKFilterClassEntry(GTV);
#endif
		}
		else if (pstrzOp[0]=='S' || pstrzOp[0]=='s') {
#ifdef DEBUG_SEGMENTATION
			ITKSegClass<ITKPIXELTYPE>::ITKSegmentationEntry(GTV);
#endif
		}
		else if (pstrzOp[0]=='R' || pstrzOp[0]=='r') {
#ifdef DEBUG_REGISTRATION
			ITKRegClass<ITKPIXELTYPE>::ITKRegistrationEntry(GTV);
#endif
		}
		else
		{
#ifdef DEBUG_FILTERING
			ITKFilterClass<ITKPIXELTYPE>::ITKFilterClassEntry(GTV);
#endif
#ifdef DEBUG_SEGMENTATION
			ITKSegClass<ITKPIXELTYPE>::ITKSegmentationEntry(GTV);
#endif
#ifdef DEBUG_REGISTRATION
			ITKRegClass<ITKPIXELTYPE>::ITKRegistrationEntry(GTV);
#endif
			aboutMATITK();
			if (pstrzOp[0]!='?')
				mexErrMsgTxt("Unknown Opcode");
		}
		return;
	}
};

