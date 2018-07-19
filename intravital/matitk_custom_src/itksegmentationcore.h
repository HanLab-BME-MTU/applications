/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: itksegmentationcore.h,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef ITKSEGMENTATIONCORE_H
#define ITKSEGMENTATIONCORE_H

#include "matitk_custom.h"
#include "seedcontainer.h"
#include "ParameterContainer.h"
#include "MATITKTemplatedVariables.h"
#include "itkcore.h"

/**************************************************************************

                 Custom Added Filters

**************************************************************************/

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "itkOtsuThresholdImageFilter.h"
extern const char *ITKOTSUTHRESHOLDIMAGEFILTERDESC;

/**************************************************************************/

template <class ITKPIXELTYPE>
class ITKSegClass
{

public:

	#include "typedefs.inl"
	static typename ImageType::PixelContainerPointer pixelContainer;
	static SeedContainer<MATSEEDTYPE> seedsIndex;
	static std::vector< ITKPIXELTYPEArray<ITKPIXELTYPE> > pixelContainers;
	static typename ImportFilterType::Pointer importFilter[2];
	const static unsigned int Dimension = DIMENSION;

	static void ITKSegmentationEntry(MATITKTemplatedVariables<ITKPIXELTYPE>& GTV);

	// Add a static member function for each filter here
	static void segmentOtsuThreshold();
};

template <class ITKPIXELTYPE>
typename ITKSegClass<ITKPIXELTYPE>::ImageType::PixelContainerPointer ITKSegClass<ITKPIXELTYPE>::pixelContainer;
template <class ITKPIXELTYPE>
SeedContainer<MATSEEDTYPE> ITKSegClass<ITKPIXELTYPE>::seedsIndex;
template <class ITKPIXELTYPE>
	std::vector< ITKPIXELTYPEArray<ITKPIXELTYPE> > ITKSegClass<ITKPIXELTYPE>::pixelContainers;
template <class ITKPIXELTYPE>
typename ITKSegClass<ITKPIXELTYPE>::ImportFilterType::Pointer ITKSegClass<ITKPIXELTYPE>::importFilter[2];

#include "itksegmentationcore.inl"

#endif
