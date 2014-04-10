/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: MATITKTemplatedVariables.h,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef MATITKTEMPLATEDVARIABLES_H
#define MATITKTEMPLATEDVARIABLES_H
#include "seedcontainer.h"

template<class ITKPIXELTYPE>
class MATITKTemplatedVariables
{
	#include "typedefs.inl"
	
	public:

		std::vector< ITKPIXELTYPEArray<ITKPIXELTYPE> > pixelContainers;

		typename ImageType::PixelContainerPointer pixelContainer;	

		MATPARAMTYPE* pParameters;	

		typename ImportFilterType::Pointer importFilter[2];	

		SeedContainer<MATSEEDTYPE> seedsIndex;

};

#endif