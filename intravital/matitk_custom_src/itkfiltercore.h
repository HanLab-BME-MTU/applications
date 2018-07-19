/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: itkfiltercore.h,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef ITKFILTERCORE_H
#define ITKFILTERCORE_H

#include "matitk_custom.h"
#include "ParameterContainer.h"
#include "MATITKTemplatedVariables.h"
#include "itkcore.h"

/**************************************************************************
                 Custom Added Filters
**************************************************************************/

#include "itkDerivativeImageFilter.h"
extern const char *ITKDERIVATIVEIMAGEFILTERDESC;

#include "itkGradientMagnitudeImageFilter.h"
extern const char *ITKGRADIENTMAGNITUDEIMAGEFILTERDESC;

#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
extern const char *ITKGRADIENTMAGNITUDERECURSIVEGAUSSIANIMAGEFILTERDESC;

#include "itkSignedMaurerDistanceMapImageFilter.h"
extern const char *ITKSIGNEDMAURERDISTANCEMAPIMAGEFILTERDESC;

#include "itkAntiAliasBinaryImageFilter.h"
extern const char *ITKANTIALIASBINARYIMAGEFILTERDESC;

#include "itkMedianImageFilter.h"
extern const char *ITKMEDIANIMAGEFILTERDESC;

#include "itkGradientAnisotropicDiffusionImageFilter.h"
extern const char *ITKGRADIENTANISOTROPICDIFFUSIONIMAGEFILTERDESC;

#include "itkKittlerIllingworthThresholdImageFilter.h"
extern const char *ITKKITTLERILLINGWORTHTHRESHOLDIMAGEFILTERDESC;

#include "itkOtsuThresholdImageFilter.h"
extern const char *ITKOTSUTHRESHOLDIMAGEFILTERDESC2;

#include "itkMaximumEntropyThresholdImageFilter.h"
extern const char *ITKMAXIMUMENTROPYTHRESHOLDIMAGEFILTERDESC;

#include "itkLiThresholdImageFilter.h"
extern const char *ITKLITHRESHOLDIMAGEFILTERDESC;

#include "itkShanbhagThresholdImageFilter.h"
extern const char *ITKSHANBHAGTHRESHOLDIMAGEFILTERDESC;

#include "itkYenThresholdImageFilter.h"
extern const char *ITKYENTHRESHOLDIMAGEFILTERDESC;

#include "itkIsoDataThresholdImageFilter.h"
extern const char *ITKISODATATHRESHOLDIMAGEFILTERDESC;

#include "itkMomentsThresholdImageFilter.h"
extern const char *ITKMOMENTSTHRESHOLDIMAGEFILTERDESC;

#include "itkHuangThresholdImageFilter.h"
extern const char *ITKHUANGTHRESHOLDIMAGEFILTERDESC;

#include "itkRenyiEntropyThresholdImageFilter.h"
extern const char *ITK_RENYI_ENTROPY_THRESHOLD_IMAGE_FILTER_DESC;

#include "itkTriangleThresholdImageFilter.h"
extern const char *ITK_TRIANGLE_THRESHOLD_IMAGE_FILTER_DESC;

#include "itkIntermodesThresholdImageFilter.h"
extern const char *ITK_INTERMODES_THRESHOLD_IMAGE_FILTER_DESC;

#include "itkThresholdMaximumConnectedComponentsImageFilter.h"
extern const char *ITK_THRESHOLD_MAXIMUM_CONNECTED_COMPONENT_IMAGE_FILTER_DESC;

#include "itkSmoothingRecursiveGaussianImageFilter.h"
extern const char *ITK_SMOOTHING_RECURSIVE_GAUSSIAN_IMAGE_FILTER_DESC;

#include "itkLaplacianRecursiveGaussianImageFilter.h"
extern const char *ITK_LAPLACIAN_RECURSIVE_GAUSSIAN_IMAGE_FILTER_DESC;

#include "itkBinaryThinningImageFilter.h"
extern const char *ITK_BINARY_THINNING_IMAGE_FILTER_DESC;

#include "itkBinaryDilateImageFilter.h"
extern const char *ITK_BINARY_DILATE_IMAGE_FILTER_DESC;

/**************************************************************************/

template <class ITKPIXELTYPE>
class ITKFilterClass
{

public:

	#include "typedefs.inl"
	static typename ImageType::PixelContainerPointer pixelContainer;
	static std::vector< ITKPIXELTYPEArray<ITKPIXELTYPE> > pixelContainers;
	static typename ImportFilterType::Pointer importFilter[2];
	const static unsigned int Dimension = DIMENSION;

	static void ITKFilterClassEntry(MATITKTemplatedVariables<ITKPIXELTYPE>& GTV);

	// Add a static member function for each filter here
	static void DerivativeImageFilter();

	static void filterGradientMagnitude();

	static void filterGradientMagnitudeWithSmoothing();

	static void filterSignedMaurerDistanceMap();	

	static void antiAliasBinaryImageFilter();

	static void medianImageFilter();

  static void filterRecursiveGaussianSmoothing();

  static void filterLaplacianOfGaussian();

  // Thresholding 
	static void filterGradientAnisotropicDiffusion();

  static void filterKittlerMinimumErrorThreshold();

  static void filterOtsuThreshold();

  static void filterKapurMaximumEntropyThreshold();

  static void filterLiMinimumCrossEntropyThreshold();

  static void filterShanbhagFuzzyEntropicThreshold();

  static void filterYenEntropyThreshold();

  static void filterIsoDataThreshold();

  static void filterTsaiMomentsPreservingThreshold();

  static void filterHuangFuzzySimilarityThreshold();

  static void filterRenyiEntropyThreshold();

  static void filterTriangleThreshold();

  static void filterPikazTopologicalStableStateThreshold();

  static void filterPrewittIntermodesThreshold();

  // mathematical morphology
  static void filterBinaryThinning();
  static void filterBinaryDilation();

};

template <class ITKPIXELTYPE>
typename ITKFilterClass<ITKPIXELTYPE>::ImageType::PixelContainerPointer ITKFilterClass<ITKPIXELTYPE>::pixelContainer;
template <class ITKPIXELTYPE>
	std::vector< ITKPIXELTYPEArray<ITKPIXELTYPE> > ITKFilterClass<ITKPIXELTYPE>::pixelContainers;
template <class ITKPIXELTYPE>
typename ITKFilterClass<ITKPIXELTYPE>::ImportFilterType::Pointer ITKFilterClass<ITKPIXELTYPE>::importFilter[2];

#include "itkfiltercore.inl"

#endif
