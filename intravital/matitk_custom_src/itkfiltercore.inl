/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: itkfiltercore.inl,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef ITKFILTERCORE_TXX
#define ITKFILTERCORE_TXX

#include "itkfiltercore.h"

#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkBinaryBallStructuringElement.h"


template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::ITKFilterClassEntry(MATITKTemplatedVariables<ITKPIXELTYPE>& GTV)
{
	////////////////////////////////start need to append if additing a new function/////////////////////////
	pixelContainer = GTV.pixelContainer;
	pixelContainers = GTV.pixelContainers;
	importFilter[0] = GTV.importFilter[0];
	importFilter[1] = GTV.importFilter[1];
	
	const FunctionCall operations[]={
			FunctionCall("FDerivative","DerivativeImageFilter",&ITKFilterClass<ITKPIXELTYPE>::DerivativeImageFilter,ITKDERIVATIVEIMAGEFILTERDESC),
			FunctionCall("FGradientMagnitude","GradientMagnitudeFilter",&ITKFilterClass<ITKPIXELTYPE>::filterGradientMagnitude,ITKGRADIENTMAGNITUDEIMAGEFILTERDESC),
			FunctionCall("FGradientMagnitudeSmoothing","GradientMagnitudeWithSmoothingFilter",&ITKFilterClass<ITKPIXELTYPE>::filterGradientMagnitudeWithSmoothing,ITKGRADIENTMAGNITUDEIMAGEFILTERDESC),
			FunctionCall("FSignedMaurerDistancemap","SignedMaurerDistanceMapImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterSignedMaurerDistanceMap,ITKSIGNEDMAURERDISTANCEMAPIMAGEFILTERDESC),			
			FunctionCall("FAAliasBininary","AntiAliasBinaryImageFilter",&ITKFilterClass<ITKPIXELTYPE>::antiAliasBinaryImageFilter,ITKANTIALIASBINARYIMAGEFILTERDESC),
			FunctionCall("FMedian","MedianImageFilter",&ITKFilterClass<ITKPIXELTYPE>::medianImageFilter,ITKMEDIANIMAGEFILTERDESC),
			FunctionCall("FGradientAnisotropicDiffusion","GradientAnisotropicDiffusionImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterGradientAnisotropicDiffusion,ITKGRADIENTANISOTROPICDIFFUSIONIMAGEFILTERDESC),
      FunctionCall("FRecursiveGaussianSmoothing","itkSmoothRecursiveGaussianImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterRecursiveGaussianSmoothing,ITK_SMOOTHING_RECURSIVE_GAUSSIAN_IMAGE_FILTER_DESC),
      FunctionCall("FLaplacianOfGaussian","itkLaplacianRecursiveGaussianImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterLaplacianOfGaussian,ITK_LAPLACIAN_RECURSIVE_GAUSSIAN_IMAGE_FILTER_DESC),
      FunctionCall("FKittlerMinimumErrorThreshold","MinimumErrorThresholdingImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterKittlerMinimumErrorThreshold,ITKKITTLERILLINGWORTHTHRESHOLDIMAGEFILTERDESC),
      FunctionCall("FOtsuThreshold","OtsuThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterOtsuThreshold,ITKOTSUTHRESHOLDIMAGEFILTERDESC2),
      FunctionCall("FKapurMaximumEntropyThreshold","MaximumEntropyThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterKapurMaximumEntropyThreshold,ITKMAXIMUMENTROPYTHRESHOLDIMAGEFILTERDESC),
      FunctionCall("FLiMinimumCrossEntropyThreshold","MinimumCrossEntropyThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterLiMinimumCrossEntropyThreshold,ITKLITHRESHOLDIMAGEFILTERDESC),
      FunctionCall("FShanbhagFuzzyEntropicThreshold","ShanbhagThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterShanbhagFuzzyEntropicThreshold,ITKSHANBHAGTHRESHOLDIMAGEFILTERDESC),
      FunctionCall("FYenEntropyThreshold","YenThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterYenEntropyThreshold,ITKYENTHRESHOLDIMAGEFILTERDESC),
      FunctionCall("FIsoDataThreshold","IsoDataThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterIsoDataThreshold,ITKISODATATHRESHOLDIMAGEFILTERDESC),
      FunctionCall("FTsaiMomentsPreservingThreshold","MomentsThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterTsaiMomentsPreservingThreshold,ITKMOMENTSTHRESHOLDIMAGEFILTERDESC),
      FunctionCall("FHuangFuzzySimilarityThreshold","HuangThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterHuangFuzzySimilarityThreshold,ITKHUANGTHRESHOLDIMAGEFILTERDESC),
      FunctionCall("FRenyiEntropyThreshold","RenyiEntropyThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterRenyiEntropyThreshold,ITK_RENYI_ENTROPY_THRESHOLD_IMAGE_FILTER_DESC),
      FunctionCall("FTriangleThreshold","TriangleThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterTriangleThreshold,ITK_TRIANGLE_THRESHOLD_IMAGE_FILTER_DESC),
      FunctionCall("FPikazTopologicalStableStateThreshold","ThresholdMaximumConnectedComponentsImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterPikazTopologicalStableStateThreshold,ITK_THRESHOLD_MAXIMUM_CONNECTED_COMPONENT_IMAGE_FILTER_DESC),
      FunctionCall("FPrewittIntermodesThreshold","itkIntermodesThresholdImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterPrewittIntermodesThreshold,ITK_INTERMODES_THRESHOLD_IMAGE_FILTER_DESC),
      FunctionCall("FBinaryThinning","itkBinaryThinningImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterBinaryThinning,ITK_BINARY_THINNING_IMAGE_FILTER_DESC),
      FunctionCall("FBinaryDilation","itkBinaryDilateImageFilter",&ITKFilterClass<ITKPIXELTYPE>::filterBinaryDilation,ITK_BINARY_DILATE_IMAGE_FILTER_DESC)
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
		mexPrintf("\nThe following filtering functions are implemented:\n"); 	
		for (int i=0;i<nFcn;i++){
			mexPrintf("%s: %s\n",operations[i].OpCode,operations[i].OpName);
		}
		//if (pstrzOp[0]=='F') mexErrMsgTxt("Unknown Opcode");
	}

	GTV.pixelContainer = pixelContainer;
	GTV.pixelContainers = pixelContainers;
	GTV.importFilter[0] = importFilter[0];
	GTV.importFilter[1] = importFilter[1];
	
	return;	

}

/**************************************************************************
				itk::BinaryDilateImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterBinaryDilation()
{
	const char* PARAM[]={"Radius"};
	const char* SUGGESTVALUE[]={"1.0 or [2.0, 2.0, 1.0]"};

	const int nMinParam = 1;
	const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);

	if( mxIsEmpty(paramIterator.getCurrentParam(0)) || !mxIsNumeric(paramIterator.getCurrentParam(0)) )
	{
		mexErrMsgTxt( "ERROR: radius parameter should be either a scalar or an array of size 3." );
	}
	
	int numRadiusElements = (unsigned int) mxGetNumberOfElements( paramIterator.getCurrentParam(0) );

	if( !( numRadiusElements == 1 || numRadiusElements == 3 ) )
	{
		mexErrMsgTxt( "ERROR: radius parameter should be either a scalar or an array of size 3." );
	}		
	
	double *pRadiusParam = (double *) mxGetData( paramIterator.getCurrentParam(0) );

	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::BinaryBallStructuringElement<ITKPIXELTYPE,DIMENSION> StructuringElementType;
  StructuringElementType structuringElement;
  InputImageType::SizeType szRadius;

  for( int i = 0; i < DIMENSION; i++ )
  {
    if( numRadiusElements == 1)
      szRadius[i] = *pRadiusParam;
    else
      szRadius[i] = pRadiusParam[i];
  }

	structuringElement.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter< InputImageType, OutputImageType, StructuringElementType>  FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
	pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());
  pFilter->SetKernel(structuringElement);
	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 
	/////////////////////End Core Filter Code////////////////////////////////////
}

/**************************************************************************
				itk::BinaryThinningImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterBinaryThinning()
{
	const char* PARAM[]={""};
	const char* SUGGESTVALUE[]={""};

	const int nMinParam = 0;
	const int nMaxParam = 0;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);

	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::BinaryThinningImageFilter< InputImageType, OutputImageType >  FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
	pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());
	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 
	/////////////////////End Core Filter Code////////////////////////////////////
}


/**************************************************************************
				itk::LaplacianRecursiveGaussianImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterLaplacianOfGaussian()
{
	const char* PARAM[]={"Sigma"};
	const char* SUGGESTVALUE[]={"1.0 (sigma is measured in the units of image spacing)"};

	const int nMinParam = 1;
	const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);

	if( mxIsEmpty(paramIterator.getCurrentParam(0)) || !mxIsNumeric(paramIterator.getCurrentParam(0)) || mxGetNumberOfElements( paramIterator.getCurrentParam(0) ) != 1)
	{
		  mexErrMsgTxt( "ERROR: sigma parameter should be a postivie real number" );
	}
	
	double sigmaParam = (double) mxGetScalar( paramIterator.getCurrentParam(0) );

	/////////////////////Begin Core Filter Code////////////////////////////////////
	typedef itk::LaplacianRecursiveGaussianImageFilter< InputImageType, OutputImageType >  FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
	pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

	pFilter->SetSigma( sigmaParam );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 
	/////////////////////End Core Filter Code////////////////////////////////////
}

/**************************************************************************
				itk::SmoothRecursiveGaussianImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterRecursiveGaussianSmoothing()
{
	const char* PARAM[]={"Sigma"};
	const char* SUGGESTVALUE[]={"1.0 or [2.0, 2.0, 1.0]"};

	const int nMinParam = 1;
	const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);

	if( mxIsEmpty(paramIterator.getCurrentParam(0)) || !mxIsNumeric(paramIterator.getCurrentParam(0)) )
	{
		mexErrMsgTxt( "ERROR: sigma parameter should be either a scalar or an array of size 3." );
	}
	
	int numSigmaElements = (unsigned int) mxGetNumberOfElements( paramIterator.getCurrentParam(0) );

	if( !( numSigmaElements == 1 || numSigmaElements == 3 ) )
	{
		mexErrMsgTxt( "ERROR: radius parameter should be either a scalar or an array of size 3." );
	}		
	
	double *pSigmaParam = (double *) mxGetData( paramIterator.getCurrentParam(0) );

	/////////////////////Begin Core Filter Code////////////////////////////////////
	typedef itk::SmoothingRecursiveGaussianImageFilter< InputImageType, OutputImageType >  FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
	pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  itk::FixedArray<double,3> sigma;

	if( numSigmaElements == 1 )
		sigma.Fill( *pSigmaParam );
	else
	{
		for( int i = 0; i < 3; i++ )
			sigma[i] = pSigmaParam[i];
	}

	pFilter->SetSigmaArray( sigma );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 
	/////////////////////End Core Filter Code////////////////////////////////////
}

/**************************************************************************
				itk::ThresholdMaximumConnectedComponentsImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterPikazTopologicalStableStateThreshold()
{
	const char* PARAM[]={"MinimumObjectSizeInPixels"};
  const char* SUGGESTVALUE[]={"1600 (Depends on the size of the objects in your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int minimumObjectSizeInPixelsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: minimumObjectSizeInPixelsParam parameter should be a positive integer." );
	}
  minimumObjectSizeInPixelsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::ThresholdMaximumConnectedComponentsImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetInsideValue( 1 );
  pFilter->SetOutsideValue( 0 );
  pFilter->SetMinimumObjectSizeInPixels( minimumObjectSizeInPixelsParam );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nNumber of Objects: " << pFilter->GetNumberOfObjects() << std::endl;
  std::cout << "\nComptued Threshold: " << pFilter->GetThresholdValue() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////
}

/**************************************************************************
				itk::IntermodesThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterPrewittIntermodesThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::IntermodesThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////
}


/**************************************************************************
				itk::TriangleThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterTriangleThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::TriangleThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////
}

/**************************************************************************
				itk::RenyiEntropyThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterRenyiEntropyThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::RenyiEntropyThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////
}

/**************************************************************************
				itk::HuangThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterHuangFuzzySimilarityThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::HuangThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////

}

/**************************************************************************
				itk::MomentThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterTsaiMomentsPreservingThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::MomentsThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////

}

/**************************************************************************
				itk::IsoDataThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterIsoDataThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::IsoDataThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////

}

/**************************************************************************
				itk::YenThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterYenEntropyThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::YenThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////

}

/**************************************************************************
				itk::ShanbhagThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterShanbhagFuzzyEntropicThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::ShanbhagThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////

}

/**************************************************************************
				itk::LiThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterLiMinimumCrossEntropyThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::LiThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////

}

/**************************************************************************
				itk::MaximumEntropyThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterKapurMaximumEntropyThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::MaximumEntropyThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////

}

/**************************************************************************
				itk::OtsuThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterOtsuThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::OtsuThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////

}


/**************************************************************************
				itk::KittlerIllingworthThresholdImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterKittlerMinimumErrorThreshold()
{
	const char* PARAM[]={"NumberOfHistogramBins"};
  const char* SUGGESTVALUE[]={"255 (depends in the intensity range of your image)"};
  
  const int nMinParam = 1;
  const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate NumberOfHistogramBins parameter
  unsigned int numHistogramBinsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: NumberOfHistogramBins parameter should be a positive integer." );
	}
  numHistogramBinsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
  
	/////////////////////Begin Core Filter Code////////////////////////////////////
  typedef itk::KittlerIllingworthThresholdImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetNumberOfHistogramBins( numHistogramBinsParam );
  pFilter->SetInsideValue( 0 );
  pFilter->SetOutsideValue( 1 );
  pFilter->SetAutoMinimumMaximum( true );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 

  std::cout << "\nComptued Threshold: " << pFilter->GetThreshold() << std::endl;
	/////////////////////End Core Filter Code////////////////////////////////////

}


/**************************************************************************
				itk::GradientAnisotropicDiffusionImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterGradientAnisotropicDiffusion()
{
	const char* PARAM[]={"ConductanceParamter", "NumberOfIterations","TimeStep"};
	const char* SUGGESTVALUE[]={"1.0 (Typical values are in between 0.5 and 2.0)","2.0","0.125 (Should be at or below (PixelSpacing)/2^ImageDimension)"};

	const int nMinParam = 2;
	const int nMaxParam = 3;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,
													                        nMinParam,nMaxParam);

  // Get and validate conductance parameter
  double conductanceParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(0)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(0)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(0)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: Conductance parameter should be a scalar. See suggested values" );
	}
  conductanceParam = mxGetScalar( paramIterator.getCurrentParam(0) );

  // Get and validate NumberOfIterations parameter
  unsigned int numIterationsParam;
  if( mxIsEmpty(paramIterator.getCurrentParam(1)) || 
      !mxIsNumeric(paramIterator.getCurrentParam(1)) || 
      mxGetNumberOfElements(paramIterator.getCurrentParam(1)) != 1 )
	{
		  mexErrMsgTxt( "ERROR: Number of iterations parameter should be a positive integer. See suggested values" );
	}
  numIterationsParam = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(1) );

  // Get and validate Time Step parameter
  double timeStepParam;

  if( nParametersSupplied >= 3 )
  {
    if( mxIsEmpty(paramIterator.getCurrentParam(2)) || 
        !mxIsNumeric(paramIterator.getCurrentParam(2)) || 
        mxGetNumberOfElements(paramIterator.getCurrentParam(2)) != 1 )
	  {
		    mexErrMsgTxt( "ERROR: time step parameter should be a real number. See suggested values" );
	  }
    timeStepParam = mxGetScalar( paramIterator.getCurrentParam(0) );
  }
  else
  {
    typename InputImageType::SpacingType spacing = importFilter[IMPORTFILTERA]->GetOutput()->GetSpacing();

    double minSpacing = spacing[0];
    for( int i = 1; i < 3 ; i++ )
    {
      if( minSpacing < spacing[i] )
        minSpacing = spacing[i];
    }

    timeStepParam = minSpacing/8.0;
  }

	/////////////////////Begin Core Filter Code////////////////////////////////////
	typedef itk::GradientAnisotropicDiffusionImageFilter<InputImageType,
	OutputImageType> FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
  pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

  pFilter->SetConductanceParameter( conductanceParam );
  pFilter->SetNumberOfIterations( numIterationsParam );
  pFilter->SetTimeStep( timeStepParam );
    
	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 
	/////////////////////End Core Filter Code////////////////////////////////////
}

/**************************************************************************
				itk::MedianImageFilter
**************************************************************************/
template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::medianImageFilter()
{
	const char* PARAM[]={"Radius"};
	const char* SUGGESTVALUE[]={"1.0 or [2.0, 2.0, 1.0]"};

	const int nMinParam = 1;
	const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);

	if( mxIsEmpty(paramIterator.getCurrentParam(0)) || !mxIsNumeric(paramIterator.getCurrentParam(0)) )
	{
		mexErrMsgTxt( "ERROR: radius parameter should be either a scalar or an array of size 3." );
	}
	
	int numRadElements = (unsigned int) mxGetNumberOfElements( paramIterator.getCurrentParam(0) );

	if( !( numRadElements == 1 || numRadElements == 3 ) )
	{
		mexErrMsgTxt( "ERROR: radius parameter should be either a scalar or an array of size 3." );
	}		
	
	double *pRadParam = (double *) mxGetData( paramIterator.getCurrentParam(0) );

	/////////////////////Begin Core Filter Code////////////////////////////////////
	typedef itk::MedianImageFilter< InputImageType, OutputImageType >  FilterType;
	typename FilterType::Pointer pFilter = FilterType::New();
	pFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());

	typename FilterType::InputSizeType radius;

	if( numRadElements == 1 )
		radius.Fill( *pRadParam );
	else
	{
		for( int i = 0; i < 3; i++ )
			radius[i] = pRadParam[i];
	}

	pFilter->SetRadius( radius );

	pFilter->Update();
	pixelContainer = pFilter->GetOutput()->GetPixelContainer(); 
	/////////////////////End Core Filter Code////////////////////////////////////
}

/**************************************************************************
                 itk::DerivativeImageFilter
**************************************************************************/

template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::DerivativeImageFilter()
{	
	const char* PARAM[]={"SETORDER","SETDIRECTION"};
	const char* SUGGESTVALUE[]={"",""};
	const int nMinParam = 2;	
	const int nMaxParam = 2;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);

	unsigned int SETORDER = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
	unsigned int SETDIRECTION = (unsigned int) mxGetScalar( paramIterator.getCurrentParam(1) );
	
	///////////////begin core filter code///////////////////////
	typedef itk::DerivativeImageFilter< InputImageType, OutputImageType >  FilterType;

	typename FilterType::Pointer filter = FilterType::New();
	filter->SetOrder(SETORDER);
	filter->SetDirection( SETDIRECTION );
	filter->SetInput( importFilter[IMPORTFILTERA]->GetOutput() );
	filter->Update();
	pixelContainer = filter->GetOutput()->GetPixelContainer();
	///////////////end core filter code///////////////////////
}

template<>
void ITKFilterClass<unsigned char>::DerivativeImageFilter(){ 
	mexPrintf("This method is not supported with this data type! Try converting to double first.");
} 

/**************************************************************************

itk::GradientMagnitudeImageFilter

**************************************************************************/

template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterGradientMagnitude(){
	const char* PARAM[]={""};
	const char* SUGGESTVALUE[]={""};
	
	const int nMinParam = 0;
	const int nMaxParam = 0;
	
	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);
	///////////////////Begin Core Filter Code////////////////////////////////////
	typedef itk::GradientMagnitudeImageFilter< InternalImageType, InternalImageType >  FilterType;

	typename FilterType::Pointer filter = FilterType::New();
	filter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());
	filter->Update();
	pixelContainer = filter->GetOutput()->GetPixelContainer(); 
	///////////////////End Core Filter Code////////////////////////////////////
}

template<>
void ITKFilterClass<unsigned char>::filterGradientMagnitude(){ 
	mexPrintf("This method is not supported with this data type! Try converting to double first.");
} 

/**************************************************************************

itk::GradientMagnitudeRecursiveGaussianImageFilter.h

**************************************************************************/

template <class ITKPIXELTYPE>
void ITKFilterClass<ITKPIXELTYPE>::filterGradientMagnitudeWithSmoothing(){
	const char* PARAM[]={"Sigma"};
	const char* SUGGESTVALUE[]={""};

	const int nMinParam = 1;
	const int nMaxParam = 1;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);
	/////////////////////Begin Core Filter Code////////////////////////////////////
	double sigma=(double) mxGetScalar( paramIterator.getCurrentParam(0) );
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< InternalImageType, InternalImageType >  GradientFilterType;
	typename GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
	gradientMagnitude->SetInput(importFilter[IMPORTFILTERA]->GetOutput());
	gradientMagnitude->SetSigma(  sigma	 );
	gradientMagnitude->Update();
	pixelContainer = gradientMagnitude->GetOutput()->GetPixelContainer(); 
	/////////////////////End Core Filter Code////////////////////////////////////
}

/**************************itk::SignedMaurerDistanceMapImageFilter***************************************************************/

template <class ITKPIXELTYPE>void ITKFilterClass<ITKPIXELTYPE>::filterSignedMaurerDistanceMap()
{	
	const char* PARAM[]={""};	
	const char* SUGGESTVALUE[]={""};	
	const int nMinParam = 0;	
	const int nMaxParam = 0;	
	
	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);	
	
	///////////////////Begin Core Filter Code////////////////////////////////////	
	
	typedef itk::SignedMaurerDistanceMapImageFilter< InternalImageType, InternalImageType >  DistanceMapFilterType;	
	typename DistanceMapFilterType::Pointer pDistanceMapFilter = DistanceMapFilterType::New();	
	pDistanceMapFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());	
	pDistanceMapFilter->SquaredDistanceOff();	
	pDistanceMapFilter->Update();	
	pixelContainer = pDistanceMapFilter->GetOutput()->GetPixelContainer(); 	
	
	///////////////////End Core Filter Code////////////////////////////////////
}


/**********************************itk::AntiAlisBinaryImageFilter*********************************************************/

template <class ITKPIXELTYPE>void ITKFilterClass<ITKPIXELTYPE>::antiAliasBinaryImageFilter()
{	
	const char* PARAM[]={"MaximumRMSChange","MaximumNumberOfIterations"};	
	const char* SUGGESTVALUE[]={"0.07 (must be < 1.0)","1000"};	
	const int nMinParam = 0;	
	const int nMaxParam = 2;	
	
	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);	
	
	/////////////////////Begin Core Filter Code////////////////////////////////////	
	
	double maximumRMSChange = 0.07;	
	
	if( nParametersSupplied == 1 )	
	{		
		maximumRMSChange = (double) mxGetScalar( paramIterator.getCurrentParam(0) );		
		if( maximumRMSChange > 1.0 )		
		{			
			mexErrMsgTxt( "ERROR --- Parameter: MaximumRMSChange must be < 1.0" );		
		}			
	}	
	
	double maxIters = 1000;	
	if( nParametersSupplied == 2 )	
	{		
		maxIters =(double) mxGetScalar( paramIterator.getCurrentParam(1) );	
	}	
	typedef itk::AntiAliasBinaryImageFilter< InternalImageType, InternalImageType >  AntiAliasFilterType;	
	typename AntiAliasFilterType::Pointer  pAntiAliasFilter = AntiAliasFilterType::New();	
	pAntiAliasFilter->SetInput(importFilter[IMPORTFILTERA]->GetOutput());	
	pAntiAliasFilter->SetMaximumRMSError( maximumRMSChange );	
	pAntiAliasFilter->SetNumberOfIterations( maxIters );	
	pAntiAliasFilter->Update();	
	if(pAntiAliasFilter->GetElapsedIterations() >= maxIters)	
	{		
		std::cout << "Warning: Anti-Aliasing Filter did not converge.";	
	}	
	
	pixelContainer = pAntiAliasFilter->GetOutput()->GetPixelContainer(); 	
	
	/////////////////////End Core Filter Code////////////////////////////////////
}

#endif