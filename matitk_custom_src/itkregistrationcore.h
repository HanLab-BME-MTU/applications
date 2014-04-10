/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: itkregistrationcore.h,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef ITKREGISTRATIONCORE_H
#define ITKREGISTRATIONCORE_H

#include "matitk_custom.h"
#include "itkcore.h"
#include "seedcontainer.h"
#include "ParameterContainer.h"
#include "MATITKTemplatedVariables.h"
#include "itkCommand.h"

// Add Filter Includes here
// Example: #include "itkDerivativeImageFilter.h"

// Adding extern statements to filter description to eliminate out of scope errors.
extern const char *ITKDEMONSREGISTRATIONFILTERDESC;
extern const char *ITKTHINPLATESPLINEKERNELTRANSFORMDESC;
extern const char *ITKLABELEDVOLREGISTRATIONDESC;
extern const char *ITKPOINTSETVOLREGISTRATIONDESC;

template <class ITKPIXELTYPE>
class ITKRegClass
{

public:

	#include "typedefs.inl"
	static typename ImageType::PixelContainerPointer pixelContainer;

	static SeedContainer<MATSEEDTYPE> seedsIndex;
	typedef SeedContainer<MATSEEDTYPE>::SeedPointType SeedPointType;

	static std::vector< ITKPIXELTYPEArray<ITKPIXELTYPE> > pixelContainers;
	static typename ImportFilterType::Pointer importFilter[2];
	const static unsigned int Dimension = DIMENSION;

	static void ITKRegistrationEntry(MATITKTemplatedVariables<ITKPIXELTYPE>& GTV);

	// Add a static member function for each filter here
	static void registerDemon();
	static void registerThinPlateSpline();
	static void registerLabeledVolumes();
	static void registerVolumesBasedOnPointSets();

};


template <class ITKPIXELTYPE>
typename ITKRegClass<ITKPIXELTYPE>::ImageType::PixelContainerPointer ITKRegClass<ITKPIXELTYPE>::pixelContainer;
template<class ITKPIXELTYPE>
SeedContainer<MATSEEDTYPE> ITKRegClass<ITKPIXELTYPE>::seedsIndex;
template<class ITKPIXELTYPE>
	std::vector< ITKPIXELTYPEArray<ITKPIXELTYPE> > ITKRegClass<ITKPIXELTYPE>::pixelContainers;
template<class ITKPIXELTYPE>
typename ITKRegClass<ITKPIXELTYPE>::ImportFilterType::Pointer ITKRegClass<ITKPIXELTYPE>::importFilter[2];

template<>
void ITKRegClass<unsigned char>::registerThinPlateSpline(){ 
	mexPrintf("This method is not supported with this data type! Try converting to double first.");
}

template<>
void ITKRegClass<float>::registerThinPlateSpline(){ 
	mexPrintf("This method is not supported with this data type! Try converting to double first.");
}

template<>
void ITKRegClass<int>::registerThinPlateSpline(){ 
	mexPrintf("This method is not supported with this data type! Try converting to double first.");
}

template<class TOptimizer>
class CommandIterationUpdate : public itk::Command 
{
public:
	typedef  CommandIterationUpdate   Self;
	typedef  itk::Command             Superclass;
	typedef itk::SmartPointer<Self>  Pointer;
	itkNewMacro( Self );
protected:
	CommandIterationUpdate() 
	{ 
		m_currentIterationIndex = 0; 
	};

public:
	typedef TOptimizer OptimizerType;
	typedef const OptimizerType * OptimizerPointer;

	void Execute(itk::Object *caller, const itk::EventObject & event)
	{
		Execute( (const itk::Object *)caller, event);
	}

	void Execute(const itk::Object * object, const itk::EventObject & event)
	{
		OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
		if( ! itk::IterationEvent().CheckEvent( &event ) )
		{
			return;
		}

		m_currentIterationIndex++;		

		if( strcmp( optimizer->GetNameOfClass() , "LevenbergMarquardtOptimizer" ) == 0 )
		{
			typedef itk::LevenbergMarquardtOptimizer LevenbergOptimizerType;
			typedef const LevenbergOptimizerType * LevenbergOptimizerPointer;

			LevenbergOptimizerPointer pLevenbergOptimizer = dynamic_cast< LevenbergOptimizerPointer >( object );

			std::cout << "Iter " << m_currentIterationIndex << " : " ;	

			LevenbergOptimizerType::MeasureType curValue = pLevenbergOptimizer->GetCachedValue();			
		
			// mexPrintf( "[min, max, mean] = [%f, %f, %f] --- " , curValue.min_value() , curValue.max_value(), curValue.mean() );

			std::cout << curValue.mean() << " --- ";
			std::cout << pLevenbergOptimizer->GetCachedCurrentPosition() << std::endl;			
		}
		else
		{
			std::cout << "Iter " << m_currentIterationIndex << " : " ;		
			std::cout << optimizer->GetValue() << " --- ";
			std::cout << optimizer->GetCurrentPosition() << std::endl;			
		}
	}

private:
	unsigned long m_currentIterationIndex;
};

#include "itkregistrationcore.inl"
#endif
