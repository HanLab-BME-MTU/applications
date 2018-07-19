/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: itkregistrationcore.inl,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef ITKREGISTRATIONCORE_TXX
#define ITKREGISTRATIONCORE_TXX

#include "itkImageRegionIterator.h"
#include "itkImportImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkCastImageFilter.h"
#include "itkregistrationcore.h"

#include "itkWarpImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkThinPlateSplineKernelTransform.h"
#include "itkNumericTraits.h"
#include "itkResampleImageFilter.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

// Registration methods
#include "itkPointSet.h"
#include "itkImageRegistrationMethod.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkPointSetToPointSetRegistrationMethod.h"

// Optimizers
#include "itkPowellOptimizer.h" 
#include "itkLevenbergMarquardtOptimizer.h"

// Metrics
#include "itkMatchCardinalityImageToImageMetric.h"
#include "itkEuclideanDistancePointMetric.h"

// Transforms
#include "itkEuler3DTransform.h"
#include "itkAffineTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredVersorTransformInitializer.h"

// Interpolators
#include "itkNearestNeighborInterpolateImageFunction.h"

// macros
#define PI 3.14
#define RAD2DEG( rad ) (( rad * 180.0 / PI))
#define DEG2RAD( deg ) (( deg * PI / 180.0))

template <class ITKPIXELTYPE>
void ITKRegClass<ITKPIXELTYPE>::ITKRegistrationEntry(MATITKTemplatedVariables<ITKPIXELTYPE>& GTV)
{
	pixelContainer = GTV.pixelContainer;
	seedsIndex = GTV.seedsIndex;
	pixelContainers = GTV.pixelContainers;
	importFilter[0] = GTV.importFilter[0];
	importFilter[1] = GTV.importFilter[1];

	////////////////////////////////start need to append if adding a new function/////////////////////////
	const FunctionCall operations[]={
		FunctionCall("RD","registerDemon",&ITKRegClass<ITKPIXELTYPE>::registerDemon,ITKDEMONSREGISTRATIONFILTERDESC),
		FunctionCall("RTPS","registerThinPlateSpline",&ITKRegClass<ITKPIXELTYPE>::registerThinPlateSpline,ITKTHINPLATESPLINEKERNELTRANSFORMDESC),
		FunctionCall("RLVREG","registerLabeledVolumes",&ITKRegClass<ITKPIXELTYPE>::registerLabeledVolumes,ITKLABELEDVOLREGISTRATIONDESC),
		FunctionCall("RVPSET","registerVolumesBasedOnPointSets",&ITKRegClass<ITKPIXELTYPE>::registerVolumesBasedOnPointSets,ITKPOINTSETVOLREGISTRATIONDESC)
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
		mexPrintf("\nThe following registration functions are implemented:\n"); 	
		for (int i=0;i<nFcn;i++){
			mexPrintf("%s: %s\n",operations[i].OpCode,operations[i].OpName);
		}
		//if (pstrzOp[0]=='R') mexErrMsgTxt("Unknown Opcode");
	}

	GTV.pixelContainer = pixelContainer;
	GTV.pixelContainers = pixelContainers;
	GTV.importFilter[0] = importFilter[0];
	GTV.importFilter[1] = importFilter[1];
	return;
}

/**************************************************************************

                 itk::DemonsRegistrationFilter

**************************************************************************/

template <class ITKPIXELTYPE>
void ITKRegClass<ITKPIXELTYPE>::registerDemon()
{ 
	const char* PARAM[]={"NumberOfHistogramLevels","NumberOfMatchPoints","DemonNumberofIterations",
		"DemonStandardDeviations"};
	const char* SUGGESTVALUE[]={"1024","7","150","1.0"};
	const int nMinParam = 4;
	const int nMaxParam = 4;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);
	if (emptyImportFilter[IMPORTFILTERB]){
		mexErrMsgTxt("This method requires two image volumes.  Input A should be the fixed image.  Input B should be the moving image.");
	}
	/////////////////////Begin Core Filter Code////////////////////////////////////
	/*
	importFilter[IMPORTFILTERA]->Update();
	importFilter[IMPORTFILTERB]->Update();
	unsigned int NumberOfHistogramLevels=(unsigned int) mxGetScalar( paramIterator.getCurrentParam(0) );
	unsigned int NumberOfMatchPoints=(unsigned int) mxGetScalar( paramIterator.getCurrentParam(1) );
	unsigned int DemonNumberofIterations=(unsigned int) mxGetScalar( paramIterator.getCurrentParam(2) );
	double DemonStandardDeviations= mxGetScalar( paramIterator.getCurrentParam(3) );
	
	typedef itk::HistogramMatchingImageFilter<InternalImageType,InternalImageType> MatchingFilterType;
	typename MatchingFilterType::Pointer matcher = MatchingFilterType::New();
	matcher->SetInput( importFilter[IMPORTFILTERB]->GetOutput());
	matcher->SetReferenceImage( importFilter[IMPORTFILTERA]->GetOutput() );
	matcher->SetNumberOfHistogramLevels( NumberOfHistogramLevels );
	matcher->SetNumberOfMatchPoints( NumberOfMatchPoints );
	matcher->ThresholdAtMeanIntensityOn();
	
	typedef itk::Vector< float, Dimension >    VectorPixelType;
	typedef itk::Image<  VectorPixelType, Dimension > DeformationFieldType;	
	typedef itk::DemonsRegistrationFilter<InternalImageType,InternalImageType,DeformationFieldType>   RegistrationFilterType;
	typename RegistrationFilterType::Pointer filter = RegistrationFilterType::New();
	filter->SetFixedImage( importFilter[IMPORTFILTERA]->GetOutput());
	filter->SetMovingImage( matcher->GetOutput() );
	filter->SetNumberOfIterations( DemonNumberofIterations );
	filter->SetStandardDeviations( DemonStandardDeviations);
	filter->Update();
	
	typedef itk::WarpImageFilter<InternalImageType, InternalImageType,DeformationFieldType  > WarperType;
	typedef itk::LinearInterpolateImageFunction<InternalImageType,double>  InterpolatorType;
	typename WarperType::Pointer warper = WarperType::New();
	typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
	warper->SetInput( importFilter[IMPORTFILTERB]->GetOutput());
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( importFilter[IMPORTFILTERA]->GetOutput()->GetSpacing() );
	warper->SetOutputOrigin( importFilter[IMPORTFILTERA]->GetOutput()->GetOrigin() );
	warper->SetDeformationField( filter->GetOutput() );
	warper->Update();
	
	pixelContainer=warper->GetOutput()->GetPixelContainer(); 
	*/
}
/**************************************************************************

                 itk::ThinPlateSplineRegistrationFilter

**************************************************************************/
template <class ITKPIXELTYPE>
void ITKRegClass<ITKPIXELTYPE>::registerThinPlateSpline(){ 
		const char* PARAM[]={""};
		const char* SUGGESTVALUE[]={""};
		//const int nParam = 0;//sizeof(PARAM)/sizeof(*PARAM);

		const int nMinParam = 4;
		const int nMaxParam = 4;


		//ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);
		if (seedsIndex.getNumberOfSeeds()%2==1 || seedsIndex.getNumberOfSeeds()==0) 
			mexErrMsgTxt("This method requires landmarks.  Each landmark should be 3-dimensional, and there should be even number of landmarks (source->target)");
		/*if (emptyImportFilter[IMPORTFILTERB]){
			mexErrMsgTxt("This method requires two image volumes.  Input A should be the source image.  Input B should be a input A's target result.  ");
		}*/
		/////////////////////Begin Core Filter Code////////////////////////////////////
		const double epsilon = 1e-10;
		typedef itk::ThinPlateSplineKernelTransform<ITKPIXELTYPE, 3>   TPSTransform3DType;
		typedef typename TPSTransform3DType::InputPointType PointType3D;
		typedef typename TPSTransform3DType::PointsIterator Points3DIteratorType;
		typedef   itk::LinearInterpolateImageFunction< 
                       InputImageType, double >  InterpolatorType;
		PointType3D sourcePoint3D;
		PointType3D targetPoint3D;
		typename TPSTransform3DType::Pointer tps3D = TPSTransform3DType::New();

		typedef   typename TPSTransform3DType::PointSetType      PointSetType;
		typedef   typename PointSetType::PointIdentifier  PointIdType;
		PointIdType id = itk::NumericTraits< PointIdType >::Zero;
		typename PointSetType::Pointer sourceLandMarks = PointSetType::New();
		typename PointSetType::Pointer targetLandMarks = PointSetType::New();
		typename PointSetType::PointsContainer::Pointer sourceLandMarkContainer = 
                                   sourceLandMarks->GetPoints();
		typename PointSetType::PointsContainer::Pointer targetLandMarkContainer = 
                                   targetLandMarks->GetPoints();

		/*tps3D->GetTargetLandmarks()->GetPoints()->Reserve( seedsIndex.getNumberOfSeeds()/2 );
		tps3D->GetSourceLandmarks()->GetPoints()->Reserve( seedsIndex.getNumberOfSeeds()/2);
		// Create landmark sets
		Points3DIteratorType tps3Ds = tps3D->GetSourceLandmarks()->GetPoints()->Begin();
		Points3DIteratorType tps3Dt = tps3D->GetTargetLandmarks()->GetPoints()->Begin();
		Points3DIteratorType tps3DsEnd  = tps3D->GetSourceLandmarks()->GetPoints()->End();*/
		//double x1, y1, z1,x2, y2, z2;
		for (int i=0; i< seedsIndex.getNumberOfSeeds(); i+=2){
			sourcePoint3D[0] = seedsIndex.getIndex(i)[0];
			sourcePoint3D[1] = seedsIndex.getIndex(i)[1];
			sourcePoint3D[2] = seedsIndex.getIndex(i)[2];
			//tps3Ds.Value() = sourcePoint3D;
			targetPoint3D[0] = seedsIndex.getIndex(i+1)[0];;
			targetPoint3D[1] = seedsIndex.getIndex(i+1)[1];;
			targetPoint3D[2] = seedsIndex.getIndex(i+1)[2];;
			/*tps3Dt.Value() = targetPoint3D;
			tps3Ds++;
			tps3Dt++;*/
			sourceLandMarkContainer->InsertElement( id, sourcePoint3D );
			targetLandMarkContainer->InsertElement( id++, targetPoint3D );
		}
		tps3D->SetSourceLandmarks(sourceLandMarks);
		tps3D->SetTargetLandmarks(targetLandMarks);
		tps3D->ComputeWMatrix();
		typedef typename itk::ResampleImageFilter< InternalImageType, InternalImageType >    ResampleImageFilterType;
		typename ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
		typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
		resample->SetInterpolator( interpolator );
		//resample->SetTransform( tps3D );
		/*
		resample->SetSize(  importFilter[IMPORTFILTERB]->GetOutput() ->GetLargestPossibleRegion().GetSize() );
		resample->SetOutputOrigin(   importFilter[IMPORTFILTERB]->GetOutput() ->GetOrigin() );
		resample->SetOutputSpacing(  importFilter[IMPORTFILTERB]->GetOutput() ->GetSpacing() );*/
		importFilter[IMPORTFILTERA]->Update();
		typename InputImageType::SpacingType spacing = importFilter[IMPORTFILTERA]->GetOutput()->GetSpacing();
		typename InputImageType::PointType   origin  = importFilter[IMPORTFILTERA]->GetOutput()->GetOrigin();
		typename InputImageType::RegionType region = importFilter[IMPORTFILTERA]->GetOutput()->GetBufferedRegion();
		typename InputImageType::SizeType   size =  region.GetSize();

		resample->SetOutputSpacing( spacing );
		resample->SetOutputOrigin(  origin  );
		resample->SetSize( size );
		resample->SetTransform( tps3D );
		resample->SetOutputStartIndex(  region.GetIndex() );
		resample->SetInput(  importFilter[IMPORTFILTERA]->GetOutput() );


		//resample->SetDefaultPixelValue( 0 );
		resample->Update();
		pixelContainer = resample->GetOutput()->GetPixelContainer();

	}

/**************************************************************************

				Registration of two labeled volumes

**************************************************************************/

template <typename TRegistration>
class RegistrationInterfaceCommand_RLVREG : public itk::Command 
{
	public:

		typedef  RegistrationInterfaceCommand_RLVREG   Self;
		typedef  itk::Command                   Superclass;
		typedef  itk::SmartPointer<Self>        Pointer;
		itkNewMacro( Self );

	protected:

		RegistrationInterfaceCommand_RLVREG() {}

	public:

		typedef   TRegistration          RegistrationType;
		typedef   RegistrationType *     RegistrationPointer;

        typedef   itk::PowellOptimizer   OptimizerType;
	    typedef   OptimizerType *        OptimizerPointer;

	void Execute(itk::Object * object, const itk::EventObject & event)
	{
		if( !(itk::IterationEvent().CheckEvent( &event )) )
		{
			return;
		}

		RegistrationPointer pRegistration = dynamic_cast<RegistrationPointer>( object );
	    OptimizerPointer pOptimizer = dynamic_cast< OptimizerPointer >( pRegistration->GetOptimizer() );

		RegistrationType::ScheduleType regSchedule = pRegistration->GetFixedImagePyramid()->GetSchedule();
		unsigned long curLevel = pRegistration->GetCurrentLevel();
		
		std::cout << "\n\nMultiResolution Level : " << curLevel << std::endl;		

		// Modify optimizer scales for translation
		OptimizerType::ScalesType optimizerScales = pOptimizer->GetScales();
		int numParameters = optimizerScales.GetSize();

		for( int i = 0; i < 3 ; i++ )
		{
			optimizerScales[numParameters - 3 + i] *= regSchedule[curLevel][i];
		}

		// Modify step length and tolerance
		if( curLevel == 0 )
		{
			m_startStepLength = pOptimizer->GetStepLength();
		}
		else
		{
			pOptimizer->SetStepLength( m_startStepLength / (curLevel + 1) );
		}
	}

	void Execute(const itk::Object * , const itk::EventObject & )
	{ return; }

private:

	double m_startStepLength;	
};

template <class ITKPIXELTYPE>
void ITKRegClass<ITKPIXELTYPE>::registerLabeledVolumes(){ 
	const char* PARAM[]={"SetFunctionConvergenceTolerance" , "SetMaximumNumberOfIterations" , "SetStepLength" , "SetRegistrationSchedule" , "SetLabelsExcluded" };
	const char* SUGGESTVALUE[]={ "0.001" , "200" , "1.0" , "[1 1 1]" , "" };

	const int nMinParam = 0;
	const int nMaxParam = 5;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);

	// Basic typedefs
	const unsigned int Dimension = 3;
	typedef InputImageType FixedImageType;
	typedef	InputImageType MovingImageType;	
	typedef itk::MultiResolutionImageRegistrationMethod<FixedImageType,MovingImageType> RegistrationType;

	// Get parameters
	double convergenceTolerance;
	long maximumNumberOfIterations;
	double stepLength;
	double *prSchedule;
	double numRegLevels;
	double *prLabelsExcluded;
	double numLabelsExcluded;
	vnl_vector<double> vecLabelsExcluded;
	RegistrationType::ScheduleType regSchedule;

	convergenceTolerance = (double) mxGetScalar( paramIterator.getCurrentParam(0) );		
	maximumNumberOfIterations = (unsigned long) mxGetScalar( paramIterator.getCurrentParam(1) );		
	stepLength = (double) mxGetScalar( paramIterator.getCurrentParam(2) );		

	if( mxGetM( paramIterator.getCurrentParam(3) ) < 1 || mxGetN( paramIterator.getCurrentParam(3) ) != 3 )
	{
		mexErrMsgTxt ("ERROR: Invalid Registration schedule - must be of size numLevels x 3");	
	}	

	prSchedule = mxGetPr( paramIterator.getCurrentParam(3) );

	numRegLevels = mxGetM( paramIterator.getCurrentParam(3) );
	regSchedule.SetSize( numRegLevels , 3 );

	int tmp = 0;
	for( int c = 0 ; c < 3 ; c++ )
	{
		for( int r = 0 ; r < numRegLevels ; r++ )
		{
			regSchedule[r][c] = prSchedule[tmp];
			tmp++;
		}
	}

	if( mxGetM( paramIterator.getCurrentParam(4) ) > 1 ) 
	{
		mexErrMsgTxt("ERROR: SetLabelsExcluded - must be 1D array");	
	}

	numLabelsExcluded =  mxGetNumberOfElements( paramIterator.getCurrentParam(4) );

	if( numLabelsExcluded  > 0 ) 
	{
		prLabelsExcluded = mxGetPr( paramIterator.getCurrentParam(4) );
		vecLabelsExcluded.set_size( numLabelsExcluded );
		vecLabelsExcluded.set( prLabelsExcluded );
	}

	std::cout << "\n\nParameters:\n\n" << std::endl;
	std::cout << "Convergence Tolerance : " << convergenceTolerance << std::endl;
	std::cout << "Max Iterations : " << maximumNumberOfIterations << std::endl;
	std::cout << "Step Length : " << stepLength << std::endl;
	std::cout << "Num Reg Levels : " << numRegLevels << std::endl;
	std::cout << "Registration Schedule : " << regSchedule << std::endl;
	std::cout << "Num Labels Excluded : " << numLabelsExcluded << std::endl;
	std::cout << "Labels Excluded : " << vecLabelsExcluded << std::endl;

	// Fixed Image
	FixedImageType::Pointer pFixedImage; 
	FixedImageType::SizeType sizeFixedImage;
	FixedImageType::SpacingType spacingFixedImage;

	importFilter[IMPORTFILTERA]->Update();
	pFixedImage = importFilter[IMPORTFILTERA]->GetOutput();
	sizeFixedImage = pFixedImage->GetLargestPossibleRegion().GetSize();
	spacingFixedImage = pFixedImage->GetSpacing();

	std::cout << "\nFixed Image Spacing : " << pFixedImage->GetSpacing() << std::endl;
	std::cout << "Fixed Image Size : " << pFixedImage->GetLargestPossibleRegion().GetSize() << std::endl;

	// Moving Image
	MovingImageType::Pointer pMovingImage; 	

	importFilter[IMPORTFILTERB]->Update();	

	if( numLabelsExcluded > 0 )
	{
		typedef itk::ImageDuplicator<MovingImageType> MovingImageDuplicatorType;
		MovingImageDuplicatorType::Pointer pDuplicator = MovingImageDuplicatorType::New();

		pDuplicator->SetInputImage( importFilter[IMPORTFILTERB]->GetOutput() );
		pDuplicator->Update();

		pMovingImage = pDuplicator->GetOutput();

		typedef itk::ImageRegionIterator<MovingImageType> MovingImageIteratorType;
		MovingImageIteratorType itMovingImage( pMovingImage, pMovingImage->GetLargestPossibleRegion());

		/** Exclude the labels specified */
		for( itMovingImage.GoToBegin(); !itMovingImage.IsAtEnd(); ++itMovingImage )
		{
			for( int lid = 0 ; lid < numLabelsExcluded ; lid++ )
			{
				if( (int) itMovingImage.Get() == (int) vecLabelsExcluded[lid] )
				{
					itMovingImage.Set( 0 );
				}							
			}
		}          
	}	
	else
	{
		pMovingImage = importFilter[IMPORTFILTERB]->GetOutput();
	}

	std::cout << "\nMoving Image Spacing : " << pMovingImage->GetSpacing() << std::endl;
	std::cout << "Moving Image Size : " << pMovingImage->GetLargestPossibleRegion().GetSize() << std::endl;

	// Perform registration with Rigid Transform first	
	RegistrationType::Pointer pRegistrationMethod = RegistrationType::New();				
	std::cout << "\n\n--------------- Rigid Registration - Start ---------------- \n\n" << std::endl;

	itkProbesCreate();

	// Input images
	pRegistrationMethod->SetFixedImage( pFixedImage );
	pRegistrationMethod->SetMovingImage( pMovingImage );
	pRegistrationMethod->SetFixedImageRegion( pFixedImage->GetLargestPossibleRegion() );				

	// Metric
	typedef itk::MatchCardinalityImageToImageMetric<FixedImageType,MovingImageType> MetricType;
	MetricType::Pointer pMetric = MetricType::New();

	pMetric->MeasureMatchesOff();

	pRegistrationMethod->SetMetric( pMetric );

	// Interpolator 
	typedef itk::NearestNeighborInterpolateImageFunction< MovingImageType, double> InterpolatorType;
	InterpolatorType::Pointer pInterpolator = InterpolatorType::New();

	pRegistrationMethod->SetInterpolator( pInterpolator );

	// transform
	typedef itk::VersorRigid3DTransform<double> RigidTransformType;		
	RigidTransformType::Pointer pRigidTransform = RigidTransformType::New();	

	// center and rotation
	typedef itk::CenteredTransformInitializer<RigidTransformType,FixedImageType,MovingImageType> TransformInitializerType;			
	TransformInitializerType::Pointer pTransformInitializer = TransformInitializerType::New();

	pTransformInitializer->SetFixedImage( pFixedImage );
	pTransformInitializer->SetMovingImage( pMovingImage );
	pTransformInitializer->SetTransform( pRigidTransform );
	pTransformInitializer->MomentsOn();
	pTransformInitializer->InitializeTransform();	

	std::cout << "\nTransform Initializer" << std::endl;
	std::cout << "\n\tCenter : " << pRigidTransform->GetCenter() << std::endl;
	std::cout << "\n\tTranslation : " << pRigidTransform->GetTranslation() << std::endl;

	// rotation
	typedef RigidTransformType::VersorType VersorType;
	typedef VersorType::VectorType VectorType;

	VectorType rot_axis;

	rot_axis[0] = 0.0;
	rot_axis[1] = 0.0;
	rot_axis[2] = 1.0;

	pRigidTransform->SetRotation( rot_axis , DEG2RAD(1.0) );

	pRegistrationMethod->SetInitialTransformParameters( pRigidTransform->GetParameters() );
	pRegistrationMethod->SetTransform( pRigidTransform );	

	// Optimizer 
	typedef itk::PowellOptimizer OptimizerType;
	OptimizerType::Pointer pOptimizer = OptimizerType::New();

	// parameter scales
	typedef OptimizerType::ScalesType OptimizerScalesType;
	OptimizerScalesType optimizerScalesRigid = OptimizerScalesType( pRigidTransform->GetNumberOfParameters() );

	// versor
	optimizerScalesRigid[0] = 1.0 / cos( DEG2RAD( 89.0 ) );
	optimizerScalesRigid[1] = 1.0 / cos( DEG2RAD( 89.0 ) );
	optimizerScalesRigid[2] = 1.0 / cos( DEG2RAD( 1.0 ) );

	// translation
	optimizerScalesRigid[3] = 1.0 / ( 0.01 * spacingFixedImage[0] * sizeFixedImage[0] );
	optimizerScalesRigid[4] = 1.0 / ( 0.01 * spacingFixedImage[1] * sizeFixedImage[1] );
	optimizerScalesRigid[5] = 1.0 / ( 0.01 * spacingFixedImage[2] * sizeFixedImage[2] );

	pOptimizer->SetScales( optimizerScalesRigid );

	pOptimizer->SetMaximize( false );
	pOptimizer->SetStepLength( stepLength );
	pOptimizer->SetStepTolerance( 1e-4 );
	pOptimizer->SetValueTolerance( convergenceTolerance );
	pOptimizer->SetMaximumIteration( maximumNumberOfIterations );

	pOptimizer->SetInitialPosition( pRigidTransform->GetParameters() );
	pRegistrationMethod->SetOptimizer( pOptimizer );		

	// Observer
	typedef CommandIterationUpdate<OptimizerType> ObserverType;
	ObserverType::Pointer pObserver = ObserverType::New();

	pOptimizer->AddObserver( itk::IterationEvent() , pObserver);

	// Image Pyramids
	typedef itk::MultiResolutionPyramidImageFilter<FixedImageType, FixedImageType> FixedImagePyramidType;
	typedef itk::MultiResolutionPyramidImageFilter<MovingImageType, MovingImageType> MovingImagePyramidType;

	FixedImagePyramidType::Pointer pFixedImagePyramid = FixedImagePyramidType::New();
	MovingImagePyramidType::Pointer pMovingImagePyramid = MovingImagePyramidType::New();

	try
	{	
		pRegistrationMethod->SetFixedImagePyramid( pFixedImagePyramid );
		pRegistrationMethod->SetMovingImagePyramid( pMovingImagePyramid );
		pRegistrationMethod->SetSchedules( regSchedule , regSchedule );	
	}
	catch(itk::ExceptionObject & err)
	{
		std::cout << "\nERROR-SCHEDULES : " << err << std::endl;
	}

	// Registration UI
	typedef RegistrationInterfaceCommand_RLVREG<RegistrationType> RegUIType;
	RegUIType::Pointer pRegUI = RegUIType::New();

	pRegistrationMethod->AddObserver( itk::IterationEvent(), pRegUI );

	// Start Registration
	try
	{	
		itkProbesStart( "Rigid Registration Probe" );

		RegistrationType::ParametersType initialParameters;
		initialParameters = pRigidTransform->GetParameters();
		std::cout << "\nInitial : " << initialParameters << std::endl;

		pRegistrationMethod->StartRegistration();

		itkProbesStop( "Rigid Registration Probe" );
		itkProbesReport( std::cout );
	}
	catch(itk::ExceptionObject & err)
	{
		std::cout << "\nERROR-REGISTRATION : " << err << std::endl;
	}		

	std::cout << "\n\n--------------- Rigid Registration - End ---------------- \n\n" << std::endl;

	// Perform Affine Registration
	std::cout << "\n\n--------------- Affine Registration - Start ---------------- \n\n" << std::endl;

	// Setup Affine Transform
	typedef itk::AffineTransform<double,Dimension> AffineTransformType;
	AffineTransformType::Pointer pAffineTransform = AffineTransformType::New();

	pAffineTransform->SetCenter( pRigidTransform->GetCenter() );
	pAffineTransform->SetTranslation( pRigidTransform->GetTranslation() );
	pAffineTransform->SetMatrix( pRigidTransform->GetMatrix() );

	pRegistrationMethod->SetTransform( pAffineTransform );
	pRegistrationMethod->SetInitialTransformParameters( pAffineTransform->GetParameters() );

	// Setup Optimizer
	pOptimizer->SetInitialPosition( pAffineTransform->GetParameters() );
	pOptimizer->SetStepLength( stepLength );

	// Scales
	typedef OptimizerType::ScalesType OptimizerScalesType;
	OptimizerScalesType optimizerScalesAffine = OptimizerScalesType( pAffineTransform->GetNumberOfParameters() );

	// rotation, scaling, skew
	optimizerScalesAffine[0] = 1.0;
	optimizerScalesAffine[1] = 1.0;
	optimizerScalesAffine[2] = 1.0;
	optimizerScalesAffine[3] = 1.0;
	optimizerScalesAffine[4] = 1.0;
	optimizerScalesAffine[5] = 1.0;
	optimizerScalesAffine[6] = 1.0;
	optimizerScalesAffine[7] = 1.0;
	optimizerScalesAffine[8] = 1.0;

	// translation
	optimizerScalesAffine[9] = 1.0 / ( 0.01 * spacingFixedImage[0] * sizeFixedImage[0] );
	optimizerScalesAffine[10] = 1.0 / ( 0.01 * spacingFixedImage[1] * sizeFixedImage[1] );
	optimizerScalesAffine[11] = 1.0 / ( 0.01 * spacingFixedImage[2] * sizeFixedImage[2] );

	pOptimizer->SetScales( optimizerScalesAffine );

	// Start Registration
	try
	{
		itkProbesStart( "Affine Registration Probe" );

		RegistrationType::ParametersType initialParameters;
		initialParameters = pAffineTransform->GetParameters();
		std::cout << "\nInitial : " << initialParameters << std::endl;

		pRegistrationMethod->StartRegistration();

		itkProbesStop( "Affine Registration Probe" );
		itkProbesReport( std::cout );
	}
	catch(itk::ExceptionObject & err)
	{
		std::cout << "\nERROR-REGISTRATION : " << err << std::endl;
	}

	std::cout << "\n\n--------------- Affine Registration - End ---------------- \n\n" << std::endl;

	// Get Final Parameters
	RegistrationType::ParametersType finalParameters;

	finalParameters = pRegistrationMethod->GetLastTransformParameters();

	std::cout << "Final: " << pMetric->GetValue( finalParameters ) << " --- " << finalParameters << std::endl;

	// Build final transform
	AffineTransformType::Pointer pFinalTransform;

	pFinalTransform = pAffineTransform;
	pFinalTransform->SetParameters( finalParameters );		

	// transform the moving image using the obtained transform parameters
	typedef itk::ResampleImageFilter<MovingImageType,FixedImageType> ResampleFilterType;
	ResampleFilterType::Pointer pResampleImageFilter = ResampleFilterType::New();

	pResampleImageFilter->SetInput( importFilter[IMPORTFILTERB]->GetOutput() );
	pResampleImageFilter->SetTransform( pFinalTransform );	

	pResampleImageFilter->SetSize( pFixedImage->GetLargestPossibleRegion().GetSize() );
	pResampleImageFilter->SetOutputOrigin(  pFixedImage->GetOrigin() );
	pResampleImageFilter->SetOutputSpacing( pFixedImage->GetSpacing() );
	pResampleImageFilter->SetOutputDirection( pFixedImage->GetDirection() );
	pResampleImageFilter->SetDefaultPixelValue( 0 );
	pResampleImageFilter->SetInterpolator( pInterpolator );

	pResampleImageFilter->Update();
	pixelContainer = pResampleImageFilter->GetOutput()->GetPixelContainer();
}

/**************************************************************************

		Registration of two volumes based on point sets using ICP

**************************************************************************/
template <class ITKPIXELTYPE>
void ITKRegClass<ITKPIXELTYPE>::registerVolumesBasedOnPointSets()
{ 
	const char* PARAM[]={"SetFixedPoints" , "SetMovingPoints" , "SetFunctionConvergenceTolerance" , "SetMaximumNumberOfIterations"};
	const char* SUGGESTVALUE[]={ "" , "" , "0.001" , "200" };

	const int nMinParam = 2;
	const int nMaxParam = 4;

	ParameterContainer<ITKPIXELTYPE> paramIterator(PARAM,SUGGESTVALUE,nMinParam,nMaxParam);

	// Basic typedefs
	const unsigned int Dimension = 3;
	typedef InputImageType FixedImageType;
	typedef	InputImageType MovingImageType;	
	typedef double PointCoordType;
	typedef itk::PointSet< PointCoordType, Dimension > PointSetType;
	typedef PointSetType::PointType PointType;
	typedef PointSetType::PointsContainer PointsContainerType;

	// Get parameters
	double convergenceTolerance;
	long maximumNumberOfIterations;	
	unsigned int numFixedPoints;
	unsigned int numMovingPoints;
	double *prFixedPoints;
	double *prMovingPoints;
	vnl_matrix<double> matFixedPoints;
	vnl_matrix<double> matMovingPoints;

	if( !mxIsNumeric( paramIterator.getCurrentParam(0) ) || 
		!mxIsNumeric( paramIterator.getCurrentParam(1) ) ||
		 mxGetN( paramIterator.getCurrentParam(0) ) != Dimension ||
		 mxGetN( paramIterator.getCurrentParam(1) ) != Dimension ||
		 mxGetNumberOfElements( paramIterator.getCurrentParam(0) ) == 0 ||
		 mxGetNumberOfElements( paramIterator.getCurrentParam(1) ) == 0 	
	  )
	{
		mexErrMsgTxt("The input to SetFixedPoints() and SetMovingPoints() should be numeric arrays of size Nx3 where N > 0.\nThe fixed and moving points must be specified in respective image space coordinates.");
	}

	numFixedPoints = mxGetM( paramIterator.getCurrentParam(0) );
	prFixedPoints = mxGetPr( paramIterator.getCurrentParam(0) );
	matFixedPoints = vnl_matrix<PointCoordType>( prFixedPoints , Dimension , numFixedPoints );
	matFixedPoints = matFixedPoints.transpose();

	numMovingPoints = mxGetM( paramIterator.getCurrentParam(1) );
	prMovingPoints = mxGetPr( paramIterator.getCurrentParam(1) );
	matMovingPoints = vnl_matrix<PointCoordType>( prMovingPoints , Dimension , numMovingPoints );
	matMovingPoints = matMovingPoints.transpose();

	convergenceTolerance = (double) mxGetScalar( paramIterator.getCurrentParam(2) );		
	maximumNumberOfIterations = (unsigned long) mxGetScalar( paramIterator.getCurrentParam(3) );		

	std::cout << "\n\nParameters:\n\n" << std::endl;
	std::cout << "Convergence Tolerance : " << convergenceTolerance << std::endl;
	std::cout << "Max Iterations : " << maximumNumberOfIterations << std::endl;
	std::cout << "Num Fixed Points: " << numFixedPoints << std::endl;
	std::cout << "Num Moving Points: " << numMovingPoints << std::endl;

	// Fixed Image
	FixedImageType::Pointer pFixedImage; 
	FixedImageType::SizeType sizeFixedImage;
	FixedImageType::SpacingType spacingFixedImage;

	importFilter[IMPORTFILTERA]->Update();
	pFixedImage = importFilter[IMPORTFILTERA]->GetOutput();
	sizeFixedImage = pFixedImage->GetLargestPossibleRegion().GetSize();
	spacingFixedImage = pFixedImage->GetSpacing();

	std::cout << "\nFixed Image Spacing : " << pFixedImage->GetSpacing() << std::endl;
	std::cout << "Fixed Image Size : " << pFixedImage->GetLargestPossibleRegion().GetSize() << std::endl;

	// Moving Image
	MovingImageType::Pointer pMovingImage; 	

	importFilter[IMPORTFILTERB]->Update();	

	pMovingImage = importFilter[IMPORTFILTERB]->GetOutput();

	std::cout << "\nMoving Image Spacing : " << pMovingImage->GetSpacing() << std::endl;
	std::cout << "Moving Image Size : " << pMovingImage->GetLargestPossibleRegion().GetSize() << std::endl;

	// Fixed Point Set
	PointSetType::Pointer pFixedPointSet = PointSetType::New(); 
	PointsContainerType::Pointer pFixedPointContainer = PointsContainerType::New();
	vnl_vector_fixed<double,Dimension> centerFixedPointSetPhys;

		centerFixedPointSetPhys.fill(0.0);
		for( int i = 0 ; i < numFixedPoints ; i++ )
		{
			PointType fixedPointPhys;
			FixedImageType::IndexType fixedPointIndex;
			
			for( int dim = 0 ; dim < Dimension ; dim++ )
				fixedPointIndex[dim] = matFixedPoints[i][dim];

			if( !( pFixedImage->GetLargestPossibleRegion().IsInside( fixedPointIndex ) ) )
			{
				std::cout << "ERROR: Fixed Image Landmark Out of Bounds : " << i << " -- " << fixedPointIndex << std::endl;
				mexErrMsgTxt("The input to SetFixedPoints() and SetMovingPoints() should be numeric arrays of size Nx3 where N > 0.\nThe fixed and moving points must be specified in respective image space coordinates.");
			}
			
			pFixedImage->TransformIndexToPhysicalPoint( fixedPointIndex , fixedPointPhys );
			
			for( int dim = 0 ; dim < Dimension ; dim++ )
				centerFixedPointSetPhys[dim] += fixedPointPhys[dim];

			pFixedPointContainer->InsertElement( i , fixedPointPhys );
		}

		centerFixedPointSetPhys /= numFixedPoints;

	pFixedPointSet->SetPoints( pFixedPointContainer );

	// Moving Point Set
	PointSetType::Pointer pMovingPointSet = PointSetType::New(); 
	PointsContainerType::Pointer pMovingPointContainer = PointsContainerType::New();
	vnl_vector_fixed<double,Dimension> centerMovingPointSetPhys;

		centerMovingPointSetPhys.fill(0.0);
		for( int i = 0 ; i < numMovingPoints ; i++ )
		{
			PointType movingPointPhys;
			MovingImageType::IndexType movingPointIndex;

			for( int dim = 0 ; dim < Dimension ; dim++ )
				movingPointIndex[dim] = matMovingPoints[i][dim];
				
			if( !( pMovingImage->GetLargestPossibleRegion().IsInside( movingPointIndex ) ) )
			{
				std::cout << "ERROR: Moving Image Landmark Out of Bounds : " << i << " -- " << movingPointIndex << std::endl;
				mexErrMsgTxt("The input to SetFixedPoints() and SetMovingPoints() should be numeric arrays of size Nx3 where N > 0.\nThe fixed and moving points must be specified in respective image space coordinates.");
			}

			pMovingImage->TransformIndexToPhysicalPoint( movingPointIndex , movingPointPhys );

			for( int dim = 0 ; dim < Dimension ; dim++ )
				centerMovingPointSetPhys[dim] += movingPointPhys[dim];

			pMovingPointContainer->InsertElement( i , movingPointPhys );
		}

		centerMovingPointSetPhys /= numMovingPoints;

	pMovingPointSet->SetPoints( pMovingPointContainer );
	
	// Perform registration with Rigid Transform first	
	typedef itk::PointSetToPointSetRegistrationMethod<PointSetType,PointSetType> RegistrationType;
	RegistrationType::Pointer pRegistrationMethod = RegistrationType::New();		

	std::cout << "\n\n--------------- Rigid Registration - Start ---------------- \n\n" << std::endl;

	itkProbesCreate();

	// Input images
	pRegistrationMethod->SetFixedPointSet( pFixedPointSet );
	pRegistrationMethod->SetMovingPointSet( pMovingPointSet );	

	// Metric
	typedef itk::EuclideanDistancePointMetric<PointSetType,PointSetType> MetricType;
	MetricType::Pointer pMetric = MetricType::New();

	pRegistrationMethod->SetMetric( pMetric );

	// transform
	typedef itk::VersorRigid3DTransform<double> RigidTransformType;		
	RigidTransformType::Pointer pRigidTransform = RigidTransformType::New();	

		// center of rotation + translation
		RigidTransformType::InputPointType rotCenter;	
		RigidTransformType::OutputVectorType initTranslation;

		for( int dim = 0 ; dim < Dimension ; dim++ )
		{
			rotCenter[dim] = centerFixedPointSetPhys[dim];
			initTranslation[dim] = centerMovingPointSetPhys[dim] - centerFixedPointSetPhys[dim];
		}

		pRigidTransform->SetCenter( rotCenter );		
		pRigidTransform->SetTranslation( initTranslation );

		// rotation
		typedef RigidTransformType::VersorType VersorType;
		typedef VersorType::VectorType VectorType;

		VectorType rot_axis;

		rot_axis[0] = 0.0;
		rot_axis[1] = 0.0;
		rot_axis[2] = 1.0;

		pRigidTransform->SetRotation( rot_axis , DEG2RAD(1.0) );

	pRegistrationMethod->SetTransform( pRigidTransform );	
	pRegistrationMethod->SetInitialTransformParameters( pRigidTransform->GetParameters() );

	std::cout << "\nTransform Initializer" << std::endl;
	std::cout << "\n\tCenter : " << pRigidTransform->GetCenter() << std::endl;
	std::cout << "\n\tTranslation : " << pRigidTransform->GetTranslation() << std::endl;
	std::cout << "\n\tParameters : " << pRigidTransform->GetParameters() << std::endl;

	// Optimizer 
	typedef itk::LevenbergMarquardtOptimizer OptimizerType;
	OptimizerType::Pointer pOptimizer = OptimizerType::New();

		// parameter scales
		typedef OptimizerType::ScalesType OptimizerScalesType;
		OptimizerScalesType optimizerScalesRigid = OptimizerScalesType( pRigidTransform->GetNumberOfParameters() );

			// versor
			optimizerScalesRigid[0] = 1.0 / cos( DEG2RAD( 89.0 ) );
			optimizerScalesRigid[1] = 1.0 / cos( DEG2RAD( 89.0 ) );
			optimizerScalesRigid[2] = 1.0 / cos( DEG2RAD( 1.0 ) );

			// translation
			optimizerScalesRigid[3] = 1.0 / ( 0.01 * spacingFixedImage[0] * sizeFixedImage[0] );
			optimizerScalesRigid[4] = 1.0 / ( 0.01 * spacingFixedImage[1] * sizeFixedImage[1] );
			optimizerScalesRigid[5] = 1.0 / ( 0.01 * spacingFixedImage[2] * sizeFixedImage[2] );

		pOptimizer->SetScales( optimizerScalesRigid );

	pOptimizer->SetUseCostFunctionGradient( false );
	pOptimizer->SetEpsilonFunction( 1e-6 );
	pOptimizer->SetGradientTolerance( 1e-5 );
	pOptimizer->SetValueTolerance( convergenceTolerance );
	pOptimizer->SetNumberOfIterations( maximumNumberOfIterations );

	pRegistrationMethod->SetOptimizer( pOptimizer );		

	// Observer
	typedef CommandIterationUpdate<OptimizerType> ObserverType;
	ObserverType::Pointer pObserver = ObserverType::New();

	pOptimizer->AddObserver( itk::IterationEvent() , pObserver);

	// Start Registration
	try
	{	
		itkProbesStart( "Rigid Registration Probe" );

		RegistrationType::ParametersType initialParameters;
		initialParameters = pRigidTransform->GetParameters();
		std::cout << "\nInitial : " << initialParameters << std::endl;

		pRegistrationMethod->StartRegistration();

		itkProbesStop( "Rigid Registration Probe" );
		itkProbesReport( std::cout );
	}
	catch(itk::ExceptionObject & err)
	{
		std::cout << "\nERROR-REGISTRATION : " << err << std::endl;
	}		

	std::cout << "\n\n--------------- Rigid Registration - End ---------------- \n\n" << std::endl;

	// Perform Affine Registration
	std::cout << "\n\n--------------- Affine Registration - Start ---------------- \n\n" << std::endl;

	// Setup Affine Transform
	typedef itk::AffineTransform<double,Dimension> AffineTransformType;
	AffineTransformType::Pointer pAffineTransform = AffineTransformType::New();

	pAffineTransform->SetCenter( pRigidTransform->GetCenter() );
	pAffineTransform->SetTranslation( pRigidTransform->GetTranslation() );
	pAffineTransform->SetMatrix( pRigidTransform->GetMatrix() );

	pRegistrationMethod->SetTransform( pAffineTransform );
	pRegistrationMethod->SetInitialTransformParameters( pAffineTransform->GetParameters() );

	// Setup Optimizer

		// Scales
		typedef OptimizerType::ScalesType OptimizerScalesType;
		OptimizerScalesType optimizerScalesAffine = OptimizerScalesType( pAffineTransform->GetNumberOfParameters() );

			// rotation, scaling, skew
			optimizerScalesAffine[0] = 1.0;
			optimizerScalesAffine[1] = 1.0;
			optimizerScalesAffine[2] = 1.0;
			optimizerScalesAffine[3] = 1.0;
			optimizerScalesAffine[4] = 1.0;
			optimizerScalesAffine[5] = 1.0;
			optimizerScalesAffine[6] = 1.0;
			optimizerScalesAffine[7] = 1.0;
			optimizerScalesAffine[8] = 1.0;

			// translation
			optimizerScalesAffine[9] = 1.0 / ( 0.01 * spacingFixedImage[0] * sizeFixedImage[0] );
			optimizerScalesAffine[10] = 1.0 / ( 0.01 * spacingFixedImage[1] * sizeFixedImage[1] );
			optimizerScalesAffine[11] = 1.0 / ( 0.01 * spacingFixedImage[2] * sizeFixedImage[2] );

		pOptimizer->SetScales( optimizerScalesAffine );

	// Start Registration
	try
	{
		itkProbesStart( "Affine Registration Probe" );

		RegistrationType::ParametersType initialParameters;
		initialParameters = pAffineTransform->GetParameters();
		std::cout << "\nInitial : " << initialParameters << std::endl;

		pRegistrationMethod->StartRegistration();

		itkProbesStop( "Affine Registration Probe" );
		itkProbesReport( std::cout );
	}
	catch(itk::ExceptionObject & err)
	{
		std::cout << "\nERROR-REGISTRATION : " << err << std::endl;
	}

	std::cout << "\n\n--------------- Affine Registration - End ---------------- \n\n" << std::endl;

	// Get Final Parameters
	RegistrationType::ParametersType finalParameters;

	finalParameters = pRegistrationMethod->GetLastTransformParameters();

	std::cout << "Final: " << pMetric->GetValue( finalParameters ).mean() << " --- " << finalParameters << std::endl;

	// Build final transform
	AffineTransformType::Pointer pFinalTransform;

	pFinalTransform = pAffineTransform;
	pFinalTransform->SetParameters( finalParameters );		

	// transform the moving image using the obtained transform parameters
	typedef itk::ResampleImageFilter<MovingImageType,FixedImageType> ResampleFilterType;
	ResampleFilterType::Pointer pResampleImageFilter = ResampleFilterType::New();

		// Interpolator 
		typedef itk::NearestNeighborInterpolateImageFunction< MovingImageType, double> InterpolatorType;
		InterpolatorType::Pointer pInterpolator = InterpolatorType::New();

	pResampleImageFilter->SetInput( importFilter[IMPORTFILTERB]->GetOutput() );
	pResampleImageFilter->SetTransform( pFinalTransform );	

	pResampleImageFilter->SetSize( pFixedImage->GetLargestPossibleRegion().GetSize() );
	pResampleImageFilter->SetOutputOrigin(  pFixedImage->GetOrigin() );
	pResampleImageFilter->SetOutputSpacing( pFixedImage->GetSpacing() );
	pResampleImageFilter->SetOutputDirection( pFixedImage->GetDirection() );
	pResampleImageFilter->SetDefaultPixelValue( 0 );
	pResampleImageFilter->SetInterpolator( pInterpolator );

	pResampleImageFilter->Update();
	pixelContainer = pResampleImageFilter->GetOutput()->GetPixelContainer();
}
#endif
