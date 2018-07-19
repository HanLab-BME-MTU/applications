/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

typedef itk::Image< ITKPIXELTYPE, DIMENSION > ImageType;
typedef ImageType InternalImageType;
//the following two lines are here out of laziness - save time during conversion since example code often use In/OutputImageType
typedef ImageType InputImageType;
typedef ImageType OutputImageType;
typedef itk::ImportImageFilter< ITKPIXELTYPE, DIMENSION >   ImportFilterType;
typedef void (*pt2Function)();