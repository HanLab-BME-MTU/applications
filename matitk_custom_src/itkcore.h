/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Module:    $RCSfile: itkcore.h,v $
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef ITKCORE_H
#define ITKCORE_H

#include "matitk_custom.h"
#include "ParameterContainer.h"
#include "seedcontainer.h"
#include "itkImage.h"

#define CONTACT "Vincent Chu <vwchu@sfu.ca>, Ghassan Hamarneh <hamarneh@cs.sfu.ca>"
#define OPCOMMAND "matitk(operationName,[parameters],[inputArray1],[inputArray2],[seed(s)Array],[Image(s)Spacing])\n"
//#define CONTACT "HIDDEN TO PRESERVE ANONYMITY"

#include "itkcore.inl"
#endif
