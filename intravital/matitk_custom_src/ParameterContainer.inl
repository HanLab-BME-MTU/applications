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

#include "ParameterContainer.h"

template <class ITKPIXELTYPE>
class ParameterContainer
{

private:

	const char* const* m_ppstrzParamNames;
	const char* const* m_ppstrzParamSuggestions;
	const int m_nMinParametersExpected;
	const int m_nMaxParametersExpected;

public:

	static MATPARAMTYPE* pParameters;

	ParameterContainer(const char* const* ParameterNames, const char* const* ParameterSuggestions, int nMinParametersExpected , int nMaxParametersExpected )
		:m_ppstrzParamNames(ParameterNames),m_ppstrzParamSuggestions(ParameterSuggestions),
	     m_nMinParametersExpected(nMinParametersExpected),m_nMaxParametersExpected(nMaxParametersExpected)
	{
		if (bHelpMode || nParametersSupplied < m_nMinParametersExpected || nParametersSupplied > m_nMaxParametersExpected){
			promptParameters();
			mexPrintf("Minimum Parameters Expected = %d, Maximum Paramters expected = %d. You supplied %d.\n", m_nMinParametersExpected, m_nMaxParametersExpected, nParametersSupplied);
			mexErrMsgTxt("Correct number of parameters must be supplied.  At least one image volume has to be supplied.");
		}
	}

	void promptParameters()
	{
		mexPrintf("You must supply parameters for this function in an array, with the elements in this order:\n");

		for(int i=0; i < m_nMaxParametersExpected; i++)
		{
			char *strOptional = "";

			if( i >= m_nMinParametersExpected )
				strOptional = "(optional) ";

			if (strcmp(m_ppstrzParamSuggestions[i],"")) 
				mexPrintf("%s%s (which usually has value equal to %s)", strOptional, m_ppstrzParamNames[i], m_ppstrzParamSuggestions[i] );
			else 
				mexPrintf("%s%s", strOptional, m_ppstrzParamNames[i]);

			if( i < m_nMaxParametersExpected-1) 
				mexPrintf(",\n");
		}
		mexPrintf("\n");
	}

	MATPARAMTYPE * getCurrentParam(int i)
	{
		if (i>=nParametersSupplied) 
			mexErrMsgTxt("Out of index exception.  Programming error?!?");

		return mxGetCell( pParameters , i );		
	}

	//~ParameterContainer();
};

template <class ITKPIXELTYPE>
MATPARAMTYPE* ParameterContainer<ITKPIXELTYPE>::pParameters;