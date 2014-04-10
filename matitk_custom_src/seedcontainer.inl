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

#include <vector>

template <class ITKPIXELTYPE>
class SeedContainer
{
	#include "typedefs.inl"

public:

	typedef typename InternalImageType::IndexType SeedPointType;

	SeedContainer()
	{
		m_seedList.clear();
	}

	void set(ITKPIXELTYPE* pfseeds, mwSize numSeeds, mwSize* volumeDimension)
	{
		m_seedList.clear();

		for(int sid = 0 ; sid < numSeeds ; sid++)
		{
			SeedPointType curSeed;

			for(int dim = 0 ; dim < InternalImageType::ImageDimension ; dim++)
			{
				double inval, adjval;				

				//this is done to correct order, and the fact that matlab starts from 1, and c starts from 0
				inval = pfseeds[dim*numSeeds+sid];
				adjval = inval - 1;

				if( adjval < 0 || adjval >= volumeDimension[dim] )
				{
					mexPrintf( "\nERROR: Invalid seed %d , coord %d" , sid , dim ); 
					mexErrMsgTxt("\nLocation of seed outside volume A image space.\nNote: Matlab follows 1-based index where as ITK follows 0-based index");
				}

				curSeed[dim] = adjval;
			}

			m_seedList.push_back( curSeed );
		}		
	}

	int getNumberOfSeeds()
	{
		return m_seedList.size();
	}
	typename const SeedPointType& getIndex(int i)
	{
		if( i >= m_seedList.size() ) 
			mexErrMsgTxt("ERROR: <SeedContainer::getIndex(i)> - Seed Index Access was out of Bounds");

		return m_seedList[i];
	}

	~SeedContainer()
	{
		m_seedList.clear();
	}

private:

	typedef typename std::vector<SeedPointType> SeedListType;
	typename SeedListType m_seedList;
};