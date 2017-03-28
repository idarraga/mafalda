/*
 * 	Copyright 2009 John Idarraga
 *
 * 	This file is part of MAFalda.

    MAFalda is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    MAFalda is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MAFalda.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OutputMng_h
#define OutputMng_h

#include <iostream>
#include "AnalysisCore/AnalysisCore_defs.h"

/**
 * This class manages the output level for MPXAlgos
 *
 */

class OutputMng {

public:
	OutputMng();
	virtual ~OutputMng() { };

	OutputMng & operator<<(MSG::Level);
	OutputMng & operator<<(MSG::Endreq);
	OutputMng & operator<<(OutputMng& (*)(OutputMng&));
	//OutputMng & operator<<(std::ostream& (*) (std::ostream&));
	//OutputMng & operator<<(std::ios& (*)(std::ios&));
	//OutputMng & operator<<(std::ios_base& (*)(std::ios_base&));

	void setAlgoName(TString);
	Bool_t isActive();

	MSG::Level OutputLevel;
	MSG::Level requestedLevel;

	Int_t newOutputLine;
	Int_t newAlgoIndicator;
	TString algoName;

	ClassDef(OutputMng, 1)
};

// non member functions needed

inline void printBlank(TString nameA){

	for(Int_t spc = 0 ; spc < LEFT_MSG_TAB - nameA.Sizeof() ;spc++)
	{
		std::cout << " ";
	}
}

template <typename T>
OutputMng & operator << (OutputMng& lhs, const T& arg)  {

	//std::cout.precision(5);

	if(lhs.isActive())
	{
		if(lhs.newOutputLine)
		{
			if(lhs.newAlgoIndicator)
			{
				std::cout << "["<< lhs.algoName << "] ";
				printBlank(lhs.algoName);
				lhs.newAlgoIndicator = 0;
			}
			else
			{

			}
			lhs.newOutputLine = 0;
		}
		std::cout << arg;
	}
	return lhs;
}


#endif
