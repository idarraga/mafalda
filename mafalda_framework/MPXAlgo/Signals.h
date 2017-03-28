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

#ifndef Signals_h
#define Signals_h

#include <vector>

#include <TROOT.h>
#include <TArrow.h>
#include <TEllipse.h>
#include <TBrowser.h>
#include <TEllipse.h>
#include <TString.h>
#include <TMath.h>
#include <TLine.h>

#include "MPXAlgo/MediPixAlgo.h"
#include "MPXStoreGate/MPXStoreGate.h"
#include "MPXStoreGate/CandidateContainer.h"

#include "OutputMng.h"

using namespace std;

/**
 * Class for drawing Signalss in the Viewer.
 * I will may be move this somewhere else later. 
 *
 */

class MAFTArrow;
class MAFTEllipse;
class MAFTLine;


class Signals : public CandidateContainer {

public:

	typedef enum {
		__SIGNAL_SKIP = 0
	} signal_type;

	Signals();
	Signals(MediPixAlgo *, signal_type);
	~Signals(){};

	OutputMng Log;
	MSG::Endreq endreq;


	ClassDef(Signals, 0)
};


#endif
