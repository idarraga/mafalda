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

#ifndef Signals_cxx
#define Signals_cxx

#include <typeinfo>

#include "Signals.h"

ClassImp(Signals)


/**
 * Class for drawing Signalss in the Viewer.
 * I will may be move this somewhere else later. 
 *
 */

Signals::Signals () : CandidateContainer () {

}

Signals::Signals (MediPixAlgo * algo, signal_type) : CandidateContainer(algo, MPXDefs::VIS_SKIP) {

}

#endif
