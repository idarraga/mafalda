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

#ifndef MediPixStoreGate_h
#define MediPixStoreGate_h

#include <map>
#include <iostream>
#include <map>
#include <vector>

#include <TROOT.h>
#include <TString.h>

#include "CandidateContainer.h"
#include "AnalysisCore/AnalysisCore_defs.h"

#include "MPXAlgo/OutputMng.h"

using namespace std;

/**
 * The StoreGate 
 *  
 *   The concept is pretty standard, but here in particular I borrowed
 *   it from the ATLAS software (this one is a totally independent,
 *   made-from-scratch implementation tough).  If an algorithm
 *   produces an object, and you would like to keep it in memory for
 *   the next algo in the chain to grab and read, this piece of code
 *   provides with such interface.
 *
 *   ATTENTION !: The user shouldn't use this class directly.  A proper
 *   interface to it is available through MediPixAlgo from where all
 *   algorithms should inherit.
 * 
 */

class MediPixStoreGate {

 public:
  MediPixStoreGate();
  virtual ~MediPixStoreGate() { };
  Bool_t SaveObject(CandidateContainer *);

  Int_t GetNObjWithType(MPXDefs::SpecialObjs);
  Int_t GetNObjWithTypeLessEqualThan(MPXDefs::SpecialObjs);
  vector<CandidateContainer *> GetObjsSpecial(MPXDefs::SpecialObjs);
  vector<CandidateContainer *> GetObjsSpecialWithAuthor(MPXDefs::SpecialObjs, TString);

  Int_t GetNObjWithAuthor(TString);
  CandidateContainer * GetObjFrom(TString, Int_t);

  Int_t CleanUpAllStoreGate();
  Int_t CleanUpAllStoreGate(MPXDefs::SpecialObjs);
  Int_t CleanUpAllStoreGateExcept(MPXDefs::SpecialObjs);

  //
  void SetOutputLevel(MSG::Level l){Log.OutputLevel = l;};
  
  // log service
  OutputMng Log;
  MSG::Endreq endreq;

 private:
  
  Int_t m_nObjects;
  map<TString, vector<CandidateContainer *> > gateMap;
  map<TString, int > gateMapNumberOfSpecialObjs;

};


#endif
