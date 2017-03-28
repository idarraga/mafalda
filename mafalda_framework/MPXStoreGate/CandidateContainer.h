/*
 * 	Copyright 2016 John Idarraga
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

#ifndef CandidateContainer_h
#define CandidateContainer_h

#include <vector>

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>

#include "AnalysisCore/AnalysisCore_defs.h"

using namespace std;

class MediPixAlgo;

/** 
 * Class holding the information about each candidate 
 *
 */

class CandidateContainer : public TObject {

 public:
  CandidateContainer();
  CandidateContainer(MediPixAlgo *);
  CandidateContainer(MediPixAlgo *, MPXDefs::SpecialObjs);
  virtual ~CandidateContainer(){};

  void SetAuthor(TString);
  void SetDSN(TString);
  void SetnFrame(Int_t);
  void SetLikelihood(Double_t);
  
  inline TString GetAuthor(){
    return author;
  };
  inline TString GetDSN(){
    return DSN;
  };
  inline Int_t GetnFrame(){
    return nFrame;
  };
  inline Double_t GetLikelihood(Int_t pos){
    return Likelihood[pos];
  };
  inline Int_t GetFrameXSize(){return m_xFrameSize;};
  inline Int_t GetFrameYSize(){return m_yFrameSize;};
  inline MPXDefs::SpecialObjs GetType(){return m_type;};

 private:
  TString author;
  TString DSN;
  Int_t m_xFrameSize;
  Int_t m_yFrameSize;
  Int_t nFrame;
  vector<Double_t> Likelihood;

  MPXDefs::SpecialObjs m_type; // for special objects
 
  ClassDef(CandidateContainer, 1)
};

#endif
