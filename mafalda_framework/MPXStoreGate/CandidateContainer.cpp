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

#ifndef CandidateContainer_cpp
#define CandidateContainer_cpp

#include "CandidateContainer.h"
#include "MPXAlgo/MediPixAlgo.h"

ClassImp(CandidateContainer)

CandidateContainer::CandidateContainer()
  : TObject() {
}

CandidateContainer::CandidateContainer(MediPixAlgo * algo) 
  : TObject() {

  author = algo->GetAlgoName();
  //DSN = algo->GetDataSetNumber();
  nFrame = algo->GetFrameId();
  m_xFrameSize = algo->GetMatrixXdim();
  m_yFrameSize = algo->GetMatrixYdim();
  m_type = MPXDefs::REGULAR; // not an special object

}

CandidateContainer::CandidateContainer(MediPixAlgo * algo,
				       MPXDefs::SpecialObjs type)
  : TObject() {
  
  //CandidateContainer(algo);
  author = algo->GetAlgoName();
  //DSN = algo->GetDataSetNumber();
  nFrame = algo->GetFrameId();
  m_xFrameSize = algo->GetMatrixXdim();
  m_yFrameSize = algo->GetMatrixYdim();
  m_type = type; // some special object, see MPXDefs::SpecialObjs

}

void CandidateContainer::SetAuthor(TString author_i){
  author = author_i;
}

void CandidateContainer::SetDSN(TString DSN){
  DSN=DSN;
}

void CandidateContainer::SetnFrame(Int_t nFrame){
  nFrame=nFrame;
}

void CandidateContainer::SetLikelihood(Double_t element){
  Likelihood.push_back(element);
}


#endif
