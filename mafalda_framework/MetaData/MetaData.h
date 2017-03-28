/*
 * 	Copyright 2010 John Idarraga
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

#ifndef MetaData_h
#define MetaData_h

#include "MPXAlgo/MediPixAlgo.h"

class MetaData : public MediPixAlgo {

public:

  MetaData();
  virtual ~MetaData() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

private:

  Int_t m_fId;
  Float_t m_acqTime;
  Double_t m_startTime;
  TString m_startTimeS;
  Double_t m_HV;

  ClassDef(MetaData, 1)
};

#endif
