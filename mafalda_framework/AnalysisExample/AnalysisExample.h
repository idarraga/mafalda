/*
 * 	Copyright 2008 John Idarraga
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

#ifndef AnalysisExample_h
#define AnalysisExample_h

#include "MPXAlgo/MediPixAlgo.h"
#include "MPXAlgo/Highlighter.h"

class AnalysisExample : public MediPixAlgo {

public:

  AnalysisExample();
  virtual ~AnalysisExample() { };

  // You have to implement Init(), Execute() and Finalize()
  //  if you inherit from MediPixAlgo
  void Init();
  void Execute();
  void Finalize();

  Bool_t ThereIsSomething(Int_t, Int_t);

private:

  Int_t m_totalCharge;
  Int_t m_totalHits;
  Int_t m_xcoor;
  Int_t m_ycoor;

  Highlighter * m_arrow;

  ClassDef(AnalysisExample, 1)
};

#endif
