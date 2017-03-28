 /*
 * 	Copyright 2011 John Idarraga, Mathieu Benoit
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

#ifndef __EffiencyMap_h
#define __EffiencyMap_h

#include "MPXAlgo/MediPixAlgo.h"

// This algorithm is using an object put in the StoreGate
//  by BlobsFinder. I need to know the object.
#include "BlobsFinder/BlobsFinder.h"

class EffiencyMap : public MediPixAlgo {

public:

  EffiencyMap();
  virtual ~EffiencyMap() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

private:

  AllBlobsContainer * m_aB;
  Int_t m_maxNHits;
  bool m_firstTimeFlag;

  TH2F * m_hitMap;
  TH2F * m_equalizationMap;
  vector<vector<int> > m_perPixelTOT;

  ClassDef(EffiencyMap, 1)
};

#endif
