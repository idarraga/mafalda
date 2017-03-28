 /*
 * 	Copyright 2014 John Idarraga
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

/**
 * Created automatically with MAFalda (Fri Jul  4 10:56:08 CEST 2014)
 *
 * An example of how to use the clustering results
 * for further processing.
 */

#ifndef __Timepix3VolcanoEdge_h
#define __Timepix3VolcanoEdge_h

#include "MPXAlgo/MediPixAlgo.h"
#include "CalibrationLoader/CalibrationLoader.h"

// This algorithm is using an object put in the StoreGate
//  by BlobsFinder. I need to know the object.
#include "BlobsFinder/BlobsFinder.h"
#include "MPXAlgo/Highlighter.h"

#include <map>

#define __max_tot_range 1024

class Timepix3VolcanoEdge : public MediPixAlgo , public CalibrationLoader {

public:

  Timepix3VolcanoEdge();
  virtual ~Timepix3VolcanoEdge() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

private:

  AllBlobsContainer * m_aB;
  Int_t m_minNPixels;
  int m_minInner;

  // for output
  vector<int> m_clusterTOT;
  vector<double> m_clusterEnergy;

  map<int, vector<double> > m_TOTEnergyMap;

  ClassDef(Timepix3VolcanoEdge, 1)
};

#endif
