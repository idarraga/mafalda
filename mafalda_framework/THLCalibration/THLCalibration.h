 /*
 * 	Copyright 2013 John Idarraga
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
 * Created automatically with MAFalda (Thu Apr  4 10:01:28 CEST 2013)
 *
 * An example of how to use the clustering results
 * for further processing.
 */

#ifndef __THLCalibration_h
#define __THLCalibration_h

#include "MPXAlgo/MediPixAlgo.h"
#include "CalibrationLoader/CalibrationLoader.h"
#include "SimpleClustering/SimpleClustering.h"

// This algorithm is using an object put in the StoreGate
//  by BlobsFinder. I need to know the object.
#include "BlobsFinder/BlobsFinder.h"

#include <TVector.h>
#include <vector>
#include <map>

#include <TGraphErrors.h>
#include <TRandom1.h>

using namespace std;

class THLCalibration : public MediPixAlgo , public CalibrationLoader {

public:

  THLCalibration();
  virtual ~THLCalibration() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

  bool CheckDACs();

private:

  //SimpleClusterContainer * m_aB;
  //SimpleClusterContainer * m_aB;
  AllBlobsContainer * m_aB;

  int m_frameId;
  int m_frameGlobalId;

  // for output
  vector<int> m_clusterTOT;
  vector<double> m_clusterEnergy;

  // DACs and clock which need to be consistent over the entire run
  vector<int> m_dacs;
  TVectorF * m_DACsStore;

  // data
  map<int, vector<int> > m_THLCountsMap;

  // THL dependency
  TGraphErrors * m_g_THLCounts;

  // final data analysis
  int m_nSmoothings;

  ClassDef(THLCalibration, 1)
};

#endif
