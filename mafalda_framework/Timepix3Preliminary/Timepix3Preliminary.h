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
 * Created automatically with MAFalda (Mon Jun  9 15:22:48 CDT 2014)
 *
 * An example of how to use the clustering results
 * for further processing.
 */

#ifndef __Timepix3Preliminary_h
#define __Timepix3Preliminary_h

#include <map>

#include "MPXAlgo/MediPixAlgo.h"
#include "CalibrationLoader/CalibrationLoader.h"
#include "MPXAlgo/Highlighter.h"

// This algorithm is using an object put in the StoreGate
//  by BlobsFinder. I need to know the object.
#include "BlobsFinder/BlobsFinder.h"

class Timepix3Preliminary : public MediPixAlgo , public CalibrationLoader {

public:

  Timepix3Preliminary();
  virtual ~Timepix3Preliminary() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

  void Set_selSizeX(int sx) { m_selSizeX = sx;};
  void Set_selSizeY(int sy) { m_selSizeY = sy;};

private:

  AllBlobsContainer * m_aB;
  Int_t m_minNPixels;

  // for output
  vector<int> m_clusterTOT;
  vector<double> m_clusterEnergy;
  vector<double> m_driftingTime1;
  vector<int> m_clusterSizeX;
  vector<int> m_clusterSizeY;

  map<int, vector<double> > m_energyAverageX;
  int m_energyAverageXCntr;
  map<int, vector<double> > m_energyAverageY;
  int m_energyAverageYCntr;

  // Visualization
  Highlighter * m_highLight;
  vector<TCanvas *> m_extraCanvas;
  vector<TGraph2D *> m_graph2DVector;

  // Selection
  int m_minInner;
  int m_minDriftLimit;
  int m_maxDriftLimit;
  int m_selSizeX;
  int m_selSizeY;

  ClassDef(Timepix3Preliminary, 1)
};

#endif
