/*
 * 	Copyright 2012 John Idarraga, Mathieu Benoit
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

#ifndef __TOTCalibrationPreparation_h
#define __TOTCalibrationPreparation_h

#include "MPXAlgo/MediPixAlgo.h"

// This algorithm is using an object put in the StoreGate
//  by BlobsFinder. I need to know the object.
#include "BlobsFinder/BlobsFinder.h"
#include "SimpleClustering/SimpleClustering.h"

#include <TROOT.h>
#include <TH2I.h>
#include <TVector.h>

#include <vector>

using namespace std;

class TOTCalibrationPreparation : public MediPixAlgo {

public:

  TOTCalibrationPreparation();
  virtual ~TOTCalibrationPreparation() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

  void FillStatsPlots();
  void CleanUpVars();
  void PrepareEntriesHistograms();
  bool CheckDACs();

  list< pair < pair<int, int>, int > >::iterator FindHotestPixel(list< pair < pair<int, int>, int > > *);

private:

  // StoreGate access to clustering
  AllBlobsContainer * m_aB;
  SimpleClusterContainer * m_sC;

  int m_sizex;
  int m_sizey;

  // Map for single hit entries
  vector<int> m_SingleHitCoor;
  vector<int> m_SingleHitTOT;

  // Map for double hit entries
  vector<int> m_DoubleHitCoor;
  vector<int> m_DoubleHitTOT;

  // Map for triple hit entries
  //map<pair<int,int>, vector<double> > m_TripleHitMap;
  // Map for triple hit entries
  //map<pair<int,int>, vector<double> > m_QuadHitMap;

  long long m_frameId;

  // Number of entries histo
  TH2I * m_SingleHit_EntriesHisto;
  TH2I * m_DoubleHit_EntriesHisto;
  TH2I * m_TripleHit_EntriesHisto;
  TH2I * m_QuadHit_EntriesHisto;

  // Single,Dougle,Triple,Quad hit TOT global distributions
  TH1F * m_SingleHit_GlobalTOTDist;

  // DACs and clock which need to be consistent over the entire run
  vector<int> m_dacs;
  TVectorF * m_DACsStore;
  TVectorF * m_MetaData;

  ClassDef(TOTCalibrationPreparation, 1)
};

#endif
