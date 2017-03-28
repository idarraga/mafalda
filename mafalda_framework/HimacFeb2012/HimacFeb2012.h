/*
 * 	Copyright 2011 John Idarraga
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

#ifndef __HimacFeb2012_h
#define __HimacFeb2012_h

#include "BlobsFinder/BlobsFinder.h"
#include "MPXAlgo/Signals.h"
#include "MPXAlgo/MediPixAlgo.h"
#include "CalibrationLoader/CalibrationLoader.h"

#include <TH2I.h>
#include <TH2F.h>
#include <TList.h>
#include <TGraph.h>

#include <vector>
#include <queue>

using namespace std;

// clases
class TGraph2D;
class TCanvas;
class blob;
class AllBlobsContainer;

#define __max_in_nContours_queue         5
#define __steady_sum_in_nContours_queue  5
#define __excess_sum_in_nContours_queue  7
#define __nContour_initState             8
#define __nContour_reached_plateau      10
#define __nContour_reached_crown        20

#define __missing_height 0.9
#define __firstfew_levels  3

// prototypes
double paraboloid_f(double * xx, double * par);

class HimacFeb2012 : public MediPixAlgo , public CalibrationLoader {

public:

  HimacFeb2012();
  virtual ~HimacFeb2012() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

  //void SetCalibrationConfigFile_a(const char * f) { m_a_configFile = f; m_calibWord |= __calib_a_loaded; };
  //void SetCalibrationConfigFile_b(const char * f) { m_b_configFile = f; m_calibWord |= __calib_b_loaded; };
  //void SetCalibrationConfigFile_c(const char * f) { m_c_configFile = f; m_calibWord |= __calib_c_loaded; };
  //void SetCalibrationConfigFile_t(const char * f) { m_t_configFile = f; m_calibWord |= __calib_t_loaded; };
  //void SetMaxFrameWidth(int mw) { m_maxWidth = mw; };
  //void SetMaxFrameHeight(int mh) { m_maxHeight = mh; };
  //void SetCalibClkFactor(double f) { m_clockFactor = f;};
  //void SetCalibClk(double f) { m_timepixClk = f;};

  bool ParticleTagging(blob);
  bool ClusterInBorder(blob);
  int CheckDiscontinuityWithModifiedMatrix(blob bl, map<int,int> totmap, int discontinuityTolerance);
  vector<pair<double, double> > FindLimitPixels_MaxMinXY(blob theBlob);

  void ProcessHIMAC_MC();
  double GetBiggestPerpDistance(double, double, vector<pair<int,int> >);
  double GetSmallestPerpDistance(double, double, vector<pair<int,int> >);
  double GetBiggestPerpDistanceSpecialCaseSlopeZero_Y(double, vector<pair<int,int> >);
  double GetBiggestPerpDistanceSpecialCaseSlopeZero_X(double, vector<pair<int,int> >);
  void GetMaxMin(vector<pair<int,int> >, int *, int *, int *, int *);
  double GetSmallestDistanceTo(pair<int, int>, vector<pair<int,int> >);
  double GetBiggestDistanceTo(pair<int, int>, vector<pair<int,int> >);
  vector<TGraph *> SelectSlicesAboveMean(vector<TGraph *>, double);
  vector<TF1 * > GetKernelDensityFunctions(int, vector<TGraph *>);
  void FindAllCritialPoints(vector<TGraph * >, vector<TF1 *>);

  TString BuildCutGeo1();
  TString BuildCutGeo_Hull();

  void UseMaskInput(bool in){ m_inputMask = in; };
  void LoadMask();

  double ClusterBleeding(cluster cl, int clusterIndex, TF2 * , TH2 *);
  vector< pair<double, double> > GetCritialPoints(TF2 * , TH2 * , double);
  double ConcentricBleeding(cluster, TH2 * , double);
  double IsoVolumeAnalysisBleeding(cluster, TH2 * , double);
  double IsoEnergyContoursAnalysisBleeding(cluster cl, int clusterIndx);
  void GetInternalExternalVolume(TGraph * contour, TGraph2D * cluster, double & int_volume, double & ext_volume);
  bool PolygonIsClosed(TGraph *, double);
  void DumpPolygon(TGraph * poly);
  bool PolygonIsInside ( TGraph * polyA_is_inside, TGraph * polyB );
  TGraph2D * ExtractParOfTGraph2D ( TGraph2D *, TGraph * );
  double HeightMissingFactor();

  bool PushToNContourQueue ( int nConts );
  bool NContourQueueIsFull();

private:

  AllBlobsContainer * m_aB;
  double m_timepixClock;

  // Some control objects also rewinded on a frame per frame basis
  queue<int> m_nContoursQueue;
  int m_nContoursQueueSum;
  int m_nContoursQueueState;
  int m_nContoursCrownLevel;

  // Output
  int m_frameId;
  int m_nClusters; // cluster per frame
  vector<int> m_innerPixels;      // [pixels]
  vector<int> m_clusterSize;      // [pixels]
  vector<int> m_clusterTOT;       // [TOT]
  vector<double> m_clusterEnergy; // [MeV]
  vector<double> m_clusterEnergyReg; // [MeV]

  vector<int> m_clusterGeoX;
  vector<int> m_clusterGeoY;
  vector<int> m_pixelTOT;
  vector<double> m_circleArea;
  vector<double> m_ellipseArea;
  vector<double> m_circleOverEllipse;
  vector<double> m_hullCordsFraction;
  vector<double> m_ebleed;               // Bleeding energy per cluster
  vector<double> m_ebleed2;               // Bleeding energy per cluster
  vector<double> m_ehot;               // Bleeding energy per cluster
  vector<double> m_ehot2;               // Bleeding energy per cluster

  vector<double> m_volumeRemainsLevel0;
  vector<double> m_ebleedOverCluster;    // Bleeding energy per cluster
  vector<int> m_volcanoBit;
  vector<double> m_areaSkirtTotalRatioFactors;
  vector<double> m_clusterLength;

  // special selection
  typedef struct {
	  int minInner;
	  int maxInner;
	  int minClusterSize;
	  int maxClusterSize;
	  int minClusterTOT;
	  int maxClusterTOT;
	  double minCircleOverElipse;
	  double maxCircleOverElipse;
	  double minHullCord;
	  double maxHullCord;
	  int meshDiv;
	  int borderExclusion;
	  double multipeak_divider;
	  double bleeding_derivative_thl;
	  double bleeding_derivative_counts;
	  double bleed_min;
	  double bleed_max;
	  double bleedToClusterE_cut;
	  double bleedConcentricDiffLimit;
	  int ndivitionsTGraph2D;
	  double distanceMultiplierClosingPoly;
  } sel_type;

  // configuration values
  sel_type m_sel;
  bool m_acceptBorder;
  bool m_inputMask;
  int * m_Mask;
  int m_sizex;
  int m_sizey;
  int m_extraDrawing;

  // 3D clusters
  std::vector<TGraph2D *> m_graph2DVector;
  std::vector<TGraph *> m_graphVector;
  vector<vector<TGraph * > > m_allslices;
  vector<vector<TGraph * > > m_clusterCountours;
  vector<vector<TF1 * > > m_allkernels;

  std::vector<TF2 *> m_tf2Vector;
  std::vector<TH2 *> m_th2Vector;

  vector<TCanvas *> m_extraCanvas;

  // A resistivity map.  Average cluster size
  TH2F * m_clusterSizeMap;
  TH2I * m_clusterSizeEntries;
  TH2F * m_clusterTOTMap;
  TH2F * m_clusterEMap;


  ClassDef(HimacFeb2012, 1)
};

#endif
