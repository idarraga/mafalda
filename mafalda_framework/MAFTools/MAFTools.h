/*
 * 	Copyright 2012 John Idarraga
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

#ifndef MAFTools_h
#define MAFTools_h

#include <set>
#include <map>
#include <iostream>
#include <vector>
#include <queue>

#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "BlobsFinder/BlobsFinder.h"

//class blob;
class MediPixAlgo;
class Highlighter;
class TGraph2D;
//class AllBlobsContainer;
//typedef struct steerSpiral;

using namespace std;

enum cluster_locations { __HEAD_INDX = 0, __TAIL_INDX = 1 };
#define __fraction_of_height_range_id 0.3

namespace MAFTools {

// overloaded when one can pass the list of points in a vector
Bool_t LinearRegression(double * slope, double * cut,
		double * chisquare_OverDof, double * chisquare, double * ndf,
		double initslope, double initycross,
		std::vector< std::pair<Int_t, Int_t> > blobContent);

// overloaded where one can pass the list of points as a set
Bool_t LinearRegression(double * slope, double * cut,
		double * chisquare_OverDof, double * chisquare, double * ndf,
		double initslope, double initycross,
		std::set< std::pair<Int_t, Int_t> > blobContent);

// avoid one unique point (guaranteed by the use of a set) and perform the Linear Reggresion
Bool_t LinearRegression(double * slope, double * cut,
		double * chisquare_OverDof, double * chisquare, double * ndf,
		double initslope, double initycross,
		std::set< std::pair<Int_t, Int_t> > blobContent,
		pair<Int_t, Int_t> pairAvoid);

// Calculate the are of a triangle given 3 vertexes
double TriangleArea(double * x, double * y);
// Convert a TGraph2D into TGraph crashing the Z component with zero suppression
TGraph * GetTGraphZeroSuppression(TGraph2D * g2D);

Int_t XYtoC(pair<Int_t, Int_t>, Int_t);
Int_t XYtoC(Int_t, Int_t, Int_t);
Int_t XYtoX(Int_t, Int_t, Int_t);
Int_t XYtoX(pair<int, int>, Int_t);
Int_t ConvertXYPositionToInteger(std::pair<Int_t, Int_t> position, Int_t dimX, Int_t dimY);
Int_t ConvertXYPositionToInteger(Int_t positionX, Int_t positionY, Int_t dimX, Int_t dimY);
pair<int, int> XtoXY(int X, int dimX);

Double_t CalcDistance(std::pair<Double_t, Double_t>, std::pair<Double_t, Double_t>);
Double_t CalcDistance(std::pair<int, int>, std::pair<int, int>);
Double_t CalcDistance(Double_t, Double_t, Double_t, Double_t);
Double_t CalcPerpDistancePointToLineWithSign(Float_t slope, Float_t cut, pair<Int_t, Int_t> point);
Double_t CalcPerpDistancePointToLineWithSign(Float_t slope, Float_t cut, Int_t px, Int_t py);

Double_t CalcPerpDistanceToLine(double, double, pair<Int_t, Int_t>);
Double_t CalcPerpDistanceToLine(double, double, Int_t, Int_t);

//template<typename T> double CalcStdDev(vector<T> v);
double CalcStdDev(vector<double> v);
double CalcStdDev(vector<int> v);
double CalcMean(vector<double> v);
double CalcMean(vector<int> v);
void DumpVector(vector<double> v);
void DumpVector(vector<int> v);

/* Drawing */
void FillValuesForDisplay(Highlighter * hl, blob bl);
//void FillValuesForDisplay(Highlighter * hl, blob bl, int nextas = 0, double * extras = 0x0);
void FillValuesForDisplay(Highlighter * hl, blob bl, int nextras, vector<TString> names, vector<double> extras);
void DrawLine(MediPixAlgo *, double, double, double, double, Int_t, Int_t, EColor);
void DrawHull(MediPixAlgo *, vector< pair<int, int> > h, int lcolor, float lshift);
void DrawLocalCoordSystem(MediPixAlgo *, blob);

/* Word handling */
std::string DumpWordInBinary(UInt_t word);

// Blobs handling
vector<pair<double, double> > FindLimitPixels(blob);

// Pixel inner
bool PixelIsInner(pair<int,int>, blob);
vector<pair<int,int> > GetClusterSkeletton_InnerNegative(blob);
vector<pair<int,int> > GetClusterSkeletton_InnerNegativeAndAdjacent(blob);
vector<pair<int,int> > GetVectorOfInnerPixels(blob);
set<pair<int,int> > GetSetOfInnerPixels(blob);
bool PixelInListIsAClusterConstituent(blob, int *, int, int);
bool PixelInListIsAClusterConstituent(blob, vector<int> , int, int);

/*
 * Check if a blob appears within a certain distance from the edges
 * Parameters:
 *  - blob
 *  - vector containing farthest pixels in the blob (see FindLimitPixels(blob))
 *  - Width of sensor in units of pixels
 *  - Height of sensor in units of pixels
 */
Bool_t BlobAtADistanceFromEdges(blob, vector<pair<double, double> >, Int_t, Int_t, Int_t, Int_t);

/**
 * Clustering stuff
 */
int DigitalBlobReclustering(int, AllBlobsContainer *, blob, int, int, std::map<int,int>);
bool StartBlob(int, int, int, AllBlobsContainer *, blob *, int, int, std::map<int,int>);
bool FollowBlob(blob *, steerSpiral *, int, int, int, std::map<int,int>);
void RewindSpiral(steerSpiral *);
Int_t GetLengthOfSpiral(Int_t);
pair<Int_t, Int_t> GetNextPosition(pair<Int_t, Int_t>, steerSpiral *);
bool IsUnsafePosition(pair<Int_t, Int_t>, int, int);
void ClearOneBlobData(blob *);
set< pair<int, int> > VectorToSet(vector<pair<int, int> >);
pair<int, int> GetNextPosition(pair<int, int> current, steerSpiral & spiral);
void RewindSpiral(steerSpiral & spiral);

// Dump contents of a blob
void DumpBlobContents(blob);

int FindDiscontinuity(blob);

// More functions on Clusters
TH2F * CropCluster (blob b, long frameId, TString type = "TOT", TString extraname = "");
TGraph * ClusterProfile(blob b, long frameId, TString dir = "X", TString type = "TOT");
vector<TGraph *> GetAllProfiles(blob b, long frameId, TString dir, TString type);
TGraph2D * ConvertClusterToGraph2D(blob, int, TString, TString ctype = "TOT", TString appendname = ""); // an identifier, can be string or int
TGraph2D * ConvertClusterToGraph2D(blob, int, long, TString ctype = "TOT", TString appendname = "");
list< pair < pair<int, int>, int > >::iterator GetLowestToA(list< pair < pair<int, int>, int > > );

TH2I * ConvertClusterToTH2I(blob, long, TString appendname = "");
pair<int, int> GetBottomLeftExtreme(cluster);
pair<int, int> FindBottomLeftCorner(cluster);
int CheckDiscontinuityWithModifiedMatrix(blob bl, map<int,int> totmap, int discontinuityTolerance, int, int, MediPixAlgo * algo);
int FindSubclusters (blob b, map<int, int> totmap, double multipeak_divider, int width, int height, MediPixAlgo * algo);

vector<pair<double, double> > FindLimitPixels_MaxMinXY(blob);

blob BuildBlob(TGraph2D * g);
TH1 * ExtractSlice(blob b, int, TString dir = "X");
int IdendifyPeaks_A1(TH1 *, vector<double> *);
void SimpleLinearReggresion(queue<int>, double *);
int GetSign(double);
int SmallMin(int, int);
bool ClusterTouchesBorder(cluster, int, int, int);
double GetTOTOfHotestsPixel(cluster cl);

bool GoodCandidateGraph2D(blob);
void PrintMatrix(int **, int);
int ** RefineGridOnce(int ** c, int * , int *);

list< pair < pair<int, int>, int > >::iterator FindHotestPixel(list< pair < pair<int, int>, int > > *);
list< pair < pair<int, int>, int > >::iterator FindColdestPixel(list< pair < pair<int, int>, int > > *, int hot = 12000);

double GetTOTOfHotestPixel(cluster cl);
double GetEnergyOfHotestPixel(cluster cl);
set<pair<int, int> > FindPixelsAboveAverageTOT(cluster);
set<pair<int, int> > FindPixelsBelowOrEqualAverageTOT(cluster);
set<pair<int, int> > FindPixelsBelow(cluster, int);
int MaskInMap(map<int, int> &, vector<pair<int, int> >, int);
int MaskInMap(map<int, int> &, set<pair<int, int> >, int);
int MaskInCluster(cluster & cl, vector<pair<int, int> > mask, int dimX);

typedef enum {
	__righthanded = 0,
	__reverse
} data_direction;

TH1F * GetHistogramFromTGraph(TGraph *, data_direction);
void FindMinMaxForFit(TH1 *, double, double &, double &);

// Math services
double GausFuncAdd(double * x, double * par);
TF1 * CreateKernelDensityFunction(int id, vector<double> hist, double bandwidth);
double DerivativeFivePointsStencil(TF1 * f, double x, double h);
//double DerivativeFivePointsStencil(queue<double>);
double DerivativeFivePointsStencil(queue<double>, double h = 1.);
TGraph * GetDerivative(TGraph *);
int GetCriticalPoints(TF1 * f, int nbins, vector<double> & min, vector<double> & max);
TGraphErrors * TriangularSmooth(TGraphErrors *);


}

#endif
