 /*
 * 	Copyright 2009 John Idarraga
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

#ifndef BlobsFinder_h
#define BlobsFinder_h

#include <vector>
#include <map>
#include <set>
#include <list>
#include <utility>

#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>

#include "MPXAlgo/MediPixAlgo.h"

using namespace std;

/* ******************************* */
// valgrind
//#include <valgrind/valgrind.h>
/* ******************************* */

using namespace std;
#define __REJECTED_BY_LEVEL_1  101
#define __BLOB_GOOD_TO_STORE   201

#define MAX_BLOBS_PER_FRAME 500000
// numbering of quadrants
#define Q_FIRST  1
#define Q_SECOND 2
#define Q_THIRD  3
#define Q_FOURTH 4
// same thing but index instead
#define Q_FIRST_INDX  0
#define Q_SECOND_INDX 1
#define Q_THIRD_INDX  2
#define Q_FOURTH_INDX 3

///////////////////////////////////
// A few definitions and vars for 
//  spiral pattern flow
#define X_STEP 10
#define Y_STEP 20
#define NEIGHBORS_DISC_0_PIXELS   8
#define NEIGHBORS_DISC_1_PIXELS  24
#define MAX_DISCONTINUITY       200
#define MIN_INNERPIXEL            4

#define TOT_OVERFLOW_COUNTS    11810

// End of blob conditions
enum {
	__blob_finished = 0,
	__blob_too_big_giveupframe
};

typedef struct {
	Int_t xySwitch;
	Int_t localMax;
	Int_t xCntr;
	Int_t yCntr;
	Int_t dir;
} steerSpiral ;
///////////////////////////////////

///////////////////////////////////
// Properties of a blob
typedef struct {

	Int_t nPixels;
	Int_t nInnerPixels;
	Int_t clusterSize;
	Int_t clusterTOT;
	Int_t totalCharge; // Alias for "clusterTOT", this is a bad choice of name.
	Double_t clusterEnergy; // in keV if calibration present
	Int_t width_x;
	Int_t width_y;
	Float_t geoCenter_x;
	Float_t geoCenter_y;
	Float_t weightedCenter_x;
	Float_t weightedCenter_y;
	Int_t boxArea;

	Float_t circleArea;
	Float_t ellipseArea;
	Float_t ellipseA;
	Float_t ellipseB;
	Float_t rotAngle;
	Float_t chisquare_OverDof;
	Float_t chisquare;
	Float_t NDF;
	Float_t fitSlope;
	Float_t fitCut;
	Float_t balanceToMin;
	Float_t balanceToMax;
	Float_t minToGeoCenter;
	Float_t maxToGeoCenter;

	vector<int> lvl1;

} blobProperties ;
///////////////////////////////////


/** 
 * Tagging with the following clasification will be done by
 * PRBasicSpecies.  The blobs already contain the variable where this
 * will be inserted
 */
enum blobtype {
	_NOTYPE_ASSIGNED = 0,
	_SINGLE_HIT, // not marked
	_DOUBLE_HIT, // not marked
	_TRIPLE_HIT, // not marked // 3
	_QUAD_HIT,   // not marked // 4
	_LONG_GAMMA, // not marked // 5
	_MIP,        // green arrow
	_HEAVY_BLOB, // black circle // 7
	_HEAVY_TRACK,// green circle
	_CURLY,      // black arrow
	_NOT_A_BASIC_TYPE, // magenta arrow | EColor::kMagenta // 10 - not a basic type from here
	_SOFT_VERTEX,// blue arrow | EColor::kBlue
	_HARD_VERTEX,// red arrow | EColor::kRed
	_LONG_CURLY, // violet arrow | EColor::kViolet
	_OVERLAP,    // orange arrow | EColor::kOrange
	_UNKNOWN,
	_NTYPES
};

/**
 * One blob/cluster contents
 */

class blob {

public:

	blob(){};
	blob(const blob &);
	blob(const blob &, set<pair<int, int > >);
	blob(vector<pair<int, int > >, map<int, int>, int);
	TString GetTypeAsString();
	//void SetBlobType(blobtype);
	void SetBlobType(blobtype);
	blobProperties GetBlobProperties(){return bP;};
	blobtype GetBlobType(){return btype;};

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Access to cluster constituents methods
	set< pair<Int_t, Int_t> > GetContentSet(){return blobContentSet;};                 // OK --> Can be used to know the (x,y) constituents of a cluster
	list< pair < pair<Int_t, Int_t>, Int_t > > GetBlobContent(){return blobContent;};  // NO --> Not to be used by the user.
	list< pair < pair<Int_t, Int_t>, Int_t > > GetClusterDescription(){return clusterDescription;};                     // OK --> Full information [ (x,y), Counts ]
	list< pair < pair<Int_t, Int_t>, Int_t > > GetClusterDescriptionToA(){return clusterDescriptionToA;};                     // OK --> Full information [ (x,y), Counts ]
	list< pair < pair<Int_t, Int_t>, Int_t > > GetClusterDescriptionFastToA(){return clusterDescriptionFastToA;};                     // OK --> Full information [ (x,y), Counts ]

	list< pair < pair<Int_t, Int_t>, Int_t > > GetClusterDescriptionCalibrated(){return clusterDescriptionCalibrated;}; // OK --> Full calibrated information (if calibration loaded) [ (x,y), Calib Energy ]
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	blob & DeletePixelEntry(pair<int, int>, int padSizex, int padSizey, double chisquare_OverDof_rotation = 2.0);
	void SetContentSet(set< pair<Int_t, Int_t> > cs){blobContentSet = cs;};
	void SetBlobContent(list< pair< pair<Int_t, Int_t>, Int_t > > c){blobContent = c;};
	void SetClusterDescription(list< pair < pair<Int_t, Int_t>, Int_t > > cd){clusterDescription = cd;};

	void SetClusterEnergy(Double_t e){ bP.clusterEnergy = e; };
	void SetPixelEnergy(Double_t e, pair<int,int> pix) {
		pair < pair<int, int>, int > p = make_pair( pix, (int) e );
		clusterDescriptionCalibrated.push_back( p );
	};

	int CalculateProperties(int xdim, int ydim, double chisquare_OverDof_rotation);
	void CalcMinMaxGeoCenter(double, double);
	bool SingleHit();
	bool DoubleHit();
	bool TripleHit();
	bool QuadHit();
	void ClearOneBlobData();
	void AllDump();

	// TODO: change all these parameters by a structure
	double GetBalance(Int_t, Int_t, Int_t, Int_t, double *, double *, Int_t *, Int_t *, double *, double *);
	double CalcMeanForBlob(Int_t, Int_t, Int_t, Int_t, Int_t * , Int_t *, Int_t * , Int_t *, double*, double *);
	//TString Dump();

	// this will go private ... TODO !!!
	// properties
	// <x, y> , trackCounter --> neighbours of '*'   8  7  6
	//                                               1  *  5
	//                                               2  3  4
	// A <list> keeps entries in order
	// WARNING The second member of the pair is not the TOT(or counts).  Is
	//  the ordering in the spiral.  Do not use this if you need the TOT or Counts !
	list< pair < pair<Int_t, Int_t>, Int_t > > blobContent;
	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator blobItr;
	// This <set> is used to verify if and entry already exists in the list
	//  stdlib does better than me ;)
	set< pair<Int_t, Int_t> > blobContentSet;
	// This list really contains the X,Y,C information
	list< pair < pair<Int_t, Int_t>, Int_t > > clusterDescription;
	// This list really contains the X,Y,C information
	list< pair < pair<Int_t, Int_t>, Int_t > > clusterDescriptionToA;
	// This list really contains the X,Y,C information
	list< pair < pair<Int_t, Int_t>, Int_t > > clusterDescriptionFastToA;
	// This list really contains the X,Y,Energy information
	list< pair < pair<Int_t, Int_t>, Int_t > > clusterDescriptionCalibrated;

	blobProperties bP;
	//blobtype btype;
	blobtype btype;

};
// A "blob" can be called "cluster" too
typedef blob cluster;


/**
 * Container for all blobs
 * Going to the store gate
 *
 */

/* *************************************** */
class AllBlobsContainer : public CandidateContainer {

public:

	AllBlobsContainer(MediPixAlgo *);
	virtual ~AllBlobsContainer(){};
	vector<blob> GetBlobsVector(){return allBlobs;};
	void SetBlobType(Int_t, blobtype);
	blobtype GetBlobType(Int_t);
	TString GetTypeAsString(Int_t);
	void PushBack(blob b){ allBlobs.push_back(b); };
	void SetClusterEnergy(Double_t e, int index){
		allBlobs[index].SetClusterEnergy( e );
	};
	void SetPixelEnergy(Double_t e, int index, pair<int,int> pixel){
		allBlobs[index].SetPixelEnergy(e, pixel);
	};

public:
	// should get to be private but I need to change things in the code
	vector<blob> allBlobs;

	ClassDef(AllBlobsContainer, 1)
};
/* *************************************** */

class BlobsFinder : public MediPixAlgo {

public:

	BlobsFinder();
	virtual ~BlobsFinder() { };

	void Init();
	void Execute();
	void Finalize();

	void SetBorderExclusion(Int_t);
	void SetDiscontinuityTolerance(Int_t);
	void SetMaxOcc(Double_t max) { m_maxOcc = max; };
	Int_t GetLengthOfSpiral(Int_t);

	int StartBlob(Int_t, Int_t);
	int FollowBlob();
	pair<Int_t, Int_t> GetNextPosition(pair<Int_t, Int_t>);
	bool IsUnsafePosition(pair<Int_t, Int_t>);
	void ClearOneBlobData();
	//int CalculateProperties();
	pair<Int_t, Int_t> RotateVector2D(Int_t, Int_t, double);
	void DumpInfo();
	void DumpProperties();
	void RewindSpiral();

	void FillClusterDescription();

private:

	// Mask handler

	void CopyPropertiesToNtupleData();
	//Bool_t LinearRegression(Float_t *, Float_t *, Float_t *, Float_t, Float_t);


	// one blob
	blob oneBlob;
	// all blobs
	AllBlobsContainer * aB;
	// spiral manager
	steerSpiral m_spiral;
	// spiral ends at
	Int_t m_spiralEnds;
	Int_t m_discontinuityTolerance;
	Int_t m_lvl1Cut;
	double m_chisquare_OverDof_rotation;

	// set limits for search in pad
	typedef struct {
		Int_t colmin;
		Int_t colmax;
		Int_t rowmin;
		Int_t rowmax;
	} limits;
	limits m_padLimits;
	Int_t m_excludeBorder;

	// 	Occupancy limit
	double m_maxOcc;
	double m_maxClusterSize;

	// README !!!
	// Vars to save in ntuple:
	// In the ntuple I need to save a value for each blob.  I am
	//  guessing that I won't have more than MAX_BLOBS_PER_FRAME blobs
	//  per frame.  If that happens the algo aborts and you will have
	//  to change that setting and recompile.
	// When arrays go to an ntuple you need to have fixed size.
	//  that's why this copy of the same structure above is
	//  here.  It won't get replicated in this algorithm.  The waste of
	//  space remains static.  Harmless.
	// FIXME
	// use std::vectors here

	Int_t m_bP_nPixels[MAX_BLOBS_PER_FRAME];
	Int_t m_bP_clusterSize[MAX_BLOBS_PER_FRAME];
	Int_t m_bP_totalCharge[MAX_BLOBS_PER_FRAME];
	Int_t m_bP_width_x[MAX_BLOBS_PER_FRAME];
	Int_t m_bP_width_y[MAX_BLOBS_PER_FRAME];
	Float_t m_bP_geoCenter_x[MAX_BLOBS_PER_FRAME];
	Float_t m_bP_geoCenter_y[MAX_BLOBS_PER_FRAME];
	Float_t m_bP_weightedCenter_x[MAX_BLOBS_PER_FRAME];
	Float_t m_bP_weightedCenter_y[MAX_BLOBS_PER_FRAME];
	Int_t m_bP_boxArea[MAX_BLOBS_PER_FRAME];
	Float_t m_bP_circleArea[MAX_BLOBS_PER_FRAME];
	Float_t m_bP_rotAngle[MAX_BLOBS_PER_FRAME];
	Float_t m_bP_ellipseArea[MAX_BLOBS_PER_FRAME];
	Float_t m_bP_balanceToMin[MAX_BLOBS_PER_FRAME];
	Float_t m_bP_balanceToMax[MAX_BLOBS_PER_FRAME];
	Int_t m_bP_nBlobs;
	Int_t m_bP_rejectedByLvl1;

	ClassDef(BlobsFinder, 1)
};

// "BlobsFinder" is called "ClusterFinder" too
typedef BlobsFinder ClusterFinder;

enum { 
	__CENTER = 0,
	__W,
	__SW,
	__S,
	__SE,
	__E,
	__NE,
	__N,
	__NW, // --> 8
};

#endif
