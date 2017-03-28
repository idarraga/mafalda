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

#ifndef FEI3Mips_h
#define FEI3Mips_h

#include "MPXAlgo/MediPixAlgo.h"

// This algorithm is using an object put in the StoreGate
//  by BlobsFinder. I need to know the object.
#include "BlobsFinder/BlobsFinder.h"

#include <TGraph2D.h>

#include <vector>

using namespace std;

// corner masks
#define __LEFT_CORNER_MASK_1  0x6  //  110
#define __LEFT_CORNER_MASK_2  0x2  //  010
#define __RIGHT_CORNER_MASK_1 0x5  //  101
#define __RIGHT_CORNER_MASK_2 0x5  //  101

#define __HEAD_INDX 0
#define __TAIL_INDX 1

class TGraph2D;

class FEI3Mips : public MediPixAlgo {

public:

	FEI3Mips();
	virtual ~FEI3Mips() { };

	// You ought to implement Init(), Execute() and Finalize()
	//  when you inherit from MediPixAlgo.  This model gives you
	//  direct access to data and services.
	void Init();
	void Execute();
	void Finalize();

	vector<pair<double, double> > FindLimitPixels(blob);
	vector<pair<double, double> > FindLimitPixels_MaxMinXY(blob);
	//-i
	//pair<Float_t, Float_t> GetParallelLineAt(Float_t, Float_t, pair<Int_t, Int_t>, bool *);
	pair<Float_t, Float_t> GetParallelLineAt(Float_t, Float_t, pair<Float_t, Float_t>, bool *);
	Int_t NumberOfBitsOnBeforeDoubleZero(UInt_t);
	Int_t GuessNumberOfDeltaRaysFromUnbalance(Int_t &, Int_t &);
	vector<Float_t> ExtractDeltaRayEnergy(blob);
	void ClearThisAlgo();

private:

	AllBlobsContainer * m_aB;
	Int_t m_minNPixels;
	Float_t m_minMipLength;
	Int_t m_nDivisions;
	Float_t m_pixSizeX;
	Float_t m_pixSizeY;
	Float_t m_dRayBalSearch;
	Int_t m_guardDistanceX;
	Int_t m_guardDistanceY;
	int m_maxYWidth;
	Float_t m_FEI3_PIXELSIZE_YX_FRACTION;
	int m_newMeshDiv;

	// For output
	vector<Float_t> m_chargeWeights;
	vector<Float_t> m_lvl1Weights;

	// charge per section
	vector< vector<Int_t> > m_chargePerSection;
	vector< Float_t > m_meanChargePerSection;
	vector< Float_t > m_sectionCenter_x;
	vector< Float_t > m_sectionCenter_y;

	// store the coordinates too
	vector< vector< pair<Int_t, Int_t> > > m_XYPerSection;

	vector<Float_t> m_deltaRayTOT;
	vector<Float_t> m_mipLength;
	double m_totalMIPsLength;

	Int_t m_nDeltaRay; // counter over the whole run
	Int_t m_nDeltaRaySoft; // counter over the whole run
	Int_t m_nDeltaRayPerFrame; // to store frame per frame
	Int_t m_nDeltaRaySoftPerFrame; // to store frame per frame

	Int_t m_distanceCutHot;
	Int_t m_distancePerpHot;

	Int_t m_clusterSize;
	Double_t m_alphaAngle;
	Int_t m_frameId;

	// graphs
	// 3D clusters
	vector<TGraph2D *> m_graph2DVector;
	TCanvas * m_extraCanvas;

	ClassDef(FEI3Mips, 1)
};

#endif
