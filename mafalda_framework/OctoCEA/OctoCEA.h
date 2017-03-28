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

#ifndef OctoCEA_h
#define OctoCEA_h

#include "MPXAlgo/MediPixAlgo.h"

// This algorithm is using an object put in the StoreGate
//  by BlobsFinder. I need to know the object.
#include "BlobsFinder/BlobsFinder.h"

//#define __HEAD_INDX 0
//#define __TAIL_INDX 1

#define __RECO_UNDEFINED    0
#define __RECO_MIP    5
#define __RECO_CLOUD 10

#define __chip0_rightedge 256
#define __chip1_rightedge 512
#define __chip2_rightedge 768
#define __chip3_rightedge 1024

class OctoTracksContainer;

class OctoCEA : public MediPixAlgo {

public:

	OctoCEA();
	virtual ~OctoCEA() { };

	// You ought to implement Init(), Execute() and Finalize()
	//  when you inherit from MediPixAlgo.  This model gives you
	//  direct access to data and services.
	void Init();
	void Execute();
	void Finalize();

	vector<double> CalculateResiduals(set<pair<int, int> > completeAssembly, double slope, double cut);


private:

	AllBlobsContainer * m_aB;

	// Config
	Int_t m_minNPixels;
	Int_t m_minNPixelsClouds;
	Int_t m_minClusterSize;
	Float_t m_joinCutPixelTolerance;
	Float_t m_joinSlopeAngleTolerance;

	double m_cutAngleValue;
	double m_occupancy;

	// Ntuple contents
	vector<Float_t> m_segmentsAngles;
	vector<Float_t> m_distanceToCenter;
	vector<Float_t> m_alphaAngle;
	vector<Float_t> m_chi2;
	vector<Float_t> m_chi2Prob;

	vector<double> m_perpDist;
	vector<double> m_perpDistSegments;
	vector<double> m_residual;
	vector<double> m_residualAligned;
	vector<double> m_odf;


	int nSegments;
	vector<int> clusterSize;


	vector<Int_t> m_clusterSize;
	Int_t m_typeOfRecoFrame;

	vector<int> m_assemblyTime;
	Int_t m_frameId;

	// Track description to send to StoreGate
	OctoTracksContainer * m_tracksPtr;

	ClassDef(OctoCEA, 1)
};



/*
 *  OctoClusterContainer is a class containing all the
 *  new OctoClusters
 */

enum octoclustertype {
	_NOTYPE = 0,
	_RECONSTRUCTED_MIP = 1,
};

class OctoTracksContainer : public CandidateContainer {

public:

	OctoTracksContainer(MediPixAlgo *);
	virtual ~OctoTracksContainer(){};

	// Store a new entry
	void PushOneTrackIndxs(vector<int>);
	// Store characteristics
	void PushOneTrackCuts(vector<double>);
	// Get the number of tracks stored
	//int GetNumberOfTracks(){return m_numberOfTracks;};
	int GetNumberOfTracks(){return (int)m_trackCandidatesIndexes.size();};
	// Get the vector of indexes which match the indexes in the ClusterContainer !
	vector<int> GetOneTrackIndxs(int indx){return m_trackCandidatesIndexes[indx];};
	vector<double> GetOneTrackCuts(int indx){return m_trackCandidatesCuts[indx];};

private:
	// Maps track indx to all the indexes of its constituents in the ClusterConatiner
	// Coming from BlobsFinder
	map<int, vector<int> > m_trackCandidatesIndexes;
	// Store some characteristics of the candidates
	map<int, vector<double> > m_trackCandidatesCuts;

	int m_numberOfTracks;

	ClassDef(OctoTracksContainer, 1)
};


#endif
