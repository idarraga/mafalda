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

#ifndef __Timepix3VolcanoEdge_cpp
#define __Timepix3VolcanoEdge_cpp

#include "Timepix3VolcanoEdge.h"
#include "MAFTools.h"

using namespace MSG;

ClassImp(Timepix3VolcanoEdge)

Timepix3VolcanoEdge::Timepix3VolcanoEdge() : MediPixAlgo(), CalibrationLoader(this) {

	// This value will be overridden by the configuration because it'll registered
	//  as a ConfigurationValue in the Init member of this class.
	m_minNPixels = 5;
	// selection
	m_minInner = 2;

	// Map TOT --> energy
	m_TOTEnergyMap.clear();
	// prepare for 10 bits --> max TOT 1024
	for(int i = 1 ; i <= __max_tot_range ; i++) m_TOTEnergyMap[i] = vector<double>(); // initialize with m_selSizeX zeroes

}

void Timepix3VolcanoEdge::Init(){

	Log << MSG::INFO << "Init function !" << endreq;


	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("clusterTOT", &m_clusterTOT);
	getMyTree()->Branch("clusterEnergy", &m_clusterEnergy);

	// A configuration value that can be tuned from the Viewer
	RegisterConfigurationValue(&m_minNPixels, "minNPixels");

}

void Timepix3VolcanoEdge::Execute(){

	// Ask the store gate if the previous algorithm (BlobsFinder --> reponsible for clustering)
	//  sent any objects to the StoreGate.
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");
	if(lastObject == 0)
		return;

	// If so, get the pointer to the last object.  BlobsFinder, as it comes out of the box, sends
	//  a single object containing all clusters.
	m_aB = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject-1);

	// AllBlobsContainer is a box full of Blobs(Clusters). Now we can iterate over all clusters inside.
	vector<blob> blobsVector = m_aB->GetBlobsVector();
	Log << MSG::INFO << "Number of blobs from clustering = " << (Int_t) blobsVector.size() << endreq;
	vector<blob>::iterator blobsItr = blobsVector.begin(); //allBlobs.begin();

	cluster cl;
	for( ; blobsItr != blobsVector.end() ; blobsItr++) {

		cl = *blobsItr;

		// Limit all this to clusters with a minimum size.
		// Note that m_minNPixels can be configured through the Viewer
		//  so you can reprocess and check results online.
		if(cl.bP.nPixels < m_minNPixels)
			continue;

		if(cl.bP.nInnerPixels < m_minInner)
			continue;

		// Finish selection.  Signal the cluster
		Highlighter * highLight = new Highlighter(cl.bP.geoCenter_x,
				cl.bP.geoCenter_y,
				"arrow", this);
		highLight->SetLineColor(kBlack);
		highLight->SetLineWidth(1);
		FillValuesForDisplay(highLight, cl);
		PullToStoreGateAccess(highLight, MPXDefs::DO_NOT_SERIALIZE_ME);

		// If the cluster passes the minimum requirements loop over the
		// the constituents of the cluster --> ClusterDescription
		list< pair < pair<int, int>, int > > cl_des = cl.GetClusterDescription();
		list< pair < pair<int, int>, int > >::iterator i = cl_des.begin();

		// Store the cluster TOT for output
		m_clusterTOT.push_back( cl.bP.clusterTOT );

		double calib_edep = 0.0, clusterEdep = 0.;
		int tot = 0;
		pair<int, int> pix;

		for ( ; i != cl_des.end() ; i++) {

			// pixel coordinates and tot
			pix = (*i).first;
			tot = (*i).second;

			// Use calibration to obtain E = Surrogate(TOT) for this pixel
			calib_edep = CalculateAndGetCalibEnergy(pix, tot);

			// Fill the TOT Vs. Energy map
			m_TOTEnergyMap[tot].push_back( calib_edep );

			// Calculate the energy of the cluster
			clusterEdep += calib_edep;

		}

		// Store the cluster Energy calculated in the previous loop
		m_clusterEnergy.push_back( clusterEdep );

	}

	// Fill the output tree of this algorithm
	getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_clusterEnergy.clear();
	m_clusterTOT.clear();

}

void Timepix3VolcanoEdge::Finalize() {

	Log << MSG::INFO << "Finalize function !" << endreq;

	// X axis
	vector<double> xg(__max_tot_range, 0);
	vector<double> xg_err(__max_tot_range, 0);

	// Calculate the averages for each bin
	vector<double> Y(__max_tot_range, 0);
	vector<double> Y_err(__max_tot_range, 0);

	for(int i = 1 ; i <= __max_tot_range ; i++) {

		double mean;
		double stdv_sel;

		Log << MSG::INFO << "[" << i << "] --> " << m_TOTEnergyMap[i].size() << " entries" << endreq;
		//MAFTools::DumpVector( m_TOTEnergyMap[i] );

		if( !m_TOTEnergyMap[i].empty() ) {
			mean = MAFTools::CalcMean( m_TOTEnergyMap[i] );
			// As a measurement of the dispersion here I will consider only the point within 1 stdev from the mean
			stdv_sel = MAFTools::CalcStdDev( m_TOTEnergyMap[i] );
		} else {
			mean = 0.;
			stdv_sel = 0.;
		}

		// X
		xg[i-1] = i;

		// Y
		Y[i-1] = mean;
		Y_err[i-1] = stdv_sel;
		Log << MSG::INFO << "mean = " << mean << ", stdv_sel = " << stdv_sel << endreq;
	}

	TGraphErrors * gy = new TGraphErrors(__max_tot_range, &xg[0], &Y[0], &xg_err[0], &Y_err[0] );
	getMyROOTFile()->cd();
	gy->Write("TOTEnergy");

}

#endif
