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

#ifndef __Timepix3Preliminary_cpp
#define __Timepix3Preliminary_cpp

#include "Timepix3Preliminary.h"
#include "MAFTools.h"
#include "MPXAlgo/Signals.h"

using namespace MSG;

ClassImp(Timepix3Preliminary)

Timepix3Preliminary::Timepix3Preliminary() : MediPixAlgo(), CalibrationLoader(this) {

	// This value will be overridden by the configuration because it'll registered
	//  as a ConfigurationValue in the Init member of this class.
	m_minNPixels = 5;
	m_aB = 0x0;
	m_highLight = 0x0;

	// selection
	m_minInner = 10;
	m_selSizeX = 20;
	m_selSizeY = 20;

	m_extraCanvas.clear();
	m_minDriftLimit = 0;
	m_maxDriftLimit = 1000;

	m_clusterSizeX.clear();
	m_clusterSizeY.clear();

	m_energyAverageXCntr = 0;
	for(int i = 0 ; i < m_selSizeX ; i++) m_energyAverageX[i] = vector<double>(); // initialize with m_selSizeX zeroes
	m_energyAverageYCntr = 0;
	for(int i = 0 ; i < m_selSizeY ; i++) m_energyAverageY[i] = vector<double>(); // initialize with m_selSizeY zeroes

}

void Timepix3Preliminary::Init(){

	Log << MSG::INFO << "Init function !" << endreq;


	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("clusterTOT", &m_clusterTOT);
	getMyTree()->Branch("clusterEnergy", &m_clusterEnergy);
	getMyTree()->Branch("driftingTime1", &m_driftingTime1);
	getMyTree()->Branch("clusterSizeX", &m_clusterSizeX);
	getMyTree()->Branch("clusterSizeY", &m_clusterSizeY);

	// A configuration value that can be tuned from the Viewer
	RegisterConfigurationValue(&m_minNPixels, "minNPixels");
	RegisterConfigurationValue(&m_minInner, "minInner");
	RegisterConfigurationValue(&m_minDriftLimit, "minDriftLimit");
	RegisterConfigurationValue(&m_maxDriftLimit, "maxDriftLimit");
	RegisterConfigurationValue(&m_selSizeX, "selSizeX");
	RegisterConfigurationValue(&m_selSizeY, "selSizeY");


}

void Timepix3Preliminary::Execute(){

	// Ask the store gate if the previous algorithm (BlobsFinder --> reponsible for clustering)
	//  sent any objects to the StoreGate.
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");
	if(lastObject == 0)
		return;

	// Clean up stuff
	vector<TCanvas *>::iterator ican = m_extraCanvas.begin();
	for( ; ican != m_extraCanvas.end() ; ican++ ) {
		delete (*ican);
	}
	m_extraCanvas.clear();
	m_graph2DVector.clear();

	// If so, get the pointer to the last object.  BlobsFinder, as it comes out of the box, sends
	//  a single object containing all clusters.
	m_aB = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject-1);

	// AllBlobsContainer is a box full of Blobs(Clusters). Now we can iterate over all clusters inside.
	vector<blob> blobsVector = m_aB->GetBlobsVector();
	Log << MSG::INFO << "Number of blobs from clustering = " << (Int_t) blobsVector.size() << endreq;
	vector<blob>::iterator blobsItr = blobsVector.begin(); //allBlobs.begin();

	cluster cl;
	int cntr = 0;
	bool extraDraw = false;
	bool skipFlag = false;

	for( ; blobsItr != blobsVector.end() ; blobsItr++) {

		cl = *blobsItr;

		// Limit all this to clusters with a minimum size.
		// Note that m_minNPixels can be configured through the Viewer
		//  so you can reprocess and check results online.
		if(cl.bP.nPixels < m_minNPixels)
			continue;

		if(cl.bP.nInnerPixels < m_minInner)
			continue;

		// Draw the cluster
		if(extraDraw) {
			m_graph2DVector.push_back(
					MAFTools::ConvertClusterToGraph2D(cl, 1,
							GetFrameId() + cntr, "ToA", "test")
			);
			extraDraw = false;
		}

		// If selected same some info about the cluster
		m_clusterSizeX.push_back( cl.bP.width_x );
		m_clusterSizeY.push_back( cl.bP.width_y );


		// If the cluster passes the minimum requirements loop over the
		// the constituents of the cluster --> ClusterDescription
		list< pair < pair<int, int>, int > > cl_des = cl.GetClusterDescription();
		list< pair < pair<int, int>, int > >::iterator i = cl_des.begin();

		// Store the cluster TOT for output
		m_clusterTOT.push_back( cl.bP.clusterTOT );

		double calib_edep = 0.0, clusterEdep = 0.;
		int tot = 0;
		pair<int, int> pix;

		unsigned int minToA = 0xFFFFFFFF;
		unsigned int maxToA = 0;
		unsigned int ToA = 0;
		list< pair < pair<int, int>, int > >::iterator minToAItr;
		list< pair < pair<int, int>, int > >::iterator maxToAItr;

		for ( ; i != cl_des.end() ; i++) {

			// pixel coordinates and tot
			pix = (*i).first;
			tot = (*i).second;

			ToA = GetToA(pix);
			//FastToA = GetFastToA(pix);

			if (minToA > ToA ) {
				minToA = ToA;
				minToAItr = i;
			}
			if (maxToA < ToA ) {
				maxToA = ToA;
				maxToAItr = i;
			}

			// Use calibration to obtain E = Surrogate(TOT) for this pixel
			calib_edep = CalculateAndGetCalibEnergy(pix, tot);

			// Calculate the energy of the cluster
			clusterEdep += calib_edep;

		}

		Log << MSG::ALWAYS << endreq <<  "At (" << cl.bP.geoCenter_x << "," << cl.bP.geoCenter_y << ")" << endreq;
		Log << MSG::ALWAYS << "ToA(max, min) : diff | (" << maxToA << ", " << minToA << ") : " << maxToA - minToA << endreq;

		// Too big differences are things that don't correspond to the same cluster.
		if ( maxToA - minToA > 100 ) continue;

		// Calculate now the drifting time
		double driftingTime = 0.; // ( maxToA - minToA ) * 25.0;

		// Deal with the case where all ToA are the same
		//  and we are bound to rely only on the FastToA
		if ( minToA == maxToA ) {

			// Search again now only for ToAFast
			unsigned int FastToA_min = 0xFFFFFFFF;
			unsigned int FastToA_max = 0x0;
			unsigned int FastToA = 0;
			Log << MSG::ALWAYS << "Zero difference !!!" << endreq;
			for ( i = cl_des.begin() ; i != cl_des.end() ; i++) {
				pix = (*i).first;
				FastToA = GetFastToA( pix );
				if( FastToA_min > FastToA ) FastToA = FastToA_min;
				if( FastToA_max < FastToA ) FastToA = FastToA_max;
			}

			driftingTime = (FastToA_max - FastToA_min) * (25.0/16.0);

		} else {

			Log << MSG::ALWAYS << "FastToA(max, min) : (" << GetFastToA( (*maxToAItr).first ) << ", " << GetFastToA( (*minToAItr).first ) << ")" << endreq;

			driftingTime  = ( maxToA*25.0 ) - ( GetFastToA( (*maxToAItr).first )*(25.0/16.0) );
			driftingTime -= ( minToA*25.0 ) - ( GetFastToA( (*minToAItr).first )*(25.0/16.0) );

			//driftingTime += GetFastToA( (*maxToAItr).first ) * (25.0/16.0);
			//driftingTime -= GetFastToA( (*minToAItr).first ) * (25.0/16.0);

		}

		Log << MSG::ALWAYS << "Difference = " <<  driftingTime << endreq << endreq;

		// Store the cluster Energy calculated in the previous loop
		m_clusterEnergy.push_back( clusterEdep );
		m_driftingTime1.push_back( driftingTime );

		if ( driftingTime < m_maxDriftLimit && driftingTime > m_minDriftLimit )  {
			// Finish selection.  Signal the cluster
			m_highLight = new Highlighter(cl.bP.geoCenter_x,
					cl.bP.geoCenter_y,
					"arrow", this);
			m_highLight->SetLineColor(kBlack);
			m_highLight->SetLineWidth(1);

			FillValuesForDisplay(m_highLight, cl);
			PullToStoreGateAccess(m_highLight, MPXDefs::DO_NOT_SERIALIZE_ME);


		}

		// Make some selection based on the geometry of the cluster
		if ( cl.bP.width_x == m_selSizeX && cl.bP.width_y == m_selSizeY
				//&&
				//cl.bP.width_x > 10 && cl.bP.width_y > 10
		) {

			Log << MSG::DEBUG << "Selected cluster" << endreq;
			// Finish selection.  Signal the cluster
			Highlighter * highSel = new Highlighter(cl.bP.geoCenter_x,
					cl.bP.geoCenter_y,
					"arrow", this);
			highSel->SetLineColor(kRed);
			highSel->SetLineWidth(1);

			FillValuesForDisplay(highSel, cl);
			PullToStoreGateAccess(highSel, MPXDefs::DO_NOT_SERIALIZE_ME);


			// Select the column and the row to consider
			int column = cl.bP.geoCenter_x - 1;
			int row = cl.bP.geoCenter_y;

			int inity = row - TMath::Floor( cl.bP.width_x/2.);

			// Check for zeroes first
			bool containsZeroes = false;
			for (int yi = inity ; yi < inity + m_selSizeY ; yi++) {
				if ( GetCalibEnergy(column, yi) == 0. ) {
					containsZeroes = true;
					break;
				}
			}
			if(containsZeroes) break;

			int itr = 0;
			for (int yi = inity ; yi < inity + m_selSizeY ; yi++) {
				m_energyAverageY[itr].push_back( GetCalibEnergy(column, yi) );
				Log << MSG::DEBUG << "InColumn | " << column << "," << yi << " | " << GetCalibEnergy(column, yi) << endreq;
				itr++;
			}
			m_energyAverageYCntr++;

		}

		cntr++;
	}

	// See if I want to skip visualization
	if ( skipFlag ) {
		Signals * sig = new Signals(this, Signals::__SIGNAL_SKIP);
		PullToStoreGateAccess(sig, MPXDefs::DO_NOT_SERIALIZE_ME);
	}

	//m_extraCanvas.push_back( DrawInSeparateWindow(m_graph2DVector, MSG::DEBUG) );

	// Fill the output tree of this algorithm
	getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_clusterEnergy.clear();
	m_clusterTOT.clear();
	m_driftingTime1.clear();
	m_clusterSizeX.clear();
	m_clusterSizeY.clear();

}

void Timepix3Preliminary::Finalize() {

	Log << MSG::INFO << "Finalize function !" << endreq;

	// Calculate the averages for each bin
	vector<double> energyAverageY;
	vector<double> energyAverageY_err;

	for(int i = 0 ; i < m_selSizeY ; i++) {

		int binCount = (int) m_energyAverageY[i].size();
		Log << MSG::INFO << "Average in bin  " << i << " with " << binCount << " entries " << endreq;
		MAFTools::DumpVector( m_energyAverageY[i] );
		double mean = MAFTools::CalcMean( m_energyAverageY[i] );

		// As a measurement of the dispersion here I will consider only the point within 1 stdev from the mean
		double stdv_sel = MAFTools::CalcStdDev( m_energyAverageY[i] );
		// select them
		vector<double> subset;
		vector<double>::iterator itr = m_energyAverageY[i].begin();
		for( ; itr != m_energyAverageY[i].end() ; itr++){
			if ( fabs((*itr) - mean) < stdv_sel ) {
				subset.push_back(*itr);
			}
		}
		// Now I'll use the stdv of the subset
		double stdv = MAFTools::CalcStdDev( subset ) / 2.; // divided by 2. 'cause TGraph will plot +/- this value
 		Log << MSG::INFO << "mean = " << mean << " | stddev = " << stdv << endreq;

		energyAverageY.push_back( mean );
		energyAverageY_err.push_back( stdv );

	}

	// X axis
	vector<double> xg(m_selSizeY, 0);
	vector<double> xg_err(m_selSizeY, 0);

	int sfrac = TMath::Floor(m_selSizeY/2.);
	int cntri = 0;
	for (int i = -1*sfrac ; i < sfrac ; i++) {
		xg[cntri] = i;
		xg_err[cntri] = 0;
		cntri++;
	}
	if( m_selSizeY % 2 != 0) { // if m_selSizeY is odd there's one entry missing
		xg[cntri] = sfrac;
		xg_err[cntri] = 0;
	}

	TGraphErrors * gy = new TGraphErrors(m_selSizeY, &xg[0], &energyAverageY[0], &xg_err[0], &energyAverageY_err[0] );
	getMyROOTFile()->cd();
	gy->Write("energyAverageY");

	int ainit[] = { 125, 50, 100, 101, 75 };
	vector<int> a ( ainit, ainit + sizeof(ainit)/sizeof(ainit[0]) );

	Log << MSG::ALWAYS << "mean  = " << MAFTools::CalcMean(a) << " | stddev = " << MAFTools::CalcStdDev(a) << endreq;

}

#endif
