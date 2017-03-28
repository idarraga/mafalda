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

#ifndef __THLCalibration_cpp
#define __THLCalibration_cpp

#include "THLCalibration.h"
#include "MAFTools.h"

using namespace MSG;

ClassImp(THLCalibration)

THLCalibration::THLCalibration() : MediPixAlgo(), CalibrationLoader(this) {

	// This value will be overridden by the configuration because it'll registered
	//  as a ConfigurationValue in the Init member of this class.
	m_frameGlobalId = 0;
	m_nSmoothings = 3;

}

void THLCalibration::Init(){

	Log << MSG::INFO << "Init function !" << endreq;


	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("clusterTOT", &m_clusterTOT);
	getMyTree()->Branch("clusterEnergy", &m_clusterEnergy);

	// A configuration value that can be tuned from the Viewer
	RegisterConfigurationValue(&m_nSmoothings, "nSmoothings");

}

void THLCalibration::Execute(){

	// Before anything else check that the DACs and the clock are consistent over the entire run
	if ( !CheckDACs() ) {
		Log << MSG::ERROR << "Some of the DACs or the MpxClock setting changed in this particular frame." << endreq;
		Log << MSG::ERROR << "This is a calibration run and the DACs need to be consistent over the run." << endreq;
		Log << MSG::ERROR << "Aborting run.  Contact John Idarraga <idarraga@cern.ch>" << endreq;
		exit(1);
	}

	m_frameId = GetFrameId();
	m_frameGlobalId = GetGlobalCounters();

	// Get the THL
	vector<int> dacs = GetDAQs();
	int thl = dacs[6];

	//int lastObject = GetNumberOfObjectsWithAuthor("SimpleClustering");
	int lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");

	int nsh = 0;

	if ( lastObject != 0 ) {

		//m_aB = static_cast<SimpleClusterContainer *> ( GetObjectFromAuthor("SimpleClustering", lastObject-1) );
		//vector<pair<int, int> > sh = m_aB->GetSingleHits();
		//vector<pair<int, int> >::iterator itr = sh.begin();
		//int nsh = (int) sh.size();

		m_aB = static_cast<AllBlobsContainer *> ( GetObjectFromAuthor("BlobsFinder", lastObject-1) );
		vector<blob> sh = m_aB->GetBlobsVector();

		// Use the number of clusters as counting
		// It can be used when the clusters are well separated.  After a certain occupancy this is not a good approach.
		// The Counting will artificially drop at low THL values
		nsh = (int) sh.size();

		if(m_frameId % 100 == 0) {
			Log << MSG::INFO << "Frame " << m_frameGlobalId << "/" << GetNFrames() << " | THL = " << thl << ". Number of single hit clusters = " << nsh << endreq;
		}


	} else {

		nsh = GetHitsInPad();

		if(m_frameId % 100 == 0) {
			Log << MSG::INFO << "Frame " << m_frameGlobalId << "/" << GetNFrames() << " | THL = " << thl << ". Number of hits (not clusters) = " << nsh << endreq;
		}

	}

	//nsh = GetHitsInPad();

	// Use the number of hits
	// Only good for low energies and necesary for high occupancies


	//	// Check is this cluster has a constituent in the extended mask
	//	vector<pair<int, int> >::iterator itr = sh.begin();
	//	if ( UsingExtendedMask() ) {
	//
	//		for ( ; itr != sh.end() ; itr++ ) {
	//			// in this case reduce the number of single hits by one
	//			if( PixelInExtendedMask( *itr, GetWidth() ) ) {
	//				//Log << MSG::INFO << " ++++++++++++++++++++++++++++++++++++++++++++++++++++" << endreq;
	//				nsh--;
	//			}
	//		}
	//
	//	}


	// Push back the number of single hits per second
	double acqtime = GetAcqTime();
	m_THLCountsMap[thl].push_back( nsh / acqtime );

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_clusterEnergy.clear();
	m_clusterTOT.clear();

}

#include <time.h>

void THLCalibration::Finalize() {

	Log << MSG::INFO << "Finalize function !" << endreq;

	// Make the plot of number of SingleHits cluster as a function of THL
	map<int, vector<int> >::iterator i = m_THLCountsMap.begin();
	vector<int>::iterator iv;

	double nshMean = 0.; // number of single hits mean
	double nshSD = 0.; // number of single hits StandardDeviation
	int thl = 0;
	vector<double> x;
	vector<double> x_err;
	vector<double> y;
	vector<double> y_err;
	int N = (int) m_THLCountsMap.size();

	TF1 * f1 = new TF1("f1","gaus",0,500);
	f1->SetParameter(0, 1000);
	f1->SetParameter(1, 250);
	f1->SetParameter(2, 50);

	// Seed the random number generator manually
	time_t rawtime;
	time(&rawtime);
	long myseed = long(rawtime);
	Log << MSG::INFO << "The random seed (localtime): " << myseed << endreq;
	TRandom1 * rand = new TRandom1(myseed);

	for( ; i != m_THLCountsMap.end() ; i++ ) {

		thl = (*i).first;
		nshMean = MAFTools::CalcMean( (*i).second );
		nshSD = MAFTools::CalcStdDev( (*i).second );

		x.push_back( thl );
		x_err.push_back( 0. );
		y.push_back( nshMean );
		y_err.push_back( nshSD );

		// fake a gaussian
		//x.push_back( thl );
		//x_err.push_back( 0. );
		//y.push_back( f1->Eval(thl) + 50*(rand->Rndm()*2 - 1.) );
		//y_err.push_back( 0. );

	}

	vector<TGraph *> graphs;
	m_g_THLCounts = new TGraphErrors( (int)x.size(), &x[0], &y[0], &x_err[0], &y_err[0] );
	graphs.push_back( m_g_THLCounts );
	m_g_THLCounts->SetName("CountsTHL");

	// Smooth the function.  Applying triangular smoothing.
	//TGraphErrors * gs1 = MAFTools::TriangularSmooth(m_g_THLCounts);
	//graphs.push_back( gs1 );

	// keep smoothing a few times
	TGraphErrors * gsmoothed = m_g_THLCounts;
	for ( int sm = 0 ; sm < m_nSmoothings ; sm++ ) {
		TGraphErrors * gnew = MAFTools::TriangularSmooth(gsmoothed);
		TString name = "CountsTHL_";
		name += sm;
		gnew->SetName(name);
		if ( gnew == 0x0 ) break; // this means the function could not smooth any more, working with last result
		graphs.push_back( gnew );
		gsmoothed = gnew;
	}

	// Derivative
	vector<TGraph *> graphfinal;
	TGraph * gder = MAFTools::GetDerivative( (TGraph *) gsmoothed );
	gder->SetMarkerStyle(7);
	gder->SetName("spectrum");
	graphfinal.push_back( gder );


	TCanvas * c1 = DrawInSeparateWindow( graphs , "CountsVsTHL");
	TCanvas * c2 = DrawInSeparateWindow( graphfinal , "d(Counts)");

	// Save all curves and the final derivative
	getMyROOTFile()->cd();
	vector<TGraph *>::iterator gi = graphs.begin();
	for( ; gi != graphs.end() ; gi++) (*gi)->Write();
	// and write the derivative
	graphfinal[0]->Write();

}

bool THLCalibration::CheckDACs(){

	vector<int> dacs = GetDAQs();
	int clock1 = GetMpxClock();
	// Append the clock at 14th element
	dacs.push_back( clock1 );

	// The very first time copy the dacs to a member TVectorF to store
	if( m_dacs.empty() ) {

		m_dacs = GetDAQs();
		m_dacs.push_back( clock1 );

		// TVectorT<int> can not be serialized by default.  Creating a vector<Float> to
		//  be used in the constructor of the serializable object
		vector<Float_t> vc;
		for(int i = 0 ; i < (int) m_dacs.size() ; i++) {
			vc.push_back( (Float_t) m_dacs[i] );
		}

		// Finally create the object to store
		m_DACsStore = new TVectorF( (int)vc.size(), &vc[0] );
		//m_DACsStore->
		m_DACsStore->Write("DACs");

	} else { // for the rest of the run make sure this doesn't change

		bool dacsok = true;
		for (int i = 0; i < (int)m_dacs.size() ; i++) {

			// For the THL calibration the THL is alloud to change
			if ( i != 6 ) {
				if( dacs[i] != m_dacs[i] ) dacsok = false; // if any of the dacs changed
			}

		}

		if( !dacsok ) {
			return false;
		}
	}

	return true;
}

#endif
