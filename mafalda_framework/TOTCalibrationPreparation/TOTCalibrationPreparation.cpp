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

#ifndef __TOTCalibrationPreparation_cpp
#define __TOTCalibrationPreparation_cpp

#include "TOTCalibrationPreparation.h"
#include "MAFTools.h"

using namespace MSG;

ClassImp(TOTCalibrationPreparation)

TOTCalibrationPreparation::TOTCalibrationPreparation(){

	m_sizex = 256;
	m_sizey = 256;

	m_SingleHit_EntriesHisto = 0x0;
	m_DoubleHit_EntriesHisto = 0x0;
	m_TripleHit_EntriesHisto = 0x0;
	m_QuadHit_EntriesHisto = 0x0;


}

void TOTCalibrationPreparation::Init(){

	getMyTree()->Branch("m_SingleHitCoor", &m_SingleHitCoor);
	getMyTree()->Branch("m_SingleHitTOT", &m_SingleHitTOT);
	getMyTree()->Branch("m_DoubleHitCoor", &m_DoubleHitCoor);
	getMyTree()->Branch("m_DoubleHitTOT", &m_DoubleHitTOT);

	// Global TOT distributions
	m_SingleHit_GlobalTOTDist = new TH1F("SingleHitTOTDist","SingleHitTOTDist",100,0,100);

	// Clean up variables in init
	CleanUpVars();

}

void TOTCalibrationPreparation::Execute(){


	PrepareEntriesHistograms();

	// Before anything else check that the DACs and the clock are consistent over the entire run
	if(!CheckDACs()) {
		Log << MSG::ERROR << "Some of the DACs or the MpxClock setting changed in this particular frame." << endreq;
		Log << MSG::ERROR << "This is a calibration run and the DACs need to be consistent over the run." << endreq;
		Log << MSG::ERROR << "Aborting run.  Contact John Idarraga <idarraga@cern.ch>" << endreq;
		exit(1);
	}

	m_frameId = GetFrameId();

	int lastObject = GetNumberOfObjectsWithAuthor("SimpleClustering");
	if(lastObject == 0)
		return;

	m_sC = static_cast<SimpleClusterContainer *> ( GetObjectFromAuthor("SimpleClustering", lastObject-1) );
	vector<pair<int, int> > sh = m_sC->GetSingleHits();
	vector<pair<int, int> >::iterator shItr = sh.begin();

	if(m_frameId % 1000 == 0) { // print only every 1000 frames
		Log << MSG::INFO << "Frame " << m_frameId << ". Number of clusters = " << (int) sh.size() << endreq;
	}

	int framewidth = GetWidth();
	int tempTOT;

	for ( ; shItr != sh.end() ; shItr++) {

		// Store coordinate and TOT to Ntuple
		// <X> and <TOT>
		tempTOT = GetMatrixElement( *shItr );
		m_SingleHitCoor.push_back( MAFTools::XYtoX( *shItr , framewidth ) );
		m_SingleHitTOT.push_back( tempTOT );

		// Stats plots. Entries and not TOT.
		m_SingleHit_EntriesHisto->Fill( (*shItr).first, (*shItr).second ); // increase by one

		// Global TOT distribution
		m_SingleHit_GlobalTOTDist->Fill( tempTOT );

	}


	/*

	return;

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
	if(m_frameId % 1000 == 0) { // print only every 1000 frames
		Log << MSG::INFO << "Frame " << m_frameId << ". Number of clusters = " << (Int_t) blobsVector.size() << endreq;
	}
	vector<blob>::iterator blobsItr = blobsVector.begin(); //allBlobs.begin();

	//int framewidth = GetWidth();

	cluster c;
	list< pair < pair<int, int>, int > > cdes; // Cluster description
	list< pair < pair<int, int>, int > >::iterator hot; // Extra iterator for hotest pixel in the case of clusterSize > 1
	pair<int, int> pix;
	// Considering single, double, triple and quad clusters only, as follows :
	// Single :  *
	// Double :  ** || *  || *  ||  *
	//                 *      *    *
	// Triple :  ** || ** || *  ||  *
	//            *    *     **    **
	// Quad   :  **
	//           **
	for( ; blobsItr != blobsVector.end() ; blobsItr++) {

		c = *blobsItr;
		cdes = c.GetClusterDescription();

		if (c.btype == _SINGLE_HIT) {         // process single hits

			// Store coordinate and TOT to Ntuple
			// <X> and <TOT>
			m_SingleHitCoor.push_back( MAFTools::XYtoX( (*cdes.begin()).first , framewidth ) );
			m_SingleHitTOT.push_back( (*cdes.begin()).second );

			//Log << MSG::INFO << "Single --> " << (*cdes.begin()).first.first << "," << (*cdes.begin()).first.second << endreq;

			// Stats plots. Entries and not TOT.
			pix = (*cdes.begin()).first;
			m_SingleHit_EntriesHisto->Fill( pix.first, pix.second ); // increase by one

			// Global TOT distribution
			m_SingleHit_GlobalTOTDist->Fill(c.bP.clusterTOT);

		} else if (c.btype == _DOUBLE_HIT) { // process double hits

			// Associate the whole charge to the hotest pixel
			hot = FindHotestPixel( &cdes );
			m_DoubleHitCoor.push_back( MAFTools::XYtoX( (*hot).first , framewidth ) );
			m_DoubleHitTOT.push_back( c.bP.totalCharge );

			m_DoubleHit_EntriesHisto->Fill( (*hot).first.first, (*hot).first.second );

		} else if (c.btype == _TRIPLE_HIT) { // process triple hits

			// Associate the whole charge to the hotest pixel
			hot = FindHotestPixel( &cdes );
			//m_TripleHitMap[ (*hot).first ].push_back( c.bP.totalCharge );
			m_TripleHit_EntriesHisto->Fill( (*hot).first.first, (*hot).first.second );

		} else if (c.btype == _QUAD_HIT) {  // process quad hits

			// Associate the whole charge to the hotest pixel
			hot = FindHotestPixel( &cdes );
			//m_QuadHitMap[ (*hot).first ].push_back( c.bP.totalCharge );
			m_QuadHit_EntriesHisto->Fill( (*hot).first.first, (*hot).first.second );

		} else {
			// Other than single, double, triple, quad, are clusters not considered in calibration
		}

	}

	 */

	Log << MSG::DEBUG << m_SingleHitCoor.size() << " single hits" << endreq;

	// Fill the output tree of this algorithm
	getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	CleanUpVars();

}

bool TOTCalibrationPreparation::CheckDACs(){

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
			if( dacs[i] != m_dacs[i] ) dacsok = false; // if any of the dacs changed
		}

		if( !dacsok ) {
			return false;
		}
	}

	// Check also the size of the matrix
	vector<Float_t> metaInfo;
	metaInfo.push_back( GetMatrixWidth() );
	metaInfo.push_back( GetMatrixHeight() );
	m_MetaData = new TVectorF ( (int)metaInfo.size(), &metaInfo[0] );
	m_MetaData->Write("MetaData");

	return true;
}

void TOTCalibrationPreparation::CleanUpVars(){

	// single hit
	m_SingleHitCoor.clear();
	m_SingleHitTOT.clear();

}

void TOTCalibrationPreparation::Finalize() {

	// Stats plots
	//BuildStatsPlots();
	m_SingleHit_EntriesHisto->Write();
	m_DoubleHit_EntriesHisto->Write();
	m_TripleHit_EntriesHisto->Write();
	m_QuadHit_EntriesHisto->Write();

	/*

	// Build full spectrum per pixel
	TString histoName;

	map<pair<int,int>, vector<double> >::iterator itr = m_SingleHitMap.begin();
	vector<double> * criticalPoints;
	int nPeaks = 0;

	for( ; itr != m_SingleHitMap.end() ; itr++) {

		if((*itr).second.size() > 10) { // at least 10 hits

			int x = (*itr).first.first;
			int y = (*itr).first.second;

			Log << MSG::INFO << "Cleared for pixel " << x << "," << y << endreq;

			histoName = "spec_single_";
			histoName += x;
			histoName += "_";
			histoName += y;
			TH1D * h = new TH1D(histoName, histoName, 100, 0, 100); // 100TOT ~ 200keV
			h->GetBin(kBlue);

			// selected pixel
			pair<int,int> pix = make_pair(x, y);
			vector<double> pixinfo = m_SingleHitMap[pix];
			vector<double>::iterator i = pixinfo.begin();
			for( ; i != pixinfo.end() ; i++){
				h->Fill( *i );
			}
			h->Write();

			// Operation on the distribution
			//nPeaks = MAFTools::IdendifyPeaks_A1( h , criticalPoints );
			Log << MSG::DEBUG << "nEntries : " << h->GetEntries() << endreq;

		}

	}

	 */

	// Global TOT distributions
	m_SingleHit_GlobalTOTDist->Write();

}

void TOTCalibrationPreparation::FillStatsPlots(){

	// Build Entries maps as TH2I hitograms too look at the overall occupancy
	// Each pixel should have a sensible amount of entries to build a spectrum
	// per pixel.


	//map<pair<int,int>, vector<double> >::iterator i;

	// Build entries Single
	/*
	vector<int>::iterator i = m_SingleHitTOT.begin();
	vector<int>::iterator coor = m_SingleHitCoor.begin();
	pair<int, int> pix;

		for ( ; i != m_SingleHitTOT.end() ; ) {

		pix = MAFTools::XtoXY( *coor , sizex );

		cout << " : " << pix.first << " , " << pix.second << endl;
		m_SingleHit_EntriesHisto->Fill( pix.first, pix.second, *i );
		coor++; i++; // iterators

	}
	 */

	/*
	m_SingleHit_EntriesHisto = new TH2I("entriesSingle","entriesSingle", sizex, 0, sizex, sizey, 0, sizey);
	map<pair<int,int>, vector<double> >::iterator i = m_SingleHitMap.begin();
	for ( ; i != m_SingleHitMap.end() ; i++) {
		m_SingleHit_EntriesHisto->Fill( (*i).first.first, (*i).first.second, (*i).second.size() );
	}
	m_SingleHit_EntriesHisto->Write();
	 */

	/*
	// Build entries Double
	m_DoubleHit_EntriesHisto = new TH2I("entriesDouble","entriesDouble", sizex, 0, sizex, sizey, 0, sizey);
	for (i = m_DoubleHitMap.begin() ; i != m_DoubleHitMap.end() ; i++) {
		m_DoubleHit_EntriesHisto->Fill( (*i).first.first, (*i).first.second, (*i).second.size() );
	}
	m_DoubleHit_EntriesHisto->Write();

	// Build entries Triple
	m_TripleHit_EntriesHisto = new TH2I("entriesTriple","entriesTriple", sizex, 0, sizex, sizey, 0, sizey);
	for (i = m_TripleHitMap.begin() ; i != m_TripleHitMap.end() ; i++) {
		m_TripleHit_EntriesHisto->Fill( (*i).first.first, (*i).first.second, (*i).second.size() );
	}
	m_TripleHit_EntriesHisto->Write();

	// Build entries Quad
	m_QuadHit_EntriesHisto = new TH2I("entriesQuad","entriesQuad",sizex, 0, sizex, sizey, 0, sizey);
	for (i = m_QuadHitMap.begin() ; i != m_QuadHitMap.end() ; i++) {
		m_QuadHit_EntriesHisto->Fill( (*i).first.first, (*i).first.second, (*i).second.size() );
	}
	m_QuadHit_EntriesHisto->Write();

	 */

}

list< pair < pair<int, int>, int > >::iterator TOTCalibrationPreparation::FindHotestPixel(list< pair < pair<int, int>, int > > * cdes){

	list< pair < pair<int, int>, int > >::iterator i = cdes->begin();
	list< pair < pair<int, int>, int > >::iterator hot;

	int max = 0.;
	for ( ; i != cdes->end() ; i++){
		if ( (*i).second > max ) {
			max = (*i).second;
			hot = i;
		}
	}

	// return the iterator to hottest pixel
	return hot;
}

void TOTCalibrationPreparation::PrepareEntriesHistograms(){


	// Prepare entries 2D histos

	if ( ! m_SingleHit_EntriesHisto ) {
		m_sizex = GetMatrixWidth();
		m_sizey = GetMatrixHeight();
		m_SingleHit_EntriesHisto = new TH2I("entriesSingle", "entriesSingle", m_sizex, 0, m_sizex, m_sizey, 0, m_sizey);
	}

	if ( ! m_DoubleHit_EntriesHisto ) {
		m_DoubleHit_EntriesHisto = new TH2I("entriesDouble", "entriesDouble", m_sizex, 0, m_sizex, m_sizey, 0, m_sizey);
	}

	if( ! m_TripleHit_EntriesHisto ) {
		m_TripleHit_EntriesHisto = new TH2I("entriesTriple", "entriesTriple", m_sizex, 0, m_sizex, m_sizey, 0, m_sizey);
	}

	if ( ! m_QuadHit_EntriesHisto ) {
		m_QuadHit_EntriesHisto = new TH2I("entriesQuad", "entriesQuad", m_sizex, 0, m_sizex, m_sizey, 0, m_sizey);
	}

	//int binsx = m_SingleHit_EntriesHisto->GetNbinsX();
	//int binsy = m_SingleHit_EntriesHisto->GetNbinsY();

	//for(int bx = 1 ; bx <= binsx ; bx++){
	//for(int by = 1 ; by <= binsy ; by++){
	//m_SingleHit_EntriesHisto->SetBinContent(bx, by, 0);
	//}
	//}

}

//template class TVectorT<int>;

#endif



