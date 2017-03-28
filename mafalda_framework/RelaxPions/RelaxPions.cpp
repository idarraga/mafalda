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

/**
 * Created automatically with MAFalda (Fri Nov 23 19:24:07 CET 2012)
 *
 * An example of how to use the clustering results
 * for further processing.
 */


#ifndef __RelaxPions_cpp
#define __RelaxPions_cpp

#include "RelaxPions.h"
#include "MAFTools/MAFTools.h"

using namespace MSG;

ClassImp(RelaxPions)

RelaxPions::RelaxPions() : MediPixAlgo(), CalibrationLoader(this) {

	// This value will be overridden by the configuration because it'll registered
	//  as a ConfigurationValue in the Init member of this class.
	m_minNPixels = 5;
	m_sizex = 512;
	m_sizey = 256;

}

void RelaxPions::Init(){

	Log << MSG::INFO << "Init function !" << endreq;


	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("clusterTOT", &m_clusterTOT);
	getMyTree()->Branch("clusterEnergy", &m_clusterEnergy);
	getMyTree()->Branch("clusterSize", &m_clusterSize);
	getMyTree()->Branch("clusterAngle", &m_clusterAngle);

	// A configuration value that can be tuned from the Viewer
	RegisterConfigurationValue(&m_minNPixels, "minNPixels");

	// Maps
	m_globalHitsEntries = new TH2I("GlobalHitsEntries", "GlobalHitsEntries", m_sizex, 0, m_sizex, m_sizey, 0, m_sizey);
	m_globalAverageTOT = new TH2D("GlobalAverageTOT", "GlobalAverageTOT", m_sizex, 0, m_sizex, m_sizey, 0, m_sizey);

}

void RelaxPions::Execute(){

	// output
	m_nGoodCandidates = 0;
	m_nGoodCandidatesFragmented = 0;
	m_GoodCadidatesClusterSize.clear();

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

	// Ready to skip the frame if nothing of relevance found
	bool skipViewer = true;

	cluster cl;
	for( ; blobsItr != blobsVector.end() ; blobsItr++) {

		cl = *blobsItr;

		// Limit all this to clusters with a minimum size.
		// Note that m_minNPixels can be configured through the Viewer
		//  so you can reprocess and check results online.
		if(cl.bP.nPixels < m_minNPixels)
			continue;


		if( cl.btype == _MIP && ! MAFTools::ClusterTouchesBorder( cl, m_sizex, m_sizey, 1 ) ) {

			// Cluster Size
			m_clusterSize.push_back( cl.bP.clusterSize );
			// Store the cluster TOT for output
			m_clusterTOT.push_back( cl.bP.clusterTOT );
			m_clusterAngle.push_back( cl.bP.rotAngle * 180 / TMath::Pi() );

			Log << MSG::INFO << "Mip found." << endreq;
			vector< pair<double, double> > limitPixels = MAFTools::FindLimitPixels_MaxMinXY(*blobsItr);
			Log << MSG::INFO << "Limits : "
					<< limitPixels[__HEAD_INDX].first << ","
					<< limitPixels[__HEAD_INDX].second <<
					" --> " << limitPixels[__TAIL_INDX].first
					<< "," << limitPixels[__TAIL_INDX].second
					<< endreq;

			// Search for a mip going through the interface of the two chips
			if (limitPixels[__HEAD_INDX].first < 256 && limitPixels[__TAIL_INDX].first > 256 ) {

				// Found good candidate !!!
				Log << MSG::INFO << "Good candidate !" << endreq;
				m_nGoodCandidates++;

				Highlighter * fp = new Highlighter(256,
						cl.bP.geoCenter_y,
						"circle", this);
				fp->SetLineColor(kRed);
				fp->SetLineWidth(1);
				PullToStoreGateAccess(fp, MPXDefs::DO_NOT_SERIALIZE_ME);

				// do not skip this frame
				skipViewer = false;

			}


			// With any identified MIP fill the Entries map using each hit information
			list< pair < pair<int, int>, int > > des = cl.GetClusterDescription();
			list< pair < pair<int, int>, int > >::iterator i = des.begin();
			int x,y,c,update_c;
			for( ; i != des.end() ; i++) {
				x = (*i).first.first;
				y = (*i).first.second;
				c = (*i).second;
				m_globalHitsEntries->Fill(x, y, 1);
				m_globalAverageTOT->Fill(x, y, c);
			}

		}

	}

	// see if I want to skip this frame
	if(skipViewer) {
		Signals * sig = new Signals(this, Signals::__SIGNAL_SKIP);
		PullToStoreGateAccess(sig, MPXDefs::DO_NOT_SERIALIZE_ME);
	}

	// Fill the output tree of this algorithm
	getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_clusterEnergy.clear();
	m_clusterTOT.clear();
	m_GoodCadidatesClusterSize.clear();
	m_clusterSize.clear();
	m_clusterAngle.clear();

}

void RelaxPions::Finalize() {

	Log << MSG::INFO << "Finalize function !" << endreq;

	GetMyROOTFile()->cd();

	// Hits map
	m_globalHitsEntries->Write();
	// Get the average
	double c, N;
	for (int j = 1 ; j <= m_sizey ; j++) {
		for (int i = 1 ; i <= m_sizex ; i++) {
			c = m_globalAverageTOT->GetBinContent(i,j);
			N = m_globalHitsEntries->GetBinContent(i,j);
			m_globalAverageTOT->SetBinContent(i,j, c/N);
		}
	}
	m_globalAverageTOT->Write();

}

#endif
