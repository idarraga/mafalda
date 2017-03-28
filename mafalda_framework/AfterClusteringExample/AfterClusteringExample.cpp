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

/*
 *  An example of how to use the Clustering results for further processing.
 */

#ifndef AfterClusteringExample_cpp
#define AfterClusteringExample_cpp

#include "AfterClusteringExample.h"
#include "MAFTools/MAFTools.h"

using namespace MSG;

ClassImp(AfterClusteringExample)

AfterClusteringExample::AfterClusteringExample(){

}

void AfterClusteringExample::Init(){

	Log << MSG::INFO << "Init function !" << endreq;

	// This value will be overridden by the configuration since it'll be set up
	//  a few lines below as a configuration value
	m_minNPixels = 5;

	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("perpdistance", &m_perpdistance);

	// A configuration value that can be tuned from the Viewer
	RegisterConfigurationValue(&m_minNPixels, "minNPixels");

}

void AfterClusteringExample::Execute(){

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

	// A few variables needed for the arithmetics in this example
	Float_t cut = 0., slope = 0.;
	Float_t pixx = 0., pixy = 0.;
	Float_t cutPrim = 0, slopePrim = 0.;
	Float_t crossx = 0., crossy = 0.;

	for( ; blobsItr != blobsVector.end() ; blobsItr++)
	{

		// Limit all this to clusters with a minimum size.
		// Note that m_minNPixels can be configured through the Viewer
		//  so you can reprocess and check results online.
		if((*blobsItr).bP.nPixels < m_minNPixels)
			continue;

		// BlobsFinder has pre-calculated a few characteristics on each cluster.
		// If you want to use those values you can get them here.  For instance there is a LinearRegresion
		//  off all pixels belonging to the Cluster.  Here's the slope and ycut values.
		cut = (*blobsItr).bP.fitCut;
		slope = (*blobsItr).bP.fitSlope;
		cutPrim = 0.;
		slopePrim = -1./slope;

		// If you want to look at this cluster in detail, i.e. pixel per pixel, you can get the contents
		//  of the cluster in the form of a list of pairs [(X,Y), counts]
		list< pair < pair<Int_t, Int_t>, Int_t > >::iterator listItr;

		for(listItr = (*blobsItr).blobContent.begin() ; listItr != (*blobsItr).blobContent.end() ; listItr++)
		{

			// As an example I will calculate the perpendicular distance between each point in the cluster
			//  and the line described by the Linear Regression. (this could be useful to identify MIPs
			//  for instance).

			pixx = (Float_t)(*listItr).first.first;
			pixy = (Float_t)(*listItr).first.second;

			cutPrim  = pixy;
			cutPrim -= slopePrim * pixx;

			crossx = (cut - cutPrim)/(slopePrim - slope);
			crossy = crossx*slope + cut;

			// Calculate the perpendicular distance between a point and the line following the
			// linear regression, and push back to the results vector.  Note the function from MAFTools !
			m_perpdistance.push_back(  MAFTools::CalcDistance(crossx, crossy, pixx, pixy)  );

			// You can naturally get the TOT and MC information if available
			//  Also lvl1 information is available for those chips who have it.
			//  See MediPixAlgo::GetLVL1(...)
			Log << MSG::DEBUG << "TOT (" << pixx << ", " << pixy << ") = "
					<< GetMatrixElement(pixx, pixy) << endreq ;
			if ( IsMCData() ) {
				Log << MSG::DEBUG << "MC energy = " << GetMatrixElementMCEdep(pixx, pixy) << endreq;
			}

		}

	}

	// Fill the output tree of this algorithm
	getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_perpdistance.clear();

}

void AfterClusteringExample::Finalize() {

	Log << MSG::INFO << "Finalize function !" << endreq;

}

#endif
