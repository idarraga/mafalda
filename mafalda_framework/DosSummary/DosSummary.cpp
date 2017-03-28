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

#ifndef __DosSummary_cpp
#define __DosSummary_cpp

#include "DosSummary.h"

using namespace MSG;

ClassImp(DosSummary)

DosSummary::DosSummary(){


}

void DosSummary::Init(){

	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("DosSummaryReport", & m_dos_summary_report);

	Log << MSG::INFO << "Init function !" << endreq;

}

void DosSummary::Execute(){

	// REMINDER: DosSummary::Execute() runs once per frame
	//  you may need ro reinitialize variables.

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

	cluster c;
	for( ; blobsItr != blobsVector.end() ; blobsItr++) {

		c = *blobsItr;
		list< pair < pair<int, int>, int > > cdes = c.GetClusterDescription();
		vector<pair<int, int> > pixels;
		vector<int> tot;
		vector<double> energy;
		list< pair < pair<int, int>, int > >::iterator i_cdes = cdes.begin();
		for( ; i_cdes != cdes.end() ; i_cdes++ ) {
			pixels.push_back( (*i_cdes).first ); // pixel coordinates (x,y)
			tot.push_back( (*i_cdes).second ) ;  // pixel tot
			energy.push_back( 0. );              // pixel calib energy
		}

		DosCluster dos_cluster;
		dos_cluster.SetClusterPixelComponents( pixels );
		dos_cluster.SetClusterTOTComponents( tot );
		dos_cluster.SetClusterCalibEnergyComponents( energy );

		dos_cluster.SetClusterDose( 0., "LET", DosCluster::__LET_BASED_DOSE);
		dos_cluster.SetClusterInnerPixels( c.bP.nInnerPixels );

		m_dos_summary_report.PushBackDosCluster( dos_cluster );

	}

	getMyTree()->Fill();

	// clean up
	m_dos_summary_report.Clear();

}

void DosSummary::Finalize(){

	Log << MSG::INFO << "Finalize function !" << endreq;

}

#endif
