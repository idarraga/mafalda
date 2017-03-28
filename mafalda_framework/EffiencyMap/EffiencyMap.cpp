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

#ifndef __EffiencyMap_cpp
#define __EffiencyMap_cpp

#include "EffiencyMap.h"
#include "MAFTools.h"

using namespace MSG;


ClassImp(EffiencyMap)

EffiencyMap::EffiencyMap(){

}

void EffiencyMap::Init(){

	m_maxNHits = 5;

	// A configuration value that can be tuned from the Viewer
	RegisterConfigurationValue(&m_maxNHits, "maxNHits");

}

void EffiencyMap::Execute(){

	Int_t width = GetWidth();
	Int_t height = GetHeight();

	if( ! m_firstTimeFlag ) {

		// Maps
		// hit map --> number of hits per pixel
		m_hitMap = new TH2F("hitmap", "hitmap",  GetWidth(), 0, GetWidth(),  GetHeight(), 0, GetHeight());
		// eq map --> mean tot per pixel
		m_equalizationMap = new TH2F("eqmap","eqmap",  GetWidth(), 0, GetWidth(),  GetHeight(), 0, GetHeight());

		// Create all the empty vectors per pixel to store the TOT
		for(Int_t xx = 0; xx < width-1 ; xx++) {
			for(Int_t yy = 0; yy < height-1 ; yy++) {
				vector<int> tmpV;
				m_perPixelTOT.push_back(tmpV);
			}
		}

		m_firstTimeFlag = true;
	}

	// Ask the store gate if the previous algorithm (BlobsFinder --> responsible for clustering)
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

	for( ; blobsItr != blobsVector.end() ; blobsItr++)
	{

		// Look for mips
		blobtype bT = (*blobsItr).GetBlobType();
		int nhits = GetHitsInPad();

		if ( bT == _SINGLE_HIT && nhits <= m_maxNHits ) { // only single hits && max hits cut

			set< pair<Int_t, Int_t> > bCSet = (*blobsItr).GetContentSet();
			set< pair<Int_t, Int_t> >::iterator bCItr = bCSet.begin();

			// clean up a bit more.  Get rid of noise columns
			//Log << MSG::INFO << "col : " << (*bCItr).first << endreq;
			//Log << MSG::INFO << "row : " << (*bCItr).second << endreq;

			vector<int> rowV = GetRowAsVector((*bCItr).second);
			vector<int> colV = GetColAsVector((*bCItr).first);

			// Calculate the number of hits in the row and column
			//  where the single hit was found
			int rowCntr = 0, colCntr = 0;
			vector<int>::iterator it = rowV.begin();
			for( ; it != rowV.end() ; it++){ if ((*it) > 0){rowCntr++;} }
			it = colV.begin();
			for( ; it != colV.end() ; it++){ if ((*it) > 0){colCntr++;} }

			//if ( (*bCItr).first == 31 || colCntr > 1 || rowCntr > 1 )
			//continue;

			m_hitMap->Fill( (*bCItr).first, (*bCItr).second, 1);
			m_equalizationMap->Fill( (*bCItr).first, (*bCItr).second, GetMatrixElement(*bCItr) );

			// landau pixel per pixel
			m_perPixelTOT[MAFTools::XYtoC((*bCItr), width)].push_back(GetMatrixElement(*bCItr));

		}

	}

	// Fill the output tree of this algorithm
	//getMyTree()->Fill();

}

void EffiencyMap::Finalize() {

	Log << MSG::INFO << "Finalize function !" << endreq;
	// histograms
	m_equalizationMap->Divide(m_hitMap);

	// write to file
	m_hitMap->Write();
	m_equalizationMap->Write();

}

#endif
