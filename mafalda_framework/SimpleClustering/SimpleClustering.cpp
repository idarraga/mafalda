/*
 * 	Copyright 2012 John Idarraga
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

#ifndef __SimpleClustering_cpp
#define __SimpleClustering_cpp

#include "SimpleClustering.h"

using namespace MSG;

ClassImp(SimpleClustering)
ClassImp(SimpleClusterContainer)

SimpleClustering::SimpleClustering(){

	m_maxOcc = 0.5;
}

void SimpleClustering::Init(){

	RegisterConfigurationValue(&m_maxOcc, "maxOcc");

}

void SimpleClustering::Execute(){

	// Check excessive occupancy
	if ( (double) GetHitsInPad() / double ( GetWidth() * GetHeight() ) > m_maxOcc ) {
		Log << MSG::DEBUG << "Frame " << GetFrameId() << " with too high occupancy, was excluded." << endreq;
		return;
	}

	int xDim = GetMatrixXdim();
	int yDim = GetMatrixYdim();
	int rowItr = 0;
	int colItr = 0;
	unsigned int word = 0;

	// create container for this event
	m_scc = new SimpleClusterContainer(this);

	/////////////////////////////////////////////////////////////////////////////////////
	// Scan the whole matrix without the border first.
	// Border and corners are special cases
	for(colItr = 1; colItr < xDim-1 ; colItr++) {

		for(rowItr = 1; rowItr < yDim-1 ; rowItr++) {

			// See if there's nothing to look for here
			if(GetMatrixElement(colItr, rowItr) == 0) continue;

			// make word
			word = GetMatrixCropAsBitWord(colItr-1, colItr+1, rowItr-1, rowItr+1);
			if (word == __inner_single_word) { // Single hit --> 0x10
				m_scc->InsertSingleHit(colItr, rowItr);
				//Log << MSG::INFO << "Single --> " << colItr << "," << rowItr << " : " << word << endreq;
				// if I find a single then I can jump twice to the right or up,
				// jump once here, and one in the loop
				rowItr++;
			}

		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	// Scan the borders avoiding the 4 corners (corners are only 4 pixels)
	// low border --> Single 0x02
	rowItr = 0;
	for (colItr = 1; colItr < xDim-1 ; colItr++) {
		if(GetMatrixElement(colItr, rowItr) == 0) continue;
		word = GetMatrixCropAsBitWord(colItr-1, colItr+1, rowItr, rowItr+1);
		if(word == __low_border_word){
			m_scc->InsertSingleHit(colItr, rowItr);
			colItr++;
			//Log << MSG::INFO << "Single low border --> " << colItr << "," << rowItr << " : " << word << endreq;
		}
	}
	// top border --> Single 0x10
	rowItr = yDim-1;
	for (colItr = 1; colItr < xDim-1 ; colItr++) {
		if(GetMatrixElement(colItr, rowItr) == 0) continue;
		word = GetMatrixCropAsBitWord(colItr-1, colItr+1, rowItr-1, rowItr);
		if(word == __top_border_word){
			m_scc->InsertSingleHit(colItr, rowItr);
			colItr++;
			//Log << MSG::INFO << "Single top border --> " << colItr << "," << rowItr << " : " << word << endreq;
		}
	}
	// left border --> Single 0x04
	colItr = 0;
	for (rowItr = 1; rowItr < yDim-1 ; rowItr++) {
		if(GetMatrixElement(colItr, rowItr) == 0) continue;
		word = GetMatrixCropAsBitWord(colItr, colItr+1, rowItr-1, rowItr+1);
		if(word == __left_border_word){
			m_scc->InsertSingleHit(colItr, rowItr);
			rowItr++;
			//Log << MSG::INFO << "Single left border --> " << colItr << "," << rowItr << " : " << word << endreq;
		}
	}
	// right border --> Single = 0x08
	colItr = xDim-1;
	for(rowItr = 1; rowItr < yDim-1 ; rowItr++) {
		if(GetMatrixElement(colItr, rowItr) == 0) continue;
		word = GetMatrixCropAsBitWord(colItr-1, colItr, rowItr-1, rowItr+1);
		if(word == __right_border_word){
			m_scc->InsertSingleHit(colItr, rowItr);
			rowItr++;
			//Log << MSG::INFO << "Single right border --> " << colItr << "," << rowItr << " : " << word << endreq;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	// Now take care of the 4 corners
	rowItr = 0; colItr = 0;
	word = GetMatrixCropAsBitWord(colItr, colItr+1, rowItr, rowItr+1);
	if(word == __blc_word){
		m_scc->InsertSingleHit(colItr, rowItr);
	}

	rowItr = 0; colItr = xDim-1;
	word = GetMatrixCropAsBitWord(colItr-1, colItr, rowItr, rowItr+1);
	if(word == __brc_word){
		m_scc->InsertSingleHit(colItr, rowItr);
	}

	rowItr = yDim-1; colItr = 0;
	word = GetMatrixCropAsBitWord(colItr, colItr+1, rowItr-1, rowItr);
	if(word == __ulc_word){
		m_scc->InsertSingleHit(colItr, rowItr);
	}

	rowItr = yDim-1; colItr = xDim-1;
	word = GetMatrixCropAsBitWord(colItr-1, colItr, rowItr-1, rowItr);
	if(word == __urc_word){
		m_scc->InsertSingleHit(colItr, rowItr);
	}


	Log << MSG::DEBUG << m_scc->GetNSingleHits() << " single hits found in frame " << GetFrameId() << endreq;

	// pull all blobs to store gage
	PullToStoreGateAccess(m_scc, MPXDefs::DO_NOT_SERIALIZE_ME);

	//getMyTree()->Fill();

}

void SimpleClustering::Finalize(){

}

SimpleClusterContainer::SimpleClusterContainer(MediPixAlgo * algo) :
				CandidateContainer(algo) {

}

#endif
