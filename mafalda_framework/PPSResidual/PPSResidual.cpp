/*
 * 	Copyright 2010 John Idarraga, Mathieu Benoit
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


#ifndef PPSResidual_cpp
#define PPSResidual_cpp

#include "PPSResidual.h"
#include "BlobsFinder/BlobsFinder.h"
#include "MPXAlgo/Highlighter.h"
#include "MAFTools/MAFTools.h"
#include "TVector3.h"
#include "TF1.h"

#include <vector>

using namespace std;
using namespace MSG;

ClassImp(PPSResidual)

PPSResidual::PPSResidual(){


}

void PPSResidual::Init(){

	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	//getMyTree()->Branch("myhits", &myhits , "totalHits/I:totalCharge");

	m_pixelSizeX = 0.4; // mm FE-I3
	m_pixelSizeY = 0.05; // mm FE-I3
	m_nClusters = 0;

	getMyTree()->Branch("deviation_single", &m_deviation_sh);
	getMyTree()->Branch("deviation_double", &m_deviation_dh);
	getMyTree()->Branch("vertex", &m_vertex.x, "vertexX/D:vertexY");
	getMyTree()->Branch("nClusters", &m_nClusters, "nClusters/I");

	RegisterConfigurationValue(&m_pixelSizeX, "pixelSizeX");
	RegisterConfigurationValue(&m_pixelSizeY, "pixelSizeY");
}

void PPSResidual::Execute(){

	// REMINDER: PPSResidual::Execute() runs once per frame
	//  you may need ro reinitialize variables.

	Float_t mwidth = (Float_t)GetMatrixWidth();
	Float_t mheight = (Float_t)GetMatrixHeight();

	// Iterate over blobs to find MIP tracks
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");
	m_allBlobs = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject - 1);

	vector<blob> blobs = m_allBlobs->GetBlobsVector();
	vector<blob>::iterator blobsItr = blobs.begin();

	set< pair<Int_t, Int_t > > content;
	set< pair<Int_t, Int_t > >::iterator itrContent;

	for( ; blobsItr != blobs.end() ; blobsItr++)
	{


		blobtype bt = (*blobsItr).GetBlobType();
		Log << MSG::INFO << "blobtype is: " << bt << endreq;
		Log << MSG::INFO << "nvertex : " << GetPrimaryMCVertex_N() << endreq;

		Double_t Vx = GetPrimaryMCVertex_X(0);
		Double_t Vy = GetPrimaryMCVertex_Y(0);
		Double_t Bx = 0., By = 0., dist = 0.;
		m_vertex.x = Vx;
		m_vertex.y = Vy;

		if(bt == _SINGLE_HIT){

			Bx = (*blobsItr).GetBlobProperties().geoCenter_x * m_pixelSizeX;
			Bx -= ( ( mwidth/2.) * m_pixelSizeX ) ;
			Bx *= -1.;

			By = (*blobsItr).GetBlobProperties().geoCenter_y * m_pixelSizeY;
			By -= ( ( mheight/2.) * m_pixelSizeY );
			By *= -1.;

			dist = MAFTools::CalcDistance(Vx, Vy, Bx, By);
			m_deviation_sh.push_back(dist);

			Log << MSG::INFO << "Distance single hit = " << dist << endreq;


		}else if(bt == _DOUBLE_HIT){

			Bx = (*blobsItr).GetBlobProperties().weightedCenter_x * m_pixelSizeX;
			Bx -= ( ( mwidth/2.) * m_pixelSizeX ) ;
			Bx *= -1.;

			By = (*blobsItr).GetBlobProperties().weightedCenter_y * m_pixelSizeY;
			By -= ( ( mheight/2.) * m_pixelSizeY );
			By *= -1.;

			dist = MAFTools::CalcDistance(Vx,Vy,Bx,By);
			m_deviation_dh.push_back(dist);

			Log << MSG::INFO << "Distance double hit = " << dist << endreq;

		}

		m_nClusters = (Int_t)blobs.size();
		// Fill and clean
		getMyTree()->Fill();
	}



	m_deviation_sh.clear();
	m_deviation_dh.clear();
	m_vertex.x = 0;
	m_vertex.y = 0;
	m_nClusters = 0;
}


void PPSResidual::Finalize(){

	Log << MSG::INFO << "Finalize function !" << endreq;

}

#endif
