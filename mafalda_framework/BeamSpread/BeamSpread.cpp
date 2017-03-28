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

#ifndef BeamSpread_cpp
#define BeamSpread_cpp

#include "BeamSpread.h"
#include "BlobsFinder/BlobsFinder.h"
#include "MAFTools.h"

//#include "TRandom3.h"

using namespace MSG;

ClassImp(BeamSpread)

BeamSpread::BeamSpread(){

}

void BeamSpread::Init(){

	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("spreadDist", &m_spreadDist , "spreadDist/D");

	RegisterConfigurationValue(&m_pixelSizeX, "pixelSizeX");
	RegisterConfigurationValue(&m_pixelSizeY, "pixelSizeY");

}

void BeamSpread::Execute(){

	// Iterate over blobs to find MIP tracks
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");
	m_allBlobs = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject - 1);

	vector<blob> blobs = m_allBlobs->GetBlobsVector();
	vector<blob>::iterator itrBlob = blobs.begin();

	set< pair<Int_t, Int_t > > content;
	set< pair<Int_t, Int_t > >::iterator itrContent;


	if(itrBlob == blobs.end()){ // no blobs
		Log << MSG::WARNING << "No blobs !, skip event" << endreq;
		return;
	}


	int mwidth = GetMatrixWidth();
	int mheight = GetMatrixHeight();
	double centx = 0.;
	double centy = 0.;

	int nVertex = GetPrimaryMCVertex_N();

	for(int pvItr = 0 ; pvItr < nVertex ; pvItr++){

		double vx = GetPrimaryMCVertex_X(pvItr);
		double vy = GetPrimaryMCVertex_Y(pvItr);
		double vz = GetPrimaryMCVertex_Z(pvItr);

		m_spreadDist = 0.;

		for( ; itrBlob != blobs.end() ; itrBlob++)
		{

			// bring it to coordinates 0,0 at the center of the pixel pad
			//centx = (*itrBlob).GetBlobProperties().geoCenter_x*(3.6/(mwidth/2.)) + -3.6; // [mm]
			//centy = (*itrBlob).GetBlobProperties().geoCenter_y*(4.1/(mheight/2.)) + -4.1;

			centx = (double) ( (*itrBlob).GetBlobProperties().geoCenter_x - (mwidth/2) );
			centy = (double) ( (*itrBlob).GetBlobProperties().geoCenter_y - (mheight/2) );
			centx *= m_pixelSizeX; // mm
			centy *= m_pixelSizeY; // mm
			centx *= -1; // turned 180 deg
			centy *= -1;

			Log << MSG::DEBUG << "vx    : " << vx << endreq;
			Log << MSG::DEBUG << "vy    : " << vy << endreq;
			Log << MSG::DEBUG << "vz    : " << vz << endreq;
			Log << MSG::DEBUG << "centx : " << centx << endreq;
			Log << MSG::DEBUG << "centy : " << centy << endreq;

			// Calculate distance
			double newdist = MAFTools::CalcDistance(centx, centy, vx, vy);
			Log << MSG::DEBUG << "dist  : " << newdist << endreq;
			if(newdist > m_spreadDist){
				m_spreadDist = newdist;
			}

		}

	}

	getMyTree()->Fill();

}

void BeamSpread::Finalize(){


}

#endif
