/*
 * 	Copyright 2009 John Idarraga
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

#ifndef PRBasicSpecies_cpp
#define PRBasicSpecies_cpp

#include "PRBasicSpecies.h"

using namespace MSG;

ClassImp(PRBasicSpecies)

/**
 * Constructor
 *
 */

PRBasicSpecies::PRBasicSpecies() : MediPixAlgo(), CalibrationLoader(this) {

	/* SG object */
	m_aB = 0;
	/* highlighters */
	m_highCurly = 0;
	m_highHeavyBlobs = 0;
	/* to ntuple */
	m_frameId = 0;
	/* other vars */
	m_overflowPixelsCntr = 0;
	// other swtiches
	m_checkOverFlow = false;
	// type switches
	m_HeavyBlobON = true;
	m_HeavyTrackON = true;
	m_mipON = true;

	/**
	 * Selection parameters for HeavyBlobs and HeavyTracks : This set of default parameters
	 * has been obtained by the study of separated and mixed Am241 and
	 * Cs137 data at the Universite de Montreal (2009).
	 */
	m_circleToEllipseMin = 1.4;
	m_circleToEllipseMax = 5.0;
	m_longGammaMax = 3;
	m_nInnerPixelsCut = 2;
	m_maxCurlyPixels = 50;

	// MIP
	m_mipChisquareDof = 0.5;
	m_mindistanceToLine = 1.0; // min distance 1 pixels (float)
	m_fractionOfPixelsAtMindistance = 0.95; // 95% of pixels
	m_mipNPixels = 25; // minimum number of pixels in a mip

	ClearNtupleVars();
}

void PRBasicSpecies::Init(){

	// a few branches for this algorithm's Tree
	getMyTree()->Branch("fId",&m_frameId,"fId/I");
	getMyTree()->Branch("speciesCntr", &m_basicSpeciesCntr.nSingleHits,
			"nSingleHits/I:nDoubleHits:nTripleHits:nQuadHits:nLongGamma:nHeavyBlobs:nHeavyTracks:nMip:nCurly:nNotBasicType:nAll");

	getMyTree()->Branch("generalCntr", &m_basicGeneral.nHits,
			"nHits/I:nCounts");

	/**  TODO ... remove the following branches
	 *  But I need to store some of this info.
	 *  May be all blob info, as it is in storegate should be dump
	 *  to a ROOT file
	 */
	getMyTree()->Branch("clusterSize", &m_clusterSize);
	getMyTree()->Branch("clusterSizeWidth", &m_clusterSizeWidth);
	getMyTree()->Branch("clusterSizeEllipse", &m_clusterSizeEllipse);
	getMyTree()->Branch("circleArea", &m_circleArea);
	getMyTree()->Branch("counts", &m_counts);
	/* ****** */

	RegisterConfigurationValue(&m_circleToEllipseMin, "circleToEllipseMin");
	RegisterConfigurationValue(&m_circleToEllipseMax, "circleToEllipseMax");
	RegisterConfigurationValue(&m_nInnerPixelsCut, "nInnerPixels");
	RegisterConfigurationValue(&m_longGammaMax, "longGammaMax");
	//RegisterConfigurationValue(&m_mipChisquareDof, "mipChisquareDof");
	RegisterConfigurationValue(&m_mindistanceToLine, "mindistanceToLine");
	RegisterConfigurationValue(&m_fractionOfPixelsAtMindistance,"fractionOfPixelsAtMindistance");
	RegisterConfigurationValue(&m_mipNPixels, "minNPixelsMip");
	RegisterConfigurationValue(&m_maxCurlyPixels, "maxNPixelsCurly");

}

void PRBasicSpecies::Execute(){

	// Pull out the object from "BlobsFinder" containing all found blobs and
	//  their characteristics.  It is the last object in the list (I need to know that)
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");
	if(lastObject == 0) return;  // nothing from clustering
	m_aB = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject-1);

	// Check clean frame // OFF BY DEFAULT
	if ( m_checkOverFlow ) // I better check this one first alone
		if ( FrameContainsOverFlowPixels() )
			return;

	// separate in basic species
	SeparateBasicSpecies();

	m_basicGeneral.nHits = GetHitsInPad();
	m_basicGeneral.nCounts = GetChargeInPad();

	// fill ntuple
	getMyTree()->Fill();

	// clear vars
	ClearNtupleVars();
	m_frameId++;
}

void PRBasicSpecies::SeparateBasicSpecies(){

	vector<blob> blobsVector = m_aB->GetBlobsVector();
	Log << MSG::DEBUG << "Number of blobs from previous algorithm = " << (Int_t) blobsVector.size() << endreq;

	vector<blob>::iterator blobsItr = blobsVector.begin(); //allBlobs.begin();

	double clusterEnergy = 0.;
	int clusterTOT = 0;
	double cE = 0.;
	pair<int, int> pix;
	bool badCalibInCluster;

	Int_t cntr = 0;
	for ( ; blobsItr != blobsVector.end() ; blobsItr++ ) {

		// Calibrate each cluster
		// Pairs [(X,Y), counts]
		list< pair < pair<int, int>, int > >::iterator listItr;

		Log << MSG::DEBUG << "Cluster [" << (*blobsItr).bP.clusterSize
				<< "] " <<  " at : " << (*blobsItr).bP.geoCenter_x << " , "
				<< (*blobsItr).bP.geoCenter_y << " | tot : ";

		for ( listItr = (*blobsItr).blobContent.begin() ; listItr != (*blobsItr).blobContent.end() ; listItr++ ) {

			// pixel corrdinates
			pix = (*listItr).first;
			int tot = GetMatrixElement(pix);
			clusterTOT += tot;

			Log << MSG::DEBUG << tot << "  ";

			if( CalibIsOK() ) {

				cE = CalculateAndGetCalibEnergy(pix, tot);
				//Log << MSG::DEBUG << tot << " cE--> " << cE << " ";

				if (cE < 0.) { // bad pixel or bad calib
					badCalibInCluster = true;
				}

			}

			// Store the calibrated energy for this pixel
			SetCalibEnergy(pix, cE);
			clusterEnergy += cE;

		}
		Log << MSG::DEBUG << endreq;

		// Set Cluster Energy
		(*blobsItr).SetClusterEnergy( clusterEnergy );

		if(SingleHit(blobsItr))
		{
			// change type now that it has been identified
			m_aB->SetBlobType(cntr, _SINGLE_HIT);

			m_basicSpeciesCntr.nSingleHits++;
		}
		else if(DoubleHit(blobsItr))
		{
			m_aB->SetBlobType(cntr, _DOUBLE_HIT);

			m_basicSpeciesCntr.nDoubleHits++;
		}
		else if(TripleHit(blobsItr))
		{
			m_aB->SetBlobType(cntr, _TRIPLE_HIT);

			m_basicSpeciesCntr.nTripleHits++;
		}
		else if(QuadHit(blobsItr))
		{
			m_aB->SetBlobType(cntr, _QUAD_HIT);

			m_basicSpeciesCntr.nQuadHits++;
		}
		else if(LongGamma(blobsItr))
		{
			m_aB->SetBlobType(cntr, _LONG_GAMMA);

			m_basicSpeciesCntr.nLongGamma++;
		}
		else if(Mip(blobsItr) && m_mipON)
		{
			// touch type
			m_aB->SetBlobType(cntr, _MIP);

			m_basicSpeciesCntr.nMip++;
			m_highMips = new Highlighter((*blobsItr).bP.geoCenter_x,
					(*blobsItr).bP.geoCenter_y,
					"arrow", this);
			m_highMips->SetLineColor(kYellow);
			m_highMips->SetLineWidth(1);

			FillValuesForDisplay(m_highMips, cntr, *blobsItr);
			PullToStoreGateAccess(m_highMips, MPXDefs::DO_NOT_SERIALIZE_ME);

		}
		else if(HeavyTrack(blobsItr) && m_HeavyTrackON)
		{
			// touch type
			m_aB->SetBlobType(cntr, _HEAVY_TRACK);

			m_basicSpeciesCntr.nHeavyTracks++;
			m_highHeavyTracks = new Highlighter((*blobsItr).bP.geoCenter_x,
					(*blobsItr).bP.geoCenter_y,
					"circle", this);
			m_highHeavyTracks->SetLineColor(kYellow);
			m_highHeavyTracks->SetLineWidth(1);

			FillValuesForDisplay(m_highHeavyTracks, cntr, *blobsItr);
			PullToStoreGateAccess(m_highHeavyTracks, MPXDefs::DO_NOT_SERIALIZE_ME);

			// get linear reg line for this blob
			//DrawLines((*blobsItr));

		}
		else if(HeavyBlob(blobsItr) && m_HeavyBlobON)
		{
			// touch type
			m_aB->SetBlobType(cntr, _HEAVY_BLOB);

			m_basicSpeciesCntr.nHeavyBlobs++;
			m_highHeavyBlobs = new Highlighter((*blobsItr).bP.geoCenter_x,
					(*blobsItr).bP.geoCenter_y,
					"circle", this);
			m_highHeavyBlobs->SetLineWidth(1);

			FillValuesForDisplay(m_highHeavyBlobs, cntr, *blobsItr);
			PullToStoreGateAccess(m_highHeavyBlobs, MPXDefs::DO_NOT_SERIALIZE_ME);

			// FIXME
			// this cluster size must be removed from here ... misleading variable
			// Heavy blob properties // TODO ... look at this branches in the Init member
			m_clusterSize.push_back( (*blobsItr).bP.nPixels );
			m_clusterSizeWidth.push_back(
					( (*blobsItr).bP.width_x + (*blobsItr).bP.width_y) / 2.
			);
			m_clusterSizeEllipse.push_back( (*blobsItr).bP.ellipseArea );
			m_circleArea.push_back( (*blobsItr).bP.circleArea );
			m_counts.push_back( (*blobsItr).bP.totalCharge );

			// get linear reg line for this blob
			//DrawLines((*blobsItr));

		}
		else if( Curly(blobsItr) )
		{
			// touch type
			m_aB->SetBlobType(cntr, _CURLY);

			m_basicSpeciesCntr.nCurly++;
			m_highCurly = new Highlighter((*blobsItr).bP.geoCenter_x,
					(*blobsItr).bP.geoCenter_y,
					"arrow", this);
			m_highCurly->SetLineWidth(1);
			FillValuesForDisplay(m_highCurly, cntr, *blobsItr);
			PullToStoreGateAccess(m_highCurly, MPXDefs::DO_NOT_SERIALIZE_ME);

		}
		else
		{
			// touch type
			m_aB->SetBlobType(cntr, _NOT_A_BASIC_TYPE);
			m_highUnknown = new Highlighter((*blobsItr).bP.geoCenter_x,
					(*blobsItr).bP.geoCenter_y, 
					"arrow", this);
			m_highUnknown->SetLineWidth(1);
			m_highUnknown->SetLineColor(kMagenta);

			FillValuesForDisplay(m_highUnknown, cntr, *blobsItr);
			PullToStoreGateAccess(m_highUnknown, MPXDefs::DO_NOT_SERIALIZE_ME);

		}
		cntr++;
	}

	m_basicSpeciesCntr.nAll = cntr;

}

bool PRBasicSpecies::SingleHit(vector<blob>::iterator oneBlob){

	if((*oneBlob).bP.nPixels == 1)
		return true;

	return false;
}

bool PRBasicSpecies::DoubleHit(vector<blob>::iterator oneBlob){

	if(
			((*oneBlob).bP.nPixels == 2 && (*oneBlob).bP.width_x == 2 && (*oneBlob).bP.width_y == 2)
			||
			((*oneBlob).bP.nPixels == 2 && (*oneBlob).bP.width_x == 2 && (*oneBlob).bP.width_y == 1)
			||
			((*oneBlob).bP.nPixels == 2 && (*oneBlob).bP.width_x == 1 && (*oneBlob).bP.width_y == 2)
	)
		return true;

	return false;
}

bool PRBasicSpecies::TripleHit(vector<blob>::iterator oneBlob){

	if((*oneBlob).bP.nPixels == 3 && (*oneBlob).bP.width_x == 2 && (*oneBlob).bP.width_y == 2)
		return true;

	return false;
}

bool PRBasicSpecies::QuadHit(vector<blob>::iterator oneBlob){

	if((*oneBlob).bP.nPixels == 4 && (*oneBlob).bP.width_x == 2 && (*oneBlob).bP.width_y == 2)
		return true;

	return false;
}

bool PRBasicSpecies::LongGamma(vector<blob>::iterator oneBlob){

	if((*oneBlob).bP.nPixels >= 3
			&& (*oneBlob).bP.width_x * (*oneBlob).bP.width_y <= m_longGammaMax)
		return true;

	return false;
}

bool PRBasicSpecies::Mip(vector<blob>::iterator oneBlob){


	// Calculate the fraction of points at minimum distance
	// fit and cut information comes with the blob

	Float_t cut = (*oneBlob).bP.fitCut;
	Float_t slope = (*oneBlob).bP.fitSlope;

	Float_t cutPrim = 0.;
	Float_t slopePrim = -1./slope;

	Float_t crossx = 0.;
	Float_t crossy = 0.;
	Float_t pixx = 0.;
	Float_t pixy = 0.;
	Float_t fractionOfPixels = 0.;
	Float_t distancePixelLine = 0.;

	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator listItr;
	Int_t cntr = 0;

	for(listItr = (*oneBlob).blobContent.begin() ; listItr != (*oneBlob).blobContent.end() ; listItr++)
	{
		pixx = (Float_t)(*listItr).first.first;
		pixy = (Float_t)(*listItr).first.second;

		cutPrim  = pixy;
		cutPrim -= slopePrim * pixx;

		crossx = (cut - cutPrim)/(slopePrim - slope);
		crossy = crossx*slope + cut;

		// calculate the perpendicular distance between a point and the fit line
		distancePixelLine = TMath::Sqrt( (crossy - pixy)*(crossy - pixy) + (crossx - pixx)*(crossx - pixx) );

		//Log << MSG::INFO << distancePixelLine << "-->" << " pixx , pixy " << pixx << " , " << pixy << endreq;

		if(distancePixelLine < m_mindistanceToLine)
			fractionOfPixels += 1.;

		cntr++;
	}

	// get the fraction
	fractionOfPixels /= (Float_t)cntr;


	if(//(*oneBlob).bP.chisquare_OverDof < m_mipChisquareDof &&
			(*oneBlob).bP.nPixels > m_mipNPixels &&
			fractionOfPixels > m_fractionOfPixelsAtMindistance)
		return true;

	return false;
}

//bool PRBasicSpecies::LongGamma(blob oneBlob){


//  return false;
//}

bool PRBasicSpecies::HeavyBlob(vector<blob>::iterator oneBlob){

	if(
			(*oneBlob).bP.nInnerPixels >= m_nInnerPixelsCut
			&&
			(*oneBlob).bP.circleArea/(*oneBlob).bP.ellipseArea < m_circleToEllipseMax
	)
		return true;

	return false;
}

bool PRBasicSpecies::HeavyTrack(vector<blob>::iterator oneBlob){

	if(
			(*oneBlob).bP.nInnerPixels >= m_nInnerPixelsCut // same as heavy blob plus ...
			&&
			(*oneBlob).bP.minToGeoCenter < 1 // no blank spaces in between
			&&
			(*oneBlob).bP.circleArea/(*oneBlob).bP.ellipseArea > m_circleToEllipseMin
			&&
			(*oneBlob).bP.circleArea/(*oneBlob).bP.ellipseArea < m_circleToEllipseMax
	)
		return true;

	return false;
}

bool PRBasicSpecies::Curly(vector<blob>::iterator oneBlob){

	if( (*oneBlob).bP.nPixels < m_maxCurlyPixels ) // shorter than a mip
		return true;

	return false;
}



void PRBasicSpecies::Finalize(){


}

void PRBasicSpecies::ClearNtupleVars(){

	m_basicSpeciesCntr.nSingleHits = 0;
	m_basicSpeciesCntr.nDoubleHits = 0;
	m_basicSpeciesCntr.nTripleHits = 0;
	m_basicSpeciesCntr.nQuadHits = 0;
	m_basicSpeciesCntr.nLongGamma = 0;
	m_basicSpeciesCntr.nHeavyBlobs = 0;
	m_basicSpeciesCntr.nHeavyTracks = 0;
	m_basicSpeciesCntr.nMip = 0;
	m_basicSpeciesCntr.nCurly = 0;
	m_basicSpeciesCntr.nNotBasicType = 0;
	m_basicSpeciesCntr.nAll = 0;

	m_clusterSize.clear();
	m_clusterSizeWidth.clear();
	m_clusterSizeEllipse.clear();
	m_circleArea.clear();
	m_counts.clear();

	m_overflowPixelsCntr = 0;

	m_basicGeneral.nHits = 0;
	m_basicGeneral.nCounts = 0;

}

void PRBasicSpecies::FillValuesForDisplay(Highlighter * hl, Int_t cntr, blob bl){

	hl->UploadDisplayValue("Number of Pixels     ", bl.bP.nPixels);
	hl->UploadDisplayValue("Inner Pixels         ", bl.bP.nInnerPixels);
	hl->UploadDisplayValue("Total charge         ", bl.bP.totalCharge);
	hl->UploadDisplayValue("Width x              ", bl.bP.width_x);
	hl->UploadDisplayValue("Width y              ", bl.bP.width_y);
	hl->UploadDisplayValue("Geo center x         ", bl.bP.geoCenter_x);
	hl->UploadDisplayValue("Geo center y         ", bl.bP.geoCenter_y);
	hl->UploadDisplayValue("Weighted center x    ", bl.bP.weightedCenter_x);
	hl->UploadDisplayValue("Weighted center y    ", bl.bP.weightedCenter_y);
	hl->UploadDisplayValue("Box area             ", bl.bP.boxArea);
	hl->UploadDisplayValue("Circle area          ", bl.bP.circleArea);
	hl->UploadDisplayValue("Ellipse area         ", bl.bP.ellipseArea);
	hl->UploadDisplayValue("  circle/elipse      ", bl.bP.circleArea/bl.bP.ellipseArea);
	hl->UploadDisplayValue("Elipse A             ", bl.bP.ellipseA);
	hl->UploadDisplayValue("Elipse B             ", bl.bP.ellipseB);
	hl->UploadDisplayValue("Rotation angle [deg] ", bl.bP.rotAngle * 180 / TMath::Pi());
	hl->UploadDisplayValue("Chisquare/dof        ", bl.bP.chisquare_OverDof);
	hl->UploadDisplayValue("Fit slope            ", bl.bP.fitSlope);
	hl->UploadDisplayValue("Fit cut              ", bl.bP.fitCut);
	hl->UploadDisplayValue("Balance to minimum   ", bl.bP.balanceToMin);
	hl->UploadDisplayValue("Balance to maximum   ", bl.bP.balanceToMin);
	hl->UploadDisplayValue("Min dist to center   ", bl.bP.minToGeoCenter);
	hl->UploadDisplayValue("Max dist to center   ", bl.bP.maxToGeoCenter);


	// values updated in this algorithm ... in this case I need the 'cntr'
	// and I need the pointer to the StoreGate again
	hl->UploadDisplayValue("type of blob         ", (Int_t)m_aB->GetBlobType(cntr));

	// upload the level1
	hl->UploadDisplayValue("lvl1                 ", bl.bP.lvl1);

	//hl->UploadDisplayValue("type of blob (string)", bl.GetTypeAsString());

}

Bool_t PRBasicSpecies::FrameContainsOverFlowPixels(){

	Int_t xDim = GetMatrixXdim();
	Int_t yDim = GetMatrixYdim();
	Bool_t foundBadPixel = false;

	for(Int_t colItr = 0; colItr < xDim ; colItr++)
	{
		for(Int_t rowItr = 0; rowItr < yDim ; rowItr++)
		{
			if(GetMatrixElement(colItr, rowItr) == TOT_OVERFLOW_COUNTS)
			{
				Log << MSG::DEBUG << "pixel bruyant at " << colItr << " , "  << rowItr << endreq ;
				m_overflowPixelsCntr++;
				foundBadPixel = true;
			}
		}
	}

	if(foundBadPixel)
	{
		return true;
	}

	return false;
}

void PRBasicSpecies::DrawLines(blob oneBlob){

	// f1
	// The linear regresion from BlobsFinder
	Double_t a1 = oneBlob.bP.fitSlope;
	Double_t b1 = oneBlob.bP.fitCut; // not needed
	Double_t geox = oneBlob.bP.geoCenter_x;

	// draw stuff
	// f1
	Double_t c1 = TMath::Sqrt(oneBlob.bP.width_x*oneBlob.bP.width_x
			+
			oneBlob.bP.width_y*oneBlob.bP.width_y); // length of the line
	Double_t alp = TMath::ATan(a1);
	Double_t y_delta = c1*TMath::Sin(alp);
	Double_t x_delta = TMath::Sqrt(c1*c1 - y_delta*y_delta);
	Highlighter * high_f1 = new Highlighter(geox - x_delta,
			a1*(geox - x_delta) + b1,
			geox + x_delta,
			a1*(geox + x_delta) + b1,
			"line", this);
	PullToStoreGateAccess(high_f1, MPXDefs::DO_NOT_SERIALIZE_ME);
	high_f1->SetLineWidth(1);
	high_f1->SetLineStyle(1);
}

#endif
