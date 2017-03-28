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

#ifndef FEI3Mips_cpp
#define FEI3Mips_cpp

#include "FEI3Mips.h"
#include "MAFTools.h"
#include "MPXAlgo/Highlighter.h"

using namespace MSG;

ClassImp(FEI3Mips)

FEI3Mips::FEI3Mips() {

	m_minNPixels = 5 ;
	m_minMipLength = 3.0; // [mm]
	m_nDivisions = 4 ;
	m_pixSizeX = 0.055; // Timepix // 0.40; // FEI3 [mm]
	m_pixSizeY = 0.055; // Timepix // 0.05; // FEI3 [mm]
	m_FEI3_PIXELSIZE_YX_FRACTION = m_pixSizeY/m_pixSizeX;
	m_totalMIPsLength = 0.0;
	m_dRayBalSearch = 1.4;
	m_distanceCutHot = 3;
	m_distancePerpHot = 3;
	m_guardDistanceX = 1;
	m_guardDistanceY = 5;

	m_totalMIPsLength = 0.0;
	m_nDeltaRay = 0;
	m_nDeltaRaySoft = 0;

	m_maxYWidth = 3;
	m_newMeshDiv = 1;

	m_extraCanvas = 0x0;

}

void FEI3Mips::Init(){

	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("chargeWeights", &m_chargeWeights);
	getMyTree()->Branch("lvl1Weights", &m_lvl1Weights);
	getMyTree()->Branch("chargePerSection", &m_chargePerSection);
	getMyTree()->Branch("deltaRayTOT", &m_deltaRayTOT);
	getMyTree()->Branch("mipLength", &m_mipLength);
	getMyTree()->Branch("meanChargePerSection", &m_meanChargePerSection);
	getMyTree()->Branch("sectionCenter_x", &m_sectionCenter_x);
	getMyTree()->Branch("sectionCenter_y", &m_sectionCenter_y);
	getMyTree()->Branch("nDeltaRayExplicit", &m_nDeltaRayPerFrame);
	getMyTree()->Branch("nDeltaRaySoft", &m_nDeltaRaySoftPerFrame);

	getMyTree()->Branch("clusterSize", &m_clusterSize, "clusterSize/I");
	getMyTree()->Branch("alphaAngle", &m_alphaAngle, "alphaAngle/D");
	getMyTree()->Branch("frameId", &m_frameId, "frameId/I");

	// A configuration value that can be tuned from the Viewer
	RegisterConfigurationValue(&m_minNPixels, "minNPixels");
	RegisterConfigurationValue(&m_minMipLength, "minMipLength");
	RegisterConfigurationValue(&m_nDivisions, "nDivisions");
	RegisterConfigurationValue(&m_pixSizeX, "pixelSizeX");
	RegisterConfigurationValue(&m_pixSizeY, "pixelSizeY");
	RegisterConfigurationValue(&m_dRayBalSearch, "dRayBalSearch");
	RegisterConfigurationValue(&m_distanceCutHot, "distanceCutHot");
	RegisterConfigurationValue(&m_distancePerpHot, "distancePerpHot");
	RegisterConfigurationValue(&m_guardDistanceX, "guardDistanceX");
	RegisterConfigurationValue(&m_guardDistanceY, "guardDistanceY");
	RegisterConfigurationValue(&m_maxYWidth, "maxYWidth");
	RegisterConfigurationValue(&m_newMeshDiv, "newMeshDiv");

}

void FEI3Mips::Execute(){


	// Clear drawing objetcs
	if(!m_graph2DVector.empty()) m_graph2DVector.clear();
	if(m_extraCanvas) delete m_extraCanvas;

	// A check on parameters
	// The number of divisions can not be smaller than the minimum amount of
	// pixels required in a mip
	if(m_nDivisions > m_minNPixels) {
		Log << MSG::ERROR << "The number of divisions can not be bigger than the minimum amount of" << endreq;
		Log << MSG::ERROR << "pixels required in a mip to process.  Fix the parameters.  Giving up." << endreq;
		exit(1);
	}

	// Ask the store gate if the previous algorithm (BlobsFinder --> reponsible for clustering)
	//  sent any objects to the StoreGate.
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");

	// If so, get the pointer to the last object.  BlobsFinder, as it comes out of the box, sends
	//  a single object containing all clusters.
	m_aB = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject-1);

	// AllBlobsContainer is a box full of Blobs(Clusters). Now we can iterate over all clusters inside.
	vector<blob> blobsVector = m_aB->GetBlobsVector();

	// If nothing to do, return
	if(blobsVector.empty()) return;

	Log << MSG::INFO << "Number of blobs from clustering = " << (Int_t) blobsVector.size() << endreq;
	vector<blob>::iterator blobsItr = blobsVector.begin(); //allBlobs.begin();

	for( ; blobsItr != blobsVector.end() ; blobsItr++)
	{

		// Limit all this to clusters with a minimum size.
		// Note that m_minNPixels can be configured through the Viewer
		//  so you can reprocess and check results online.
		if((*blobsItr).bP.nPixels < m_minNPixels)
			continue;

		// Look for mips
		blobtype bT = (*blobsItr).GetBlobType();

		if(bT == _MIP) { // only Mips

			// info about the mip to store
			m_clusterSize = (*blobsItr).GetBlobProperties().nPixels;

			// Here we obtain two pairs with the pixels at the extreme ends of the Mip
			//vector< pair<Int_t, Int_t> > limitPixels = FindLimitPixels(*blobsItr);
			//i-
			vector< pair<double, double> > limitPixels = FindLimitPixels_MaxMinXY(*blobsItr);

			Log << MSG::INFO << "Limits : " << limitPixels[__HEAD_INDX].first << "," << limitPixels[__HEAD_INDX].second <<
					" --> " << limitPixels[__TAIL_INDX].first << "," << limitPixels[__TAIL_INDX].second << endreq;

			// If the cluster is too wide in Y don't processi t
			if( TMath::Abs( limitPixels[__TAIL_INDX].second - limitPixels[__HEAD_INDX].second ) > m_maxYWidth) {
				continue;
			}

			// Make sure the mip starts and ends inside the pad at a certain minimum
			//  distance from the borders
			if ( MAFTools::BlobAtADistanceFromEdges(*blobsItr, limitPixels, GetWidth(), GetHeight(),
					m_guardDistanceX, m_guardDistanceY) ) {
				continue;
			}

			// create all slots to put the weighted charge and others observables
			for(int ni = 0 ; ni < m_nDivisions ; ni++){
				m_chargeWeights.push_back(0.0);
				m_lvl1Weights.push_back(0.0);
				m_chargePerSection.push_back(vector<Int_t>()); // push back empty vectors
				m_XYPerSection.push_back(vector<pair <Int_t, Int_t> >());
				m_meanChargePerSection.push_back(0.0);
				m_sectionCenter_x.push_back(0.0);
				m_sectionCenter_y.push_back(0.0);
			}

			// If the whole mip is in the same column, don't process it.
			// Problems with reading of too many pixels activated in the
			//  same column in FE-I3.  FE-I4 too ?
			if (limitPixels[__TAIL_INDX].first - limitPixels[__HEAD_INDX].first == 0) {
				ClearThisAlgo();
				continue;
			}

			// OK, if the mip is fine, collect a few stats.
			//  This is distance, as opposed to that calculated for divOnTrack, has to have
			//  real units.
			Double_t totalLength = MAFTools::CalcDistance(limitPixels[__HEAD_INDX].first * m_pixSizeX,
					limitPixels[__HEAD_INDX].second * m_pixSizeY,
					limitPixels[__TAIL_INDX].first * m_pixSizeX,
					limitPixels[__TAIL_INDX].second * m_pixSizeY);
			m_mipLength.push_back( totalLength ); // mm
			m_totalMIPsLength += totalLength;
			Log << MSG::INFO << "Integrated mips length = " << m_totalMIPsLength << " [mm]" << endreq;

			// cut by mip length
			if ( totalLength < m_minMipLength ) {
				ClearThisAlgo();
				continue;
			}

			// Fech the linear regression information for this mip
			Float_t slope = (*blobsItr).GetBlobProperties().fitSlope;
			Float_t cut = (*blobsItr).GetBlobProperties().fitCut;

			// Draw the linear regression and the line which passes through the limits
			MAFTools::DrawLine(this, limitPixels[__HEAD_INDX].first, limitPixels[__HEAD_INDX].second,
					limitPixels[__TAIL_INDX].first, limitPixels[__TAIL_INDX].second, 2, 2, kBlack);
			MAFTools::DrawLine(this, limitPixels[__HEAD_INDX].first, limitPixels[__HEAD_INDX].first * slope + cut,
					limitPixels[__TAIL_INDX].first, limitPixels[__TAIL_INDX].first * slope + cut, 2, 2, kRed);

			// I won't work with the linear regression but with the further most pixels,
			//  otherwise the division is not even.  Change the slope and cut values.

			slope = (Float_t(limitPixels[__TAIL_INDX].second - limitPixels[__HEAD_INDX].second))/
					(Float_t(limitPixels[__TAIL_INDX].first - limitPixels[__HEAD_INDX].first));
			cut = limitPixels[__TAIL_INDX].second - slope*limitPixels[__TAIL_INDX].first;

			// Find the points for divisions as requested by m_nDivisions.
			// Lenght of each section in the track
			Double_t divOnTrack = MAFTools::CalcDistance(limitPixels[__HEAD_INDX], limitPixels[__TAIL_INDX]);
			divOnTrack /= m_nDivisions;

			// Calculate the angle of the track with respect to the X axis
			Float_t alpha = TMath::ATan(slope);
			// The angle going to ntuple has to be corrected to the size of the pixel X/Y = 0.05/0.4
			m_alphaAngle = TMath::ATan(slope * m_FEI3_PIXELSIZE_YX_FRACTION)*180.0/TMath::Pi();
			// The projection value would be
			Float_t divProj = divOnTrack * TMath::Cos(alpha);

			Log << MSG::INFO << "Track Length (arbitrary units) = " << divOnTrack*m_nDivisions << endreq;
			Log << MSG::INFO << "sectionOnTrack = " << divOnTrack << " ,  divProj = " << divProj << endreq;

			// Decide where to start, Always by the left
			Int_t x1i = limitPixels[__HEAD_INDX].first;
			Int_t x2i = limitPixels[__TAIL_INDX].first;
			Int_t startX = x1i;
			if(x2i < x1i) startX = x2i;

			// Build the rest of the point within the boundaries
			for (int i = 1 ; i < m_nDivisions ; i++) {

				Float_t newPX = startX + i*divProj;
				Float_t newPY = newPX*slope + cut;

				Log << MSG::INFO << "division point : " << newPX << " , " << newPY << endreq;

				//limitPixels.push_back( make_pair( TMath::FloorNint(newPX), TMath::FloorNint(newPY) ) );
				//i-
				limitPixels.push_back( make_pair( newPX, newPY ) );
			}

			// Now I need lines perpendicular to the line chosen (given by the extremes pixels)
			// Stored as a vector of a and b (as in ax + b) of the parallel lines
			vector< pair<Float_t, Float_t > > linePars;

			// Flag which indicates if we are dealing with a track with slope = 0.  Special case.
			bool parallelLine = false;
			for (int i = 0 ; i < m_nDivisions+1 ; i++) { // there are m_nDivisions+1 lines !

				// WARNING ! If slope == 0. a flag indicating it is raised
				linePars.push_back(  GetParallelLineAt(slope, cut, limitPixels[i], &parallelLine)  );

			}

			//////////////////////////////////////////////////////////////////////////////////
			// Ready to calculate fraction of charge in each division
			//////////////////////////////////////////////////////////////////////////////////
			// linePars contains all the parallel lines of the grid <slope, cut>
			//   the first two lines are the extremes and the rest are all those in between.
			// I need to loop over the components of the mip and decide in which section
			// they are placed.
			vector< pair<Float_t, Float_t > >::iterator separatorItr = linePars.begin();

			set< pair<Int_t, Int_t> > contentSet = (*blobsItr).GetContentSet();
			set< pair<Int_t, Int_t> >::iterator contItr = contentSet.begin();

			// Iterate over pixels in mip
			for( ; contItr != contentSet.end() ; contItr++) {

				// I will build a bitword
				UInt_t locationBitWord = 0x0;
				Int_t x = (*contItr).first;
				Int_t y = (*contItr).second;

				// Loop in the inverse order to get the bitwork little endian
				Log << MSG::DEBUG << "Section size = " << divOnTrack << " (arbitrary units !)" << endreq;

				for( separatorItr = linePars.end()-1 ; separatorItr != linePars.begin()-1 ; separatorItr-- ) {

					// deal with the slope=0 case
					Double_t perpDistance = 0.;
					if(parallelLine) {
						//Log << MSG::DEBUG << x << " --> " << (*separatorItr).first << endreq;
						perpDistance = TMath::Abs( x - (*separatorItr).first);
					} else {
						// calculate the perpendicular distance to each line
						perpDistance = MAFTools::CalcPerpDistanceToLine( (*separatorItr).first,
								(*separatorItr).second, *contItr );
					}

					Log << MSG::DEBUG << "pixel : " << x << " , " << y << " | perpDistance = " << perpDistance << " ";

					// put a 1 in the chain only if the distance is greater than the division
					if (perpDistance > divOnTrack) { // --> 1
						locationBitWord ^= 0x1;
						Log << MSG::DEBUG << " 1 ";
					} else {
						//locationBitWord ^= 0x0;
						Log << MSG::DEBUG << " 0 ";
					}
					// switch except last time
					if(separatorItr != linePars.begin()) locationBitWord = locationBitWord << 1;
					Log << MSG::DEBUG << endreq;

				}

				Log << MSG::DEBUG << "bitWord --> " << MAFTools::DumpWordInBinary(locationBitWord) << endreq;

				// Find first if the hit is in a border, which
				//  can be done using the first three bits.
				// Get the first three bits
				UInt_t firstThreeBitsMask = 0x7;
				UInt_t firstThreeBits = locationBitWord & firstThreeBitsMask;
				Log << MSG::DEBUG << "firstThree --> " << MAFTools::DumpWordInBinary(firstThreeBits) << endreq;

				if( firstThreeBits == __LEFT_CORNER_MASK_1
						|| firstThreeBits == __LEFT_CORNER_MASK_2){

					Log << MSG::DEBUG << "... left corner ... " << endreq;
					m_chargeWeights[__HEAD_INDX] += GetMatrixElement(*contItr);
					m_lvl1Weights[__HEAD_INDX] += GetLVL1(*contItr);

					// store charge for landau distribution per section
					m_chargePerSection[__HEAD_INDX].push_back( GetMatrixElement(*contItr) );
					m_XYPerSection[__HEAD_INDX].push_back( make_pair(x, y) );

					// mean charge per section
					m_meanChargePerSection[__HEAD_INDX] += GetMatrixElement(*contItr);

					// mean position
					m_sectionCenter_x[__HEAD_INDX] += x;
					m_sectionCenter_y[__HEAD_INDX] += y;


				}else if( firstThreeBits == __RIGHT_CORNER_MASK_1
						|| firstThreeBits == __RIGHT_CORNER_MASK_2){

					Log << MSG::DEBUG << "... right corner ... " << endreq;
					m_chargeWeights[__TAIL_INDX] += GetMatrixElement(*contItr);
					m_lvl1Weights[__TAIL_INDX] += GetLVL1(*contItr);

					// store charge for landau distribution per section
					m_chargePerSection[__TAIL_INDX].push_back( GetMatrixElement(*contItr) );
					m_XYPerSection[__TAIL_INDX].push_back( make_pair(x, y) );

					// mean charge per section
					m_meanChargePerSection[__TAIL_INDX] += GetMatrixElement(*contItr);

					// mean position
					m_sectionCenter_x[__TAIL_INDX] += x;
					m_sectionCenter_y[__TAIL_INDX] += y;

				}else{

					// The rest of the words build up the following way.
					// For instace.
					// second slice :  ...  1110011
					// third slice  :  ... 11100111
					// etc.

					Int_t nB = NumberOfBitsOnBeforeDoubleZero(locationBitWord);

					// Handle special case when the pixel is equidistant to two lines.
					// If for example we end up with a word  001111  for m_nDivisions = 4
					//  in fact this pixel belongs to the 3rd slice.
					if(nB == m_nDivisions){
						nB--;
					}

					Log << MSG::DEBUG << "... column " << nB << " ... " << GetFrameId() << endreq;
					m_chargeWeights[nB] += GetMatrixElement(*contItr);
					m_lvl1Weights[nB] += GetLVL1(*contItr);

					// Store charge for landau distribution per section
					m_chargePerSection[nB].push_back( GetMatrixElement(*contItr) );
					m_XYPerSection[nB].push_back( make_pair(x, y) );

					// charge per section
					m_meanChargePerSection[nB] += GetMatrixElement(*contItr);

					// mean position
					m_sectionCenter_x[nB] += x;
					m_sectionCenter_y[nB] += y;

				}

			}

			// FIXME !!!
			// finish up the mean charge per section and position
			for (int ix = 0 ; ix < (int)m_chargePerSection.size() ; ix++) {
				int nEntriesSection = m_chargePerSection[ix].size();
				m_meanChargePerSection[ix] /= nEntriesSection;
				m_sectionCenter_x[ix] /= nEntriesSection;
				m_sectionCenter_y[ix] /= nEntriesSection;
			}

			// normalize by the total charge
			vector<Float_t>::iterator itr = m_chargeWeights.begin();
			int cntr = 0 ;
			for( ; itr != m_chargeWeights.end() ; itr++) {
				(*itr) /= (*blobsItr).GetBlobProperties().totalCharge;
				(*itr) *= 100.; // percentage
				Log << MSG::INFO << "Charge [" << cntr << "] : " << (*itr) << " %" << endreq;
				cntr++;
			}

			// Try a simple guess on the number of delta rays from unbalance
			Int_t nDeltaR = 0;
			Int_t nDeltaSoft = 0;

			GuessNumberOfDeltaRaysFromUnbalance(nDeltaR, nDeltaSoft);
			m_nDeltaRay += nDeltaR;
			m_nDeltaRaySoft += nDeltaSoft;

			m_nDeltaRayPerFrame += nDeltaR; // this one is cleared on a frame-per-frame basis
			m_nDeltaRaySoftPerFrame += nDeltaSoft; // this one is cleared on a frame-per-frame basis

			// Look at the hard delta rays trying to extract an energy profile
			if (nDeltaR > 0) {
				m_deltaRayTOT = ExtractDeltaRayEnergy(*blobsItr);
			}

			// WARNING, filling once per mip !!!
			m_frameId = GetFrameId();
			getMyTree()->Fill();

			m_chargeWeights.clear();
			m_lvl1Weights.clear();
			m_chargePerSection.clear();
			m_XYPerSection.clear();
			m_mipLength.clear();
			m_meanChargePerSection.clear();
			m_sectionCenter_x.clear();
			m_sectionCenter_y.clear();

			m_alphaAngle = 0.;
			m_clusterSize = 0;
			m_nDeltaRayPerFrame = 0;
			m_nDeltaRaySoftPerFrame = 0;

		} else {
			;
		}

	}

	// Draw
	m_extraCanvas = DrawInSeparateWindow(m_graph2DVector);

	// Fill the output tree of this algorithm
	//getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_chargeWeights.clear();
	m_lvl1Weights.clear();
	m_chargePerSection.clear();
	m_XYPerSection.clear();
	m_deltaRayTOT.clear();
	m_mipLength.clear();

}

void FEI3Mips::ClearThisAlgo(){

	m_chargeWeights.clear();
	m_lvl1Weights.clear();
	m_chargePerSection.clear();
	m_XYPerSection.clear();
	m_mipLength.clear();
	m_meanChargePerSection.clear();
	m_sectionCenter_x.clear();
	m_sectionCenter_y.clear();

}

void FEI3Mips::Finalize() {

	Log << MSG::INFO << "Total mips length = " << m_totalMIPsLength << " [mm]" << endreq;
	Log << MSG::INFO << "Number of explicit delta rays = " << m_nDeltaRay << endreq;
	Log << MSG::INFO << "Number of soft delta rays = " << m_nDeltaRaySoft << endreq;

	Log << MSG::INFO << "Number of explicit delta rays per centimeter = " << m_nDeltaRay/(m_totalMIPsLength*0.1) << endreq;
	Log << MSG::INFO << "Number of soft delta rays per centimeter = " << m_nDeltaRaySoft/(m_totalMIPsLength*0.1) << endreq;
	Log << MSG::INFO << "Total Number delta rays per centimeter = " << (m_nDeltaRay+m_nDeltaRaySoft)/(m_totalMIPsLength*0.1) << endreq;

}

vector<Float_t> FEI3Mips::ExtractDeltaRayEnergy (blob bl) {

	// At this point the information in m_chargeWeights is already
	// expressed in percentage per segment
	vector<Float_t>::iterator itr = m_chargeWeights.begin();
	Int_t indxMaxSeg = 0;
	Float_t maxSeg = 0.;

	Int_t cntr = 0;
	// Find the hotest segment
	for ( ; itr != m_chargeWeights.end() ; itr++ ) {

		if ( (*itr) > maxSeg ) {
			maxSeg = (*itr);
			indxMaxSeg = cntr;
		}
		cntr++;
	}

	//  indxMaxSeg is the index of the hotest segment
	Log << MSG::INFO << "*** Delta ray Search ***" << endreq;
	Log << MSG::INFO << "Hotest segment id : [" << indxMaxSeg << "] = " << maxSeg << "%"<< endreq;

	// Now in this segment find the hotest point ( look at indxMaxSeg ! )
	vector<Int_t>::iterator cItr = m_chargePerSection[indxMaxSeg].begin();
	vector< pair<Int_t, Int_t> >::iterator xyItr = m_XYPerSection[indxMaxSeg].begin();
	Int_t highC = 0;
	Int_t highIndx = 0;
	cntr = 0;
	for ( ; cItr !=  m_chargePerSection[indxMaxSeg].end() ; ) {

		if( *cItr > highC ) {
			highC = *cItr;
			highIndx = cntr;
		}

		cntr++;
		cItr++;
		xyItr++;
	}

	Log << MSG::INFO << "Hot spot at <" << m_XYPerSection[indxMaxSeg][highIndx].first << ", "
			<< m_XYPerSection[indxMaxSeg][highIndx].second << " with TOT = "
			<< m_chargePerSection[indxMaxSeg][highIndx] << endreq;

	// Once having the hotest pixel I can search for all the charge belonging to the
	// delta ray which has been emitted at m_XYPerSection[indxMaxSeg][highIndx].
	//
	// The approach is the following:
	// 1) Pixels very close to the hotest pixel probably have a bit of a higher TOT that
	//    the mean energy in the rest of the mip track.  That extra energy should be added
	//    to the delta ray spectrum.
	// 2) If the delta ray stayed in the plane of the sensor, there will probably be a few
	//    pixels ON close to the mip track with the shape of a low energy electron (curly).
	//    This energy ought to be added too.

	vector<Float_t> deltaRaysTOT;

	// The mean TOT in all the track is
	// TODO ... in this mean I should extract the hot pixels
	Float_t meanTOT = Float_t( bl.GetBlobProperties().totalCharge ) /
			Float_t( bl.GetBlobProperties().nPixels );

	Log << MSG::INFO << "Mean TOT of track = " << meanTOT << endreq;

	// calculate perp distance in each point and also the distance to the hot point
	Float_t slope = bl.GetBlobProperties().fitSlope;
	Float_t cut = bl.GetBlobProperties().fitCut;
	// Prepare perpendicular line passing by hot pixel
	Float_t slopePrim = -1/slope;
	Float_t cutPrim = m_XYPerSection[indxMaxSeg][highIndx].second;
	cutPrim -= slopePrim * m_XYPerSection[indxMaxSeg][highIndx].first;

	list< pair < pair<Int_t, Int_t>, Int_t > > contents = bl.GetBlobContent();
	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator bcItr = contents.begin();

	Double_t distP = 0., distHot = 0., distHotPerp = 0.;
	Float_t deltaRayTotalTOT = 0.;

	// draw the perp line passing by hot pixel
	//Float_t x1p = m_XYPerSection[indxMaxSeg][highIndx].first;
	//Float_t y1p = m_XYPerSection[indxMaxSeg][highIndx].second;
	//Float_t x2p = x1p*1.5;
	//Float_t y2p = slopePrim * x2p + cutPrim;
	//MAFTools::DrawLine(this, x1p, y1p, x2p, y2p, 2, 1, kBlack);

	// Arrows to graphically identify delta rays
	bool hotSpotDrawn = false;
	for ( ; bcItr != contents.end() ; bcItr++ ) {

		// perpendicular distance to Linear Regression
		distP = MAFTools::CalcPerpDistanceToLine(slope, cut, (*bcItr).first);
		// distance to hot pixel (vertex)
		distHot = MAFTools::CalcDistance(m_XYPerSection[indxMaxSeg][highIndx], (*bcItr).first);
		// perpendicular distance to perp to Linear Regression passing by hot pixel
		distHotPerp = MAFTools::CalcPerpDistanceToLine(slopePrim, cutPrim, (*bcItr).first);

		// Select all the points surrounding the hotest spot
		if( distP < 1.0 && ( distHot < m_distanceCutHot || distHotPerp < m_distancePerpHot ) ) { // close to vertex

			if(!hotSpotDrawn) {

				Highlighter * fp = new Highlighter((*bcItr).first.first,
						(*bcItr).first.second,
						"circle", this);
				fp->SetLineColor(kRed);
				fp->SetLineWidth(1);
				PullToStoreGateAccess(fp, MPXDefs::DO_NOT_SERIALIZE_ME);
				hotSpotDrawn = true;

				// Before anything let's have a finner grid
				TGraph2D * g = MAFTools::ConvertClusterToGraph2D(bl, m_newMeshDiv, GetFrameId(), "E");
				m_graph2DVector.push_back( g );

			}

			Log << MSG::DEBUG << "perp dist = " << distP << ", distance to Hot = " << distHot << endreq;
			// Finally store the excess energy if it is over the mean
			Float_t TOTDiff = GetMatrixElement( (*bcItr).first ) - meanTOT;
			if( TOTDiff > 0 ) {
				// FIXME !, store the whole energy ! --> why ?
				deltaRayTotalTOT += TOTDiff;
			}
		}

		// At this point I can search for the rest of the delta ray if it stayis in the
		// plane of the sensor


	}

	deltaRaysTOT.push_back(deltaRayTotalTOT);

	return deltaRaysTOT;
}

Int_t FEI3Mips::GuessNumberOfDeltaRaysFromUnbalance(Int_t & nDR, Int_t & nDRSoft){

	vector<Float_t>::iterator itr = m_chargeWeights.begin();
	vector<Float_t>::iterator itrC = m_chargeWeights.begin();
	vector<Int_t> nNoticeableUnbalances(m_nDivisions);  // initialize at 0

	Float_t unbalance = 1.;
	Int_t cntr = 0, cntrE = 0;
	Int_t nExplicit = 0, nSoft = 0;

	for ( ; itr != m_chargeWeights.end() ; itr++ ) {

		cntrE = 0;

		for ( itrC = m_chargeWeights.begin() ; itrC != m_chargeWeights.end() ; itrC++ ) {

			// unbalance between the two fragments
			unbalance = (*itr)/(*itrC);
			if ( unbalance > m_dRayBalSearch ) {
				Log << MSG::INFO << "section " << cntr << " much hotter than " << cntrE << endreq;
				nNoticeableUnbalances[cntr]++;
			}

			cntrE++;
		}

		cntr++;
	}

	// based on nNoticeableUnbalances and m_nDivisions determine if there is a
	// delta ray or not

	vector<Int_t>::iterator itrU = nNoticeableUnbalances.begin();  // initialize at 0

	for( ; itrU != nNoticeableUnbalances.end() ; itrU++){

		if ( (*itrU) == m_nDivisions - 1 ) { // an explicit delta ray.  Unbalance to all other segments
			nExplicit++;
		} else if ( (*itrU) == m_nDivisions - 2 ) {
			nSoft++;
		}

	}

	Log << MSG::INFO << "Explicit : " << nExplicit << " , Soft : " << nSoft << endreq;

	nDR = nExplicit;
	nDRSoft = nSoft;

	return nExplicit;
}

vector<pair<double, double> > FEI3Mips::FindLimitPixels(blob theBlob) {

	// result vector
	vector<pair<double, double> > limits;
	limits.push_back( make_pair(0, 0) );
	limits.push_back( make_pair(0, 0) );

	set< pair<Int_t, Int_t> > contents = theBlob.GetContentSet();

	set< pair<Int_t, Int_t> >::iterator contItrA = contents.begin();
	set< pair<Int_t, Int_t> >::iterator contItrB = contents.begin();

	double max = 0., dist = 0.;

	for( ; contItrA != contents.end() ; contItrA++) {

		for(contItrB = contItrA ; contItrB != contents.end() ; contItrB++) {

			dist = MAFTools::CalcDistance( *contItrA, *contItrB );

			if( dist > max ) {

				max = dist;

				limits[0] = *contItrA;
				limits[1] = *contItrB;

				// The limit in the extreme right needs an adjustement of 1 pixel
				limits[1].first += 1; // 1 pixel

			}

		}

	}

	/*
	Highlighter * limitA = new Highlighter(limits[0].first,
			limits[0].second,
			"arrow", this);
	limitA->SetLineColor(kRed);
	limitA->SetLineWidth(1);
	PullToStoreGateAccess(limitA, MPXDefs::DO_NOT_SERIALIZE_ME);

	Highlighter * limitB = new Highlighter(limits[1].first,
			limits[1].second,
			"arrow", this);
	limitB->SetLineColor(kRed);
	limitB->SetLineWidth(1);
	PullToStoreGateAccess(limitB, MPXDefs::DO_NOT_SERIALIZE_ME);
	 */

	return limits;
}

vector<pair<double, double> > FEI3Mips::FindLimitPixels_MaxMinXY(blob theBlob) {

	// result vector
	vector<pair<double, double> > limits;
	limits.push_back( make_pair( 0 , 0) );
	limits.push_back( make_pair( 0 , 0) );

	set< pair<Int_t, Int_t> > contents = theBlob.GetContentSet();
	set< pair<Int_t, Int_t> >::iterator i = contents.begin();

	int minx = 0x0fffffff, maxx = 0, miny = 0x0fffffff, maxy = 0;

	for( ; i != contents.end() ; i++) {

		if((*i).first  < minx) minx = (*i).first;

		if((*i).second < miny) miny = (*i).second;

		if((*i).first  > maxx) maxx = (*i).first;

		if((*i).second > maxy) maxy = (*i).second;

	}

	limits[0] = make_pair( (double)minx, (double)miny );
	limits[1] = make_pair( (double)maxx, (double)maxy );

	return limits;
}


// FIXME ... shouldn't this be called "Perpendicular" ? just check the logics.
//pair<Float_t, Float_t> FEI3Mips::GetParallelLineAt(Float_t slope, Float_t cut, pair<Int_t, Int_t> point, bool * pL) {
//i-
pair<Float_t, Float_t> FEI3Mips::GetParallelLineAt(Float_t slope, Float_t /*cut*/, pair<Float_t, Float_t> point, bool * pL) {

	if ( slope == 0. ) {
		// In this trivial case, I basically return the pair 'point'.
		// I will only need the x coordinate to calculate the distance later on
		*pL = true;
		return point;

	}

	Float_t cutPrim = 0.;
	Float_t slopePrim = -1./slope;

	// I need to decide on a few points, the extreme and other few
	cutPrim  = point.second;
	cutPrim -= slopePrim * point.first;

	Highlighter * lineD = new Highlighter(0,
			cutPrim,
			18,
			18*slopePrim + cutPrim,
			"line", this);
	lineD->SetLineColor(kRed);
	lineD->SetLineWidth(1);
	PullToStoreGateAccess(lineD, MPXDefs::DO_NOT_SERIALIZE_ME);

	return make_pair(slopePrim, cutPrim);
}

Int_t FEI3Mips::NumberOfBitsOnBeforeDoubleZero(UInt_t word){


	Int_t sizeOfWord = sizeof(word)*8; // number of bits
	UInt_t oneBit = 0x0;
	Int_t nB = 0;

	for (int i = 0 ; i < sizeOfWord ; i++) { // number of bits

		oneBit = word & 0x1;

		if(oneBit == 0x1) { nB++; }
		else { return nB; }

		word = word >> 1;
	}

	return 0;
}

#endif
