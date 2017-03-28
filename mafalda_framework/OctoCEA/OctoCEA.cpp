 /*
 * 	Copyright 2011 John Idarraga
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

#ifndef OctoCEA_cpp
#define OctoCEA_cpp

#include "OctoCEA.h"
#include "MAFTools.h"

#include <map>
#include <vector>

using namespace MSG;

ClassImp(OctoCEA)

OctoCEA::OctoCEA(){
	m_minClusterSize = 10;
	m_joinSlopeAngleTolerance = 5; // deg
	m_joinCutPixelTolerance = 20; // pixels
	m_minNPixels = 5;
	m_minNPixelsClouds = 10;

	m_tracksPtr = 0x0;

	// for output
	nSegments = 0;
}

void OctoCEA::Init(){

	Log << MSG::INFO << "Init function !" << endreq;

	// This value will be overridden by the configuration since it'll be set up
	//  a few lines below as a configuration value

	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("distanceToCenter", &m_distanceToCenter);
	getMyTree()->Branch("AssemblyclusterSize", &clusterSize);
	getMyTree()->Branch("alphaAngle", &m_alphaAngle);
	getMyTree()->Branch("chi2", &m_chi2);
	getMyTree()->Branch("chi2Prob", &m_chi2Prob);
	getMyTree()->Branch("typeOfRecoFrame",&m_typeOfRecoFrame);
	getMyTree()->Branch("segmentsAngles", &m_segmentsAngles);
	getMyTree()->Branch("clusterSize", &m_clusterSize);

	getMyTree()->Branch("assemblyTime", &m_assemblyTime);
	getMyTree()->Branch("frameId", &m_frameId);
	getMyTree()->Branch("occupancy", &m_occupancy, "occupancy/D");
	getMyTree()->Branch("perpDist", &m_perpDist);
	getMyTree()->Branch("perpDistSegments", &m_perpDistSegments);
	getMyTree()->Branch("residual", &m_residual);
	getMyTree()->Branch("residualAligned", &m_residualAligned);


	getMyTree()->Branch("odf",&m_odf);
	getMyTree()->Branch("nSegments", &nSegments, "nSegments/I");


	// A configuration value that can be tuned from the Viewer
	RegisterConfigurationValue(&m_minNPixels, "minNPixels");
	RegisterConfigurationValue(&m_minNPixelsClouds, "minNPixelsClouds");
	RegisterConfigurationValue(&m_minClusterSize, "minClusterSize");
	RegisterConfigurationValue(&m_joinSlopeAngleTolerance, "joinSlopeAngleTolerance");
	RegisterConfigurationValue(&m_joinCutPixelTolerance, "joinCutPixelTolerance");
	RegisterConfigurationValue(&m_cutAngleValue, "cutAncleValue");
	//RegisterConfigurationValue(&m_occupancy, "occupancy");

}

void OctoCEA::Execute(){

	// calculate and store the occupancy
	m_occupancy = (double)GetHitsInPad()/(512. * 1024.);
	m_occupancy *= 100.;

	Log << MSG::INFO << "Occupancy = " << m_occupancy << " %"<< endreq;

	if(m_occupancy > 1) return;

	m_frameId = GetFrameId();
	m_typeOfRecoFrame = __RECO_UNDEFINED;

	Log << MSG::INFO << "Frame Id = " << m_frameId << endreq;


	// Ask the store gate if the previous algorithm (BlobsFinder --> responsible for clustering)
	//  sent any objects to the StoreGate.
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");
	if(lastObject == 0)
		return;

	// Create the container for track candidates
	m_tracksPtr = new OctoTracksContainer(static_cast<MediPixAlgo *>(this));

	// If so, get the pointer to the last object.  BlobsFinder, as it comes out of the box, sends
	//  a single object containing all clusters.
	m_aB = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject-1);

	// AllBlobsContainer is a box full of Blobs(Clusters). Now we can iterate over all clusters inside.
	vector<blob> blobsVector = m_aB->GetBlobsVector();
	Log << MSG::INFO << "Number of blobs from clustering = " << (Int_t) blobsVector.size() << endreq;
	vector<blob>::iterator blobsItr = blobsVector.begin(); //allBlobs.begin();

	// I will try to join long mip tracks isolated by the cracks between Timepix chips.
	vector<double> slopesV;
	vector<double> cutsV;
	vector<Int_t> blobIndexV;

	int blobIndex = -1;
	for( ; blobsItr != blobsVector.end() ; blobsItr++) {

		blobIndex++;


		// Limit all this to clusters with a minimum size.
		// Note that m_minNPixels can be configured through the Viewer
		//  so you can reprocess and check results online.

		if((*blobsItr).bP.nPixels < m_minClusterSize)
			continue;

		if ( (*blobsItr).GetBlobType() >= _MIP ) {

			// store the angle
			double cangle = (*blobsItr).bP.rotAngle * 180. / TMath::Pi();
			m_segmentsAngles.push_back( cangle );

			//m_cutAngleValue = m_segmentsAngles.back();

	//		if (TMath::Abs( cangle ) > m_cutAngleValue  ) {
	//		continue;
	//		}


			vector< pair<double, double> > limitPixels = MAFTools::FindLimitPixels(*blobsItr);

			float cut = (*blobsItr).bP.fitCut;
			float slope = (*blobsItr).bP.fitSlope;

			float cutOnChip=0;



			//calculate perpendicular distance from each point to trande line in each segments

			set< pair<Int_t, Int_t> > insaidbc = (*blobsItr).GetContentSet();

			Log << MSG::INFO << "Size of segment = " <<(*blobsItr).bP.nPixels<< endreq;

			set< pair<Int_t, Int_t> >::iterator bcI = insaidbc.begin();

			for(;bcI != insaidbc.end();bcI++ ){

				float perpdist_segments=  MAFTools::CalcPerpDistancePointToLineWithSign(slope, cut, *bcI);
				m_perpDistSegments.push_back(perpdist_segments);


				//Log << MSG::INFO << "Perpendicular distance segment = " <<perpdist_segments<< endreq;

			}


			//calculate slope for differences chips on first pixel in chip
			if ((*blobsItr).bP.geoCenter_x <= __chip0_rightedge) {

				cutOnChip = cut;

			} else if ( ((*blobsItr).bP.geoCenter_x > __chip0_rightedge) && ((*blobsItr).bP.geoCenter_x <= __chip1_rightedge)) {

				cutOnChip = slope*__chip0_rightedge + cut;

			} else if ( ((*blobsItr).bP.geoCenter_x > __chip1_rightedge) && ((*blobsItr).bP.geoCenter_x <= __chip2_rightedge)) {

				cutOnChip = slope*__chip1_rightedge + cut;

			} else if ( ((*blobsItr).bP.geoCenter_x > __chip2_rightedge) && ((*blobsItr).bP.geoCenter_x <= __chip3_rightedge)) {

				cutOnChip = slope*__chip2_rightedge + cut;

			} else {

				cutOnChip = cut; // should never happen

			}



			// Gather info for track assembly
			slopesV.push_back(slope);
			cutsV.push_back(cutOnChip);
			blobIndexV.push_back(blobIndex);

			Log << MSG::INFO << "limitPixels head = " << limitPixels[__HEAD_INDX].first << "," << limitPixels[__HEAD_INDX].second << endreq;
			Log << MSG::INFO << "limitPixels tail = " << limitPixels[__TAIL_INDX].first << "," << limitPixels[__TAIL_INDX].second << endreq;

			// Line Corresponding to the linear regression
			MAFTools::DrawLine(this,
					limitPixels[__HEAD_INDX].first,
					limitPixels[__HEAD_INDX].first*slope + cut,
					limitPixels[__TAIL_INDX].first,
					limitPixels[__TAIL_INDX].first*slope + cut,
					2, 2, kRed);

			// Line passing through the farthest pixels
			//MAFTools::DrawLine(this,
			//		limitPixels[__HEAD_INDX].first,
			//		limitPixels[__HEAD_INDX].second,
			//		limitPixels[__TAIL_INDX].first,
			//		limitPixels[__TAIL_INDX].second,
			//		2, 2, kBlack);


		}

	}


	// Map index of the segment with index in the blobs container
	map<int, int> matchAssemblyIndexesMap;
	map<int, double> matchAssemblyCutsMap;

	// Study a possible mip assembly
	if ( (int)slopesV.size() > 1 ) { // more than 1 MIP found

		Log << MSG::INFO << "Multiple mip clusters --> " << slopesV.size() << endreq;

		vector<double>::iterator itrS = slopesV.begin();
		vector<double>::iterator itrC = cutsV.begin();
		vector<int>::iterator itrIndex = blobIndexV.begin();
		vector<int>::iterator itrIndex2 = blobIndexV.begin();

		// Printing the information on the fragments
		Log << MSG::INFO << "List of fragments " << endreq;
		int cntr = 0;
		for ( ; itrS != slopesV.end() ; ) {
			Log << MSG::INFO << "[" << cntr << "] " << "cut : " << *itrC << " | slope : " << *itrS;
			Log << MSG::INFO << " | container index : " << *itrIndex2 << endreq;
			itrS++;
			itrC++;
			cntr++;
			itrIndex2++;
		}



		// Search of a trend.  Rewind iterators.
		vector<double>::iterator itrC2 = cutsV.begin();
		cntr = 0;
		itrIndex2 = blobIndexV.begin();
		int cntr2 = 0, matchesCntr = 0;

		for (itrC = cutsV.begin() ; itrC != cutsV.end() ; itrC++) {

			cntr2=0;
			itrIndex2 = blobIndexV.begin();
			for (itrC2 = cutsV.begin() ; itrC2 != cutsV.end() ; itrC2++) {
				if(itrC == itrC2) {cntr2++; itrIndex2++; continue;} // don't compare same element
				Log << MSG::INFO << "[" << cntr << "] --> " << "[" << cntr2 << "] ";


				if ( (*itrC2) >= (*itrC) - m_joinCutPixelTolerance &&
						(*itrC2) <= (*itrC) + m_joinCutPixelTolerance ) {

					Log << MSG::INFO << " within tolerance ";
					matchesCntr++;
					// Save indexes
					matchAssemblyIndexesMap[cntr] = (*itrIndex);
					matchAssemblyIndexesMap[cntr2] = (*itrIndex2);
					// Save calculated cuts
					matchAssemblyCutsMap[cntr] = (double) (*itrC);
					matchAssemblyCutsMap[cntr2] = (double) (*itrC2);

				}

				Log << MSG::INFO << endreq;
				cntr2++;
				itrIndex2++;
			}
			if (matchesCntr > 1) { // If at least two candidates, finish the search.  Only one assembly per frame taken.
				break;
			} else {
				matchAssemblyIndexesMap.clear();
				matchAssemblyCutsMap.clear();
			}
			matchesCntr = 0; // otherwise rewind
			cntr++;
			itrIndex++;
			Log << MSG::INFO << "--------------------" << endreq;
		}

		//////////////////////////////////////////////////////////////
		// Indexes of the track constituents !!

		if(matchesCntr > 0) {
			Log << MSG::INFO << "Found mip assembly" << endreq;
			map<int, int>::iterator mItr = matchAssemblyIndexesMap.begin();
			map<int, double>::iterator cItr = matchAssemblyCutsMap.begin();
			// Indexes in the container
			vector<int> indexesToStore;
			// Calculated cuts (with Y coordinate) to save per candidate segment
			vector<double> cutsToStore;
			for ( ; mItr != matchAssemblyIndexesMap.end() ; mItr++) {
				Log << MSG::INFO << "[" << (*mItr).first << "] -- in container --> " << (*mItr).second << endreq;
				indexesToStore.push_back( (*mItr).second ); // saving the index
				cutsToStore.push_back( (*cItr).second );    // saving the cut associated
				cItr++; // increase cuts iterator
			}
			m_tracksPtr->PushOneTrackIndxs(indexesToStore);
			m_tracksPtr->PushOneTrackCuts(cutsToStore);
		}

		// Store de number of good segments found
		nSegments = matchesCntr;

	}

	// If an assembly was found let's recalculate the LinearRegression of the big cluster
	if ( !matchAssemblyIndexesMap.empty() ) {

		// New assembly set of pixels
		set<pair<int, int> > completeAssembly;
		// New assembly with alignment
		set<pair<int, int> > completeAssemblyAligned;


		map<int, int>::iterator mItr = matchAssemblyIndexesMap.begin();
		//int totalClusterSize = 0;
		blob segment;
		// we will need the max and min in the X direction.
		int minX = GetMatrixXdim(), maxX = 0;
		for ( ; mItr != matchAssemblyIndexesMap.end() ; mItr++) {
			// Mip segment index
			int msI = (*mItr).second;
			// Cluster segment
			segment = blobsVector[msI];
			Log << MSG::INFO << "Size " << segment.bP.nPixels << endreq;
			// Get set of pixels in this segment
			set<pair<int, int> > temppixelinfo = segment.GetContentSet();
			// Put this information in the complete assembly
			set<pair<int, int> >::iterator ti = temppixelinfo.begin();

			// Preparation for the alignment
			int y_alignment = 0;
			pair<int, int> pix;
			if (segment.bP.geoCenter_x <= __chip0_rightedge){
				y_alignment = 20;

			} else if ( (segment.bP.geoCenter_x > __chip0_rightedge) &&
					(segment.bP.geoCenter_x <= __chip1_rightedge) ){
				y_alignment = 20;

			} else if ( (segment.bP.geoCenter_x > __chip1_rightedge) &&
					(segment.bP.geoCenter_x <= __chip2_rightedge) ){
				y_alignment = 7;

			} else y_alignment = 0;



			for( ; ti != temppixelinfo.end() ; ti++){
				// Build the assembly as it comes
				pix = *ti;
				completeAssembly.insert(pix);


				// save max and minX
				if((*ti).first > maxX) maxX = (*ti).first;
				if((*ti).first < minX) minX = (*ti).first;
				//Log << MSG::INFO << "Pix.second before alignment = " << pix.second << endreq;
				// Build the assembly with corrections
				pix.second += y_alignment; // Y alignment
				//Log << MSG::INFO << "Pix.second after  alignment = " << pix.second << endreq;
				completeAssemblyAligned.insert( pix );
			}
		}



		Log << MSG::INFO << "Complete Assembly size = " << completeAssembly.size() << endreq;

		// Recover the first segment cluster for Linear Regression initialization
		blob hintCluster = blobsVector[ matchAssemblyIndexesMap[0] ];
		// Redo the Linear Regression
		// return values
		double slope, cut, chi2odf, chi2, ndf;
		// init hint values
		double initslope = hintCluster.bP.fitSlope, initcut = hintCluster.bP.fitCut;
		// Perform the linear regression
		MAFTools::LinearRegression(&slope, &cut, &chi2odf, &chi2, &ndf,
				initslope, initcut, completeAssembly);

		Log << MSG::INFO << "Assembly slope = " << slope << endreq;
		Log << MSG::INFO << "Assembly cut = " << cut << endreq;

		// Store the quality of the fit
		m_odf.push_back(chi2odf);

		// Draw the new linear regression line
		MAFTools::DrawLine(this,
				minX,
				minX * slope + cut,
				maxX,
				maxX * slope + cut,
				2, 2, kBlue);


/*
		// Calculate the perpendicular distance of each point on the Assembly to the line
		// Iterate over the big Assembly
		set<pair<int, int> >::iterator assemblyItr = completeAssembly.begin();


		for( ; assemblyItr != completeAssembly.end() ; assemblyItr++){
			pair<int, int> pix = (*assemblyItr);


			//double perpdist = MAFTools::CalcPerpDistanceToLine(slope, cut, pix);

			// calculate perpendicular distance
			double perpdist=  MAFTools::CalcPerpDistancePointToLineWithSign(slope, cut, pix);

			m_perpDist.push_back(perpdist);
			Log << MSG::DEBUG << "Perpendicular distance   = " << perpdist << endreq;

		}
*/

		////////////////////////////////////////////////////////////////////
		// Calculate the residuals for un-aligned track
		// return values
		double slope_2, cut_2, chi2odf_2, chi2_2, ndf_2;
		// init hint values
		double initslope_2 = hintCluster.bP.fitSlope, initcut_2 = hintCluster.bP.fitCut;

		set<pair<int, int> >::iterator assemblyItr = completeAssembly.begin();

		for(; assemblyItr != completeAssembly.end() ;assemblyItr++ ) {


			// Perform the linear regression
			//Log << MSG::DEBUG << " | " << &completeAssembly << endreq;
			MAFTools::LinearRegression(&slope_2, &cut_2, &chi2odf_2, &chi2_2, &ndf_2,
					initslope_2, initcut_2, //fit guess
					completeAssembly,
					*assemblyItr); // Avoid this particular point


			// Calculate first distance, from this point to original Track
			double perpdist_original = MAFTools::CalcPerpDistancePointToLineWithSign(slope, cut, *assemblyItr);
			double perpdist_new = MAFTools::CalcPerpDistancePointToLineWithSign(slope_2, cut_2, *assemblyItr);
			double residual_1 = TMath::Sqrt( TMath::Abs( perpdist_original * perpdist_new ) );
			// Consider negative residual
			if(perpdist_original < 0 || perpdist_new < 0) residual_1 *= -1.;
			Log << MSG::DEBUG << "Unaligned: perpdist_original = " << perpdist_original << "   perpdist_new =  "<<perpdist_new <<endreq;
			Log << MSG::DEBUG << "         residual = " << residual_1 << endreq;

			m_residual.push_back(residual_1);


		}


/*
		////////////////////////////////////////////////////////////////////
		// Calculate the residuals for aligned track
		set<pair<int, int> >::iterator assemblyItrAligned = completeAssemblyAligned.begin();

		// Redo the Linear Regression avoiding one point at a time
		// return values
		double slope_Aligned, cut_Aligned, chi2odf_Aligned, chi2_Aligned, ndf_Aligned,
		slope_Aligned_2, cut_Aligned_2, chi2odf_Aligned_2, chi2_Aligned_2, ndf_Aligned_2;
		// init hint values
		float initslope_2 = hintCluster.bP.fitSlope, initcut_2 = hintCluster.bP.fitCut;


		for( ; assemblyItrAligned != completeAssemblyAligned.end() ; assemblyItrAligned++) {
			// Perform the linear regression
			// For all Aligned assembly
			MAFTools::LinearRegression(&slope_Aligned, &cut_Aligned, &chi2odf_Aligned, &chi2_Aligned, &ndf_Aligned,
					initslope_2, initcut_2, // fit guess
					completeAssemblyAligned);
			// For all Aligned assembly without one point
			MAFTools::LinearRegression(&slope_Aligned_2, &cut_Aligned_2, &chi2odf_Aligned_2, &chi2_Aligned_2, &ndf_Aligned_2,
					initslope_2, initcut_2, // fit guess
					completeAssemblyAligned,
					assemblyItrAligned); // Avoid this particular point

			// Calculate first distance, from this point to original Track
			double perpdist_original = MAFTools::CalcPerpDistancePointToLineWithSign(slope_Aligned, cut_Aligned, *assemblyItrAligned);
			double perpdist_new = MAFTools::CalcPerpDistancePointToLineWithSign(slope_Aligned_2, cut_Aligned_2, *assemblyItrAligned);
			double residual_1 = TMath::Sqrt( TMath::Abs( perpdist_original * perpdist_new ) );
			// Consider negative residual
			if(perpdist_original < 0 || perpdist_new < 0) residual_1 *= -1.;
			Log << MSG::DEBUG << "Aligned: perpdist_original = "<< perpdist_original <<"   perpdist_new =  "<<perpdist_new <<endreq;
			Log << MSG::DEBUG << "         residual = " << residual_1 << endreq;
			m_residualAligned.push_back(perpdist_original);

		}
*/


		// store cluster size
		clusterSize.push_back((int) completeAssembly.size());
	}



	// If an assembly was found
	if ( !matchAssemblyIndexesMap.empty() ) {

		map<int, int>::iterator mItr = matchAssemblyIndexesMap.begin();
		int totalClusterSize = 0;
		for ( ; mItr != matchAssemblyIndexesMap.end() ; mItr++) {
			// Mip segment index
			int msI = (*mItr).second;
			blob ms = blobsVector[msI];
			set< pair<Int_t, Int_t> > bc = ms.GetContentSet();
			set< pair<Int_t, Int_t> >::iterator bcI = bc.begin();

			for( ; bcI != bc.end() ; bcI++) {
				//Log << MSG::INFO << GetMatrixElement(*bcI) << endreq;
				m_assemblyTime.push_back( GetMatrixElement(*bcI) ); // saving time
			}
			// characteristics of the assembly
			float alpha = ms.bP.rotAngle;
			totalClusterSize += ms.bP.nPixels;
			m_alphaAngle.push_back(  alpha*180./TMath::Pi()  );
			m_chi2.push_back(  ms.bP.chisquare_OverDof  );
			m_chi2Prob.push_back( TMath::Prob(ms.bP.chisquare, ms.bP.NDF));
		}
		m_clusterSize.push_back(  totalClusterSize  );

		m_typeOfRecoFrame = __RECO_MIP;
	}
	else // fill it for all in the case of big clusters
	{
		// look at all blobs here
		for(blobsItr = blobsVector.begin() ; blobsItr != blobsVector.end() ; blobsItr++)
		{
			set< pair<Int_t, Int_t> > bc = (*blobsItr).GetContentSet();
			if((int)bc.size() >= m_minNPixelsClouds){
				set< pair<Int_t, Int_t> >::iterator bcI = bc.begin();
				for( ; bcI != bc.end() ; bcI++){
					//Log << MSG::INFO << GetMatrixElement(*bcI) << endreq;
					m_assemblyTime.push_back( GetMatrixElement(*bcI) ); // saving time
				}

			}
		}
		m_typeOfRecoFrame = __RECO_CLOUD;
	}

	// Store track candidates to the StoreGate
	if(m_tracksPtr->GetNumberOfTracks() > 0 ) {
		PullToStoreGateAccess(m_tracksPtr, MPXDefs::DO_NOT_SERIALIZE_ME);
	}

	// Fill the output tree of this algorithm
	getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_distanceToCenter.clear();
	clusterSize.clear();
	m_alphaAngle.clear();
	m_chi2.clear();
	m_chi2Prob.clear();
	m_assemblyTime.clear();
	m_segmentsAngles.clear();
	m_perpDist.clear();
	m_perpDistSegments.clear();
	m_residual.clear();
	m_residualAligned.clear();
	m_odf.clear();
	m_clusterSize.clear();




	nSegments = 0; // for empty frames

}

void OctoCEA::Finalize() {

	Log << MSG::INFO << "Finalize function !" << endreq;

}


/*
 * OctoClusterContainer constructor
 */		// Draw the new linear regression line for Short assembly

ClassImp(OctoTracksContainer)

OctoTracksContainer::OctoTracksContainer(MediPixAlgo * algo) : CandidateContainer(algo) {
	m_numberOfTracks = -1;
}		// Draw the new linear regression line for Short assembly


void OctoTracksContainer::PushOneTrackCuts(vector<double> cuts) {
	// Copy the cuts map which matches in order the indexes
	m_trackCandidatesCuts[m_numberOfTracks] = cuts;
}

void OctoTracksContainer::PushOneTrackIndxs(vector<int> idx) {
	m_numberOfTracks++;
	m_trackCandidatesIndexes[m_numberOfTracks] = idx;
}

#endif
