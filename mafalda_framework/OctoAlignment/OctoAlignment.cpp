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

#ifndef __OctoAlignment_cpp
#define __OctoAlignment_cpp

#include "OctoAlignment.h"
#include "MAFTools.h"

using namespace MSG;

ClassImp(OctoAlignment)

OctoAlignment::OctoAlignment(){

}

void OctoAlignment::Init(){

	Log << MSG::INFO << "Init function !" << endreq;

	//getMyTree()->Branch("distanceToCenter", &m_distanceToCenter);
	getMyTree()->Branch("angleDifferenses1_2",&m_angleDifferences1_2);
	getMyTree()->Branch("angleDifferenses2_3",&m_angleDifferences2_3);
	getMyTree()->Branch("angleDifferenses3_4",&m_angleDifferences3_4);
	getMyTree()->Branch("angleDifferenses1_3",&m_angleDifferences1_3);
	getMyTree()->Branch("angleDifferenses2_4",&m_angleDifferences2_4);
	getMyTree()->Branch("angleDifferenses1_4",&m_angleDifferences1_4);
	getMyTree()->Branch("cutDifferences1_2",&m_cutDifferences1_2);
	getMyTree()->Branch("cutDifferences2_3",&m_cutDifferences2_3);
	getMyTree()->Branch("cutDifferences3_4",&m_cutDifferences3_4);
	getMyTree()->Branch("cutDifferences1_3",&m_cutDifferences1_3);
	getMyTree()->Branch("cutDifferences2_4",&m_cutDifferences2_4);
	getMyTree()->Branch("cutDifferences1_4",&m_cutDifferences1_4);
	// A configuration value that can be tuned from the Viewer
	//RegisterConfigurationValue(&m_minNPixels, "minNPixels");

}

void OctoAlignment::Execute(){

	Log << MSG::INFO << "Frame Id = " << GetFrameId() << endreq;
	//////////////////////////////////////////////////////////////////
	// 1) Retrieve the clusters container
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");
	// If no object registered by BlobsFinder, get off !
	if(lastObject == 0) return;

	m_clusterContainerPtr = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject-1);

	// AllBlobsContainer is a box full of Blobs(Clusters). Now we can iterate over all clusters inside.
	vector<blob> blobsVector = m_clusterContainerPtr->GetBlobsVector();
	Log << MSG::INFO << "Number of blobs from clustering = " << (Int_t) blobsVector.size() << endreq;
	// If non clusters, get off !
	if((Int_t) blobsVector.size() == 0) return;
	vector<blob>::iterator blobsItr = blobsVector.begin(); //allBlobs.begin();

	//////////////////////////////////////////////////////////////////
	// 2) Retrieve the track candidate (segments) container
	lastObject = GetNumberOfObjectsWithAuthor("OctoCEA");
	// If no object registered by OctoCEA, get off !
	if(lastObject == 0) return;

	m_candidatesContainerPtr = (OctoTracksContainer *) GetObjectFromAuthor("OctoCEA", lastObject-1);
	Log << MSG::INFO << "Found cluster container at    : " << m_clusterContainerPtr << endreq;
	Log << MSG::INFO << "Found candidates container at : " << m_candidatesContainerPtr << endreq;


	// Print constituents of tracks
	int nOfTracks = m_candidatesContainerPtr->GetNumberOfTracks();
	Log << MSG::INFO << "Number of tracks reconstructed = " << nOfTracks << endreq;

	// If no tracks, get off !
	if(nOfTracks == 0 || nOfTracks > 5) return;

	// Get the information of the tracks
	for(int track_ind = 0 ; track_ind < nOfTracks  ; track_ind++){

		vector<int> segments_vector = m_candidatesContainerPtr->GetOneTrackIndxs(track_ind);
		vector<double> cuts_vector = m_candidatesContainerPtr->GetOneTrackCuts(track_ind);

		if (segments_vector.size()<2)
			return;

		Log << MSG::INFO << "Track N. " << track_ind << " has " <<segments_vector.size() << " segments <";
		vector<int>::iterator itr = segments_vector.begin();
		for(;itr!=segments_vector.end();itr++){
			cout<<*itr<< " ";
		}
		Log << MSG::INFO << ">" << endreq;

		//Work with geo centers

		set<pair<float_t, float_t> > GeoCenters;

		if (segments_vector.size()>2){
			for(itr = segments_vector.begin();itr!=segments_vector.end();itr++){

				 ;

				Log << MSG::INFO << "Geo Center =  " << blobsVector[*itr].bP.geoCenter_x<<"  "<<blobsVector[*itr].bP.geoCenter_y<<endreq;

			}

		}


		//Print Angles for segments
		for(itr=segments_vector.begin() ; itr!=segments_vector.end() ; itr++) {

			double angle = blobsVector[*itr].bP.rotAngle*180/TMath::Pi();
			Log << MSG::INFO << "Angles :"<< angle <<endreq;

		}
		//Print Cuts for segments
		vector<double>::iterator itrC = cuts_vector.begin();
		Log << MSG::INFO << "cuts vector size = " << cuts_vector.size() << endreq;
		for( ; itrC != cuts_vector.end() ; itrC++) {

			double cut = *itrC;
			Log << MSG::INFO << "Cut :"<< cut << endreq;

		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// Calculate angle and cut differences for all pair of candidate to track

		//Calculate angle and cut differences
		double previous_angle = blobsVector[*segments_vector.begin()].bP.rotAngle*180/TMath::Pi();
		double previous_cut = *(cuts_vector.begin());

		itrC = (cuts_vector.begin() + 1);
		for( itr = (segments_vector.begin() + 1) ; itr!=segments_vector.end() ; itr++) {

			//Log << MSG::INFO << "[" << *itr << "]"<< " cuts : " << *itrC << endreq;

			// Calculate the angle and cut for segment *itr
			double angle = blobsVector[*itr].bP.rotAngle*180/TMath::Pi();
			double cut = *itrC;

			// Calculate the differences
			double angle_difference = angle - previous_angle;
			previous_angle = angle;
			double cut_differents=previous_cut - cut;
			previous_cut=cut;

			if ( (blobsVector[*itr].bP.geoCenter_x>256) && (blobsVector[*itr].bP.geoCenter_x<512)&&
					(blobsVector[*(itr-1)].bP.geoCenter_x<=256)){ // segment 2
				m_cutDifferences1_2.push_back(cut_differents);
				m_angleDifferences1_2.push_back(angle_difference);
			}

			if ( (blobsVector[*itr].bP.geoCenter_x>512) && (blobsVector[*itr].bP.geoCenter_x<768)&&
					(blobsVector[*(itr-1)].bP.geoCenter_x>256)&&(blobsVector[*(itr-1)].bP.geoCenter_x<=512)){ // segment 3
				m_cutDifferences2_3.push_back(cut_differents);
				m_angleDifferences2_3.push_back(angle_difference);
			}

			if ( (blobsVector[*itr].bP.geoCenter_x>768) && (blobsVector[*itr].bP.geoCenter_x<1024)&&
					(blobsVector[*(itr-1)].bP.geoCenter_x>512)&&(blobsVector[*(itr-1)].bP.geoCenter_x<=768) ){ // segment 4
				m_cutDifferences3_4.push_back(cut_differents);
				m_angleDifferences3_4.push_back(angle_difference);
			}

			// incrementing the second iterator
			itrC++;
		}


		// Third to first, and fourth to second differences
		itrC = (cuts_vector.begin() + 2);
		for(itr=(segments_vector.begin()+2) ; itr!=segments_vector.end() ; itr++) {

			double previous_angle =blobsVector[*(itr-2)].bP.rotAngle*180/TMath::Pi();
			double previous_cut=*(itrC-2);

			double angle = blobsVector[*itr].bP.rotAngle*180/TMath::Pi();
			double cut = *itrC;

			double angle_difference = angle - previous_angle;
			double cut_differents=previous_cut - cut;


			if ( (blobsVector[*itr].bP.geoCenter_x>512) && (blobsVector[*itr].bP.geoCenter_x<768)&&
					(blobsVector[*(itr-2)].bP.geoCenter_x<=256) ){
				m_cutDifferences1_3.push_back(cut_differents);
				m_angleDifferences1_3.push_back(angle_difference);
			}

			if ( (blobsVector[*itr].bP.geoCenter_x>768) && (blobsVector[*itr].bP.geoCenter_x<1024)&&
					(blobsVector[*(itr-2)].bP.geoCenter_x>256)&&(blobsVector[*(itr-2)].bP.geoCenter_x<=512) ){
				m_cutDifferences2_4.push_back(cut_differents);
				m_angleDifferences2_4.push_back(angle_difference);
			}

			itrC++;
		}

		// for 1 and 4th candidate
		if (segments_vector.size()>3){

			//for(itr=(segments_vector.begin()) ; itr!=segments_vector.end() ; itr++) {

			itr = (segments_vector.begin()); // step on the very first segment
			itrC = cuts_vector.begin();

			double previous_angle =blobsVector[*itr].bP.rotAngle*180/TMath::Pi();
			double previous_cut=*itrC;

			double angle = blobsVector[*(itr+3)].bP.rotAngle*180/TMath::Pi();
			double cut = *(itrC+3);

			double angle_difference = angle - previous_angle;
			double cut_differents=previous_cut - cut;

			//Log << MSG::INFO << "Cut differences 1_4  :"<<cut_differents  <<endreq;
			//Log << MSG::INFO << "Geo center  :"<<blobsVector[*(itr+3)].bP.geoCenter_x <<" size :" << segments_vector.size() <<endreq;
			//Log << MSG::INFO << "Geo center  :"<<blobsVector[*segments_vector.end()].bP.fitCut <<" size :" << segments_vector.size() <<endreq;

			if ( (blobsVector[*(itr+3)].bP.geoCenter_x>768) && (blobsVector[*(itr+3)].bP.geoCenter_x<1024)&&
					(blobsVector[*(itr)].bP.geoCenter_x<=256) ){
				m_cutDifferences1_4.push_back(cut_differents);
				m_angleDifferences1_4.push_back(angle_difference);

			}
			//break;
			//}

		}


	}



	// Fill the output tree of this algorithm
	getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_angleDifferences1_2.clear();
	m_angleDifferences2_3.clear();
	m_angleDifferences3_4.clear();
	m_angleDifferences1_3.clear();
	m_angleDifferences2_4.clear();
	m_angleDifferences1_4.clear();
	m_cutDifferences1_2.clear();
	m_cutDifferences2_3.clear();
	m_cutDifferences3_4.clear();
	m_cutDifferences1_3.clear();
	m_cutDifferences2_4.clear();
	m_cutDifferences1_4.clear();
}

void OctoAlignment::Finalize() {

	Log << MSG::INFO << "Finalize function !" << endreq;

}

#endif
