 /*
 * 	Copyright 2012 John Idarraga, Mathieu Benoit
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

#ifndef BlobsFinder_cpp
#define BlobsFinder_cpp

#include "BlobsFinder.h"
#include "MAFTools/MAFTools.h"

#include "MPXAlgo/Highlighter.h"

using namespace MSG;

ClassImp(BlobsFinder)
ClassImp(AllBlobsContainer)

/* 
 * AllBlobsContainer constructor
 *
 */
AllBlobsContainer::AllBlobsContainer(MediPixAlgo * algo) : 
CandidateContainer(algo) {

}

void AllBlobsContainer::SetBlobType(Int_t id, blobtype type){
	allBlobs[id].SetBlobType(type);

}

blobtype AllBlobsContainer::GetBlobType(Int_t id) { 
	return allBlobs[id].GetBlobType();
}

TString AllBlobsContainer::GetTypeAsString(Int_t id) { 
	return allBlobs[id].GetTypeAsString();
}

void blob::SetBlobType(blobtype type) {
	btype = type;
}

/*
template<class T>
void blob<T>::SetBlobType(T type) {
//void blob<T>::SetBlobType(blobtype type) {
	btype = type;
}
 */

/** 
 * BlobsFinder constructor
 * 
 */
BlobsFinder::BlobsFinder(){

	aB = 0;
	// init parameters
	m_spiralEnds = NEIGHBORS_DISC_0_PIXELS;
	// generally safe to exclude the border line of pixels
	m_excludeBorder = 1;
	// minimum Chi2 required to perform rotation
	m_chisquare_OverDof_rotation = 2.0;
	// no discontinuity in blob's pixel members
	m_discontinuityTolerance = 0;

	// lvl1 if available
	m_lvl1Cut = -1; // -1: no lvl1 available
	// rejected by level1
	m_bP_rejectedByLvl1 = 0;

	// Max Occupancy accepted
	m_maxOcc = 0.5; // 50%
	// This is the maximum size of a single cluster
	m_maxClusterSize = 1000; // In units of pixels

}

/** 
 * Calculate length of the spiral through the formula :
 *  (2*tolerancePixels + 2)*4 + recursive_call(--tolerancePixels)
 *
 *  disc --> length
 *    0        8
 *    1       24
 *    2       48
 *    3       80
 *    and so on ...
 *
 */
Int_t BlobsFinder::GetLengthOfSpiral(Int_t tolerancePixels){

	if(tolerancePixels == -1)
		return 0;
	else
		return (2*tolerancePixels + 2)*4 +
				GetLengthOfSpiral(--tolerancePixels);

}

/** 
 * Set length of the spiral
 *
 */
void BlobsFinder::SetDiscontinuityTolerance(Int_t disc = 0){

	m_discontinuityTolerance = disc;
	m_spiralEnds = GetLengthOfSpiral(disc);

	if(disc >= MAX_DISCONTINUITY)
	{
		Log << MSG::ERROR << "Trying to set a discontinuity value probaly too high ? --> " << disc << endreq;
		Log << MSG::ERROR << "  check BlobsFinder::SetDiscontinuityTolerance." << endreq;
		exit(1);
	}

}

/**
 *  Called at the beginning of the run, before the first frame is
 *  processed.
 *
 */
void BlobsFinder::Init(){

	// a few branches for this algorithm's Tree
	getMyTree()->Branch("nBlobs", &m_bP_nBlobs , "nBlobs/I");
	getMyTree()->Branch("nRejectedLvl1", &m_bP_rejectedByLvl1, "nRejectedLvl1/I");

	getMyTree()->Branch("nPixels", m_bP_nPixels, "nPixels[nBlobs]/I");
	getMyTree()->Branch("clusterSize", m_bP_clusterSize, "clusterSize[nBlobs]/I");
	getMyTree()->Branch("totalCharge", m_bP_totalCharge, "totalCharge[nBlobs]/I");
	getMyTree()->Branch("width_x", m_bP_width_x, "width_x[nBlobs]/I");
	getMyTree()->Branch("width_y", m_bP_width_y, "width_y[nBlobs]/I");
	getMyTree()->Branch("geoCenter_x", m_bP_geoCenter_x, "geoCenter_x[nBlobs]/F");
	getMyTree()->Branch("geoCenter_y", m_bP_geoCenter_y, "geoCenter_y[nBlobs]/F");
	getMyTree()->Branch("weightedCenter_x", m_bP_weightedCenter_x, "weightedCenter_x[nBlobs]/F");
	getMyTree()->Branch("weightedCenter_y", m_bP_weightedCenter_y, "weightedCenter_y[nBlobs]/F");
	getMyTree()->Branch("boxArea", m_bP_boxArea, "boxArea[nBlobs]/I");
	getMyTree()->Branch("circleArea", m_bP_circleArea, "circleArea[nBlobs]/F");
	getMyTree()->Branch("rotAngle", m_bP_rotAngle, "rotAngle[nBlobs]/F");
	getMyTree()->Branch("ellipseArea", m_bP_ellipseArea, "ellipseArea[nBlobs]/F");
	getMyTree()->Branch("balanceToMin", m_bP_balanceToMin, "balanceToMin[nBlobs]/F");
	getMyTree()->Branch("balanceToMax", m_bP_balanceToMax, "balanceToMax[nBlobs]/F");

	// register some values for manipulation from the MPXViewer
	RegisterConfigurationValue(&m_excludeBorder, "border");
	RegisterConfigurationValue(&m_discontinuityTolerance, "discontinuity");
	RegisterConfigurationValue(&m_lvl1Cut, "lvl1cut");
	RegisterConfigurationValue(&m_chisquare_OverDof_rotation, "Chi2OverDof_rot");
	RegisterConfigurationValue(&m_maxOcc, "maxOccupancy");
	RegisterConfigurationValue(&m_maxClusterSize, "maxClusterSize");

	//RegisterConfigurationValue(&m_spiralEnds, "spiral");

}

/**
 *  Called once per frame
 *
 */
void BlobsFinder::Execute(){

	if( GetHitsInPad() == 0 ) { // do nothing if no hits
		return;
	}

	// Check excessive occupancy
	Log << MSG::INFO << "Width  = " << GetWidth() << endreq;
	Log << MSG::INFO << "Height = " << GetHeight() << endreq;
	Log << MSG::INFO << "Hits   = " << GetHitsInPad() << endreq;

	if ( (double) GetHitsInPad() / double ( GetWidth() * GetHeight() ) > m_maxOcc ) {
		return;
	}

	// Get dimensions
	Int_t xDim = GetMatrixXdim();
	Int_t yDim = GetMatrixYdim();

	// Set limits
	limits m_padLimits;
	m_padLimits.colmin = m_excludeBorder;
	m_padLimits.colmax = xDim - m_excludeBorder;
	m_padLimits.rowmin = m_excludeBorder;
	m_padLimits.rowmax = yDim - m_excludeBorder;
	Int_t rowItr = 0;

	// Recalculate lenght or spiral in case discontinuity tolerance
	//  changed.
	m_spiralEnds = GetLengthOfSpiral(m_discontinuityTolerance);

	// initialize One blob data
	ClearOneBlobData();

	// new object containing all the blobs, ready to go to StoreGate when's ready
	aB = new AllBlobsContainer((MediPixAlgo *) this);
	Log << MSG::INFO << "container --> " << aB << endreq;

	//DumpExtendedMask();

	// Loop over the whole matrix except the selected border
	for (Int_t colItr = m_padLimits.colmin ; colItr < m_padLimits.colmax ; colItr++) {

		for (rowItr = m_padLimits.rowmin ; rowItr < m_padLimits.rowmax ; rowItr++) {

			if ( GetMatrixElement(colItr, rowItr) > 0 && ! PixelInExtendedMask(colItr, rowItr, xDim) ) {

				int blobCond = StartBlob(colItr, rowItr);

				if( blobCond == __blob_finished ) {

					// Fill cluster description
					FillClusterDescription();

					// calculate one blob's properties
					int retVal = oneBlob.CalculateProperties(GetFrameWidth(), GetFrameHeight(), m_chisquare_OverDof_rotation);
					DumpProperties();
					//DumpInfo();
					// fill blobs vector ... this guy's going to StoreGate
					if(retVal == __BLOB_GOOD_TO_STORE) aB->allBlobs.push_back(oneBlob);

					// any rejection
					if(retVal == __REJECTED_BY_LEVEL_1) m_bP_rejectedByLvl1++;
					// clean up current blob
					ClearOneBlobData();

				} else if ( blobCond == __blob_too_big_giveupframe ) {
					return;
				}

			}

		}

	}

	Log << MSG::INFO << aB->allBlobs.size() << " blobs found in frame " << GetFrameId() << endreq;

	// pull all blobs to store gage
	PullToStoreGateAccess(aB, MPXDefs::DO_NOT_SERIALIZE_ME);

	// Filling up information to put into this algo's ntuple.
	m_bP_nBlobs = aB->allBlobs.size();
	if(getMyTree()->GetNbranches() > 0) CopyPropertiesToNtupleData();

	// fill ntuple
	getMyTree()->Fill();

	// clean up all blobs
	//aB->allBlobs.clear();

	// rewind
	m_bP_nBlobs = 0;

	// delete the blobs container
	//delete aB;

}

/**
 *  Called at the end of the run, after the last frame is processed.
 *
 */
void BlobsFinder::Finalize(){

}

/**
 *  Single blobs are treated as a transient structure.  Here I clean
 *  up that information after the blob was copied somewhere else.
 *
 */
void BlobsFinder::ClearOneBlobData(){

	// clean up
	oneBlob.blobContent.clear();
	oneBlob.blobContentSet.clear();
	oneBlob.clusterDescription.clear();
	oneBlob.clusterDescriptionToA.clear();
	oneBlob.clusterDescriptionFastToA.clear();

	oneBlob.bP.nPixels = 0;
	oneBlob.bP.clusterSize = 0;
	oneBlob.bP.totalCharge = 0;
	oneBlob.bP.clusterTOT = 0;
	oneBlob.bP.clusterEnergy = 0.;
	oneBlob.bP.width_x = 0;
	oneBlob.bP.width_y = 0;
	oneBlob.bP.geoCenter_x = 0;
	oneBlob.bP.geoCenter_y = 0;
	oneBlob.bP.weightedCenter_x = 0;
	oneBlob.bP.weightedCenter_y = 0;
	oneBlob.bP.boxArea = 0;
	oneBlob.bP.circleArea = 0.;
	oneBlob.bP.minToGeoCenter = 0.;
	oneBlob.bP.maxToGeoCenter = 0.;

	oneBlob.bP.lvl1.clear();

}

/**
 *  Start point of a blob
 *
 */
int BlobsFinder::StartBlob(Int_t x_ini, Int_t y_ini){

	// start point of a blob
	pair<Int_t, Int_t> startPoint = make_pair(x_ini, y_ini);

	/////////////////////////////////////////////////////////////////////////////////////
	// Check if a blob containing this point has already been saved
	//  I only need to find the 'startPoint' in 'bloblContentSet' in each blob
	vector<blob>::iterator blobsItr = aB->allBlobs.begin();
	set< pair<Int_t, Int_t> >::iterator itr;
	for( ; blobsItr != aB->allBlobs.end() ; blobsItr++)
	{
		itr = (*blobsItr).blobContentSet.find(startPoint);
		if(itr != (*blobsItr).blobContentSet.end()) // point found in previous blob !!!
			return false;
	}

	// If we hit this point we are ready to start a new blob
	// Insert in <set>, push_back in <list>, first time
	oneBlob.blobContentSet.insert( startPoint );
	oneBlob.blobContent.push_back( pair< pair<Int_t, Int_t>, Int_t > (startPoint, __CENTER));
	oneBlob.blobItr = oneBlob.blobContent.begin();

	// Follow
	return FollowBlob();

	//	return true;
}

/**
 *  Using recursion to follow a blob up to its end.
 *
 */
int BlobsFinder::FollowBlob(){

	// test insertion in the set
	pair< set < pair<Int_t, Int_t> >::iterator , bool> rtest;

	// spiral starts, i.e. search through the neighbors
	RewindSpiral();

	// spiral starts at
	pair<Int_t, Int_t> nextNeighbor = (*oneBlob.blobItr).first;

	while((*oneBlob.blobItr).second < m_spiralEnds) {

		// try current position, next neighbor
		//nextNeighbor = GetNextPosition( (*oneBlob.blobItr).first );
		nextNeighbor = GetNextPosition( nextNeighbor );

		// see if the next position does not fall in the pad
		if(IsUnsafePosition(nextNeighbor)) {
			// increment neighbor counter, and skip the rest within the while loop
			(*oneBlob.blobItr).second++;
			continue;
		}


		//Log << MSG::LOOP_DEBUG << "next: " <<  nextNeighbor.first << " , " << nextNeighbor.second << " ---> "
		//<< GetMatrixElement(nextNeighbor) << endreq;

		// if there is something
		if(GetMatrixElement(nextNeighbor) > 0)
		{
			Log << MSG::LOOP_DEBUG << "next: " <<  nextNeighbor.first << " , " << nextNeighbor.second << " ---> "
					<< GetMatrixElement(nextNeighbor) << endreq;

			// check the Set see if it existed
			rtest = oneBlob.blobContentSet.insert( nextNeighbor );
			// if insertion is possible it means this point didn't exist.  So, push_back new entry in the list
			if(rtest.second == true) {
				oneBlob.blobContent.push_back( pair< pair<Int_t, Int_t>, Int_t > (nextNeighbor, __CENTER) );
			}

		}

		// increment neighbor counter
		(*oneBlob.blobItr).second++;

		if ( oneBlob.blobContent.size() > m_maxClusterSize ) {
			Log << MSG::ERROR << "This cluster is growing too big .. giving up on the entire frame. " << endreq;
			return __blob_too_big_giveupframe;
		}

	}

	// go to the next in the list, if there are any left
	if ( ++oneBlob.blobItr != oneBlob.blobContent.end() ) {
		// Have a safety switch in case the blob gets too big
		int ret = FollowBlob();
		if ( ret == __blob_too_big_giveupframe ) { return __blob_too_big_giveupframe; }
	}

	return __blob_finished;
}

/**
 *  Check if the position we are about to step in is safe, i.e. see if
 *  it falls inside the pad.
 *
 */
bool BlobsFinder::IsUnsafePosition(pair<Int_t, Int_t> pos){

	if(pos.first < 0 || pos.second < 0)
		return true;
	if(pos.first > GetMatrixXdim()-1 || pos.second > GetMatrixYdim()-1)
		return true;

	return false;
}

/**
 *  Simple blob's contents dump
 *
 */
void BlobsFinder::DumpInfo(){

	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator blobItrLocal = oneBlob.blobContent.begin();
	Int_t totalPixels = 0;
	for( ; blobItrLocal != oneBlob.blobContent.end() ; blobItrLocal++)
	{
		Log << MSG::DEBUG << (*blobItrLocal).first.first << " , " << (*blobItrLocal).first.second
				<< " -- spiral scan --> " << (*blobItrLocal).second << " " << endreq;
		totalPixels++;
	}

	Log << MSG::DEBUG << "---------- " << totalPixels << " pixels in this blob" << endreq;

}

void BlobsFinder::DumpProperties(){

	Log << MSG::DEBUG << "center      : " << oneBlob.bP.geoCenter_x << " , "
			<< oneBlob.bP.geoCenter_y << endreq;
	Log << MSG::DEBUG << "boxArea     : " << oneBlob.bP.boxArea << endreq;
	Log << MSG::DEBUG << "circle      : " << oneBlob.bP.circleArea << endreq;
	Log << MSG::DEBUG << "balance max : " << oneBlob.bP.balanceToMax << endreq;
	Log << MSG::DEBUG << "balance min : " << oneBlob.bP.balanceToMin << endreq;
	Log << MSG::DEBUG << "---------- " << oneBlob.bP.nPixels << " pixels in this blob" << endreq;

}

/**
 *  When done looking around a pixel. Rewind spiral steering
 *  information.
 *
 */
void BlobsFinder::RewindSpiral(){

	m_spiral.xySwitch = X_STEP;
	m_spiral.localMax = 1;
	m_spiral.xCntr = 1;
	m_spiral.yCntr = 1;
	m_spiral.dir = 1;

}

/**
 *  Follows a spiral shape around a pixel
 */
pair<Int_t, Int_t> BlobsFinder::GetNextPosition(pair<Int_t, Int_t> current) {

	pair<Int_t, Int_t> newPos;

	if(m_spiral.xySwitch == X_STEP)
	{
		newPos.first = current.first - m_spiral.dir;
		newPos.second = current.second;
		m_spiral.xCntr++;

		if(m_spiral.xCntr > m_spiral.localMax)
			m_spiral.xySwitch = Y_STEP;
	}
	else if(m_spiral.xySwitch == Y_STEP)
	{
		newPos.second = current.second - m_spiral.dir;
		newPos.first = current.first;
		m_spiral.yCntr++;

		if(m_spiral.yCntr > m_spiral.localMax)
		{
			m_spiral.xySwitch = X_STEP;
			m_spiral.localMax++;
			m_spiral.xCntr = 1;
			m_spiral.yCntr = 1;
			m_spiral.dir *= -1; // change direction
		}
	}

	return newPos;
}

/** 
 * Set variable m_excludeBorder. m_excludeBorder is the number of
 * pixels that I will ignore in the border of the pixel pad.
 *
 */
void BlobsFinder::SetBorderExclusion (Int_t m_ex_i = 1) {

	/*
  if(m_ex_i > GetMatrixXdim()/2 - 1) {
    Log << MSG::ERROR << "Trying to set a border too big ... abort" << endreq;
    exit(1);
  }
	 */
	m_excludeBorder = m_ex_i;

}

////////////////////////////////////////////////////////////
// Cluster

blob::blob (vector<pair<int, int> > v, map<int, int> totmap, int dimX) {

	// init
	ClearOneBlobData();

	vector<pair<int, int> >::iterator i = v.begin();
	pair< pair<int, int>, int > en;
	int tot;

	for ( ; i != v.end() ; i++ ) {


		// Filling  clusterDescription
		tot = totmap[ MAFTools::XYtoC(*i, dimX) ];
		en = make_pair( *i , tot );
		//clusterDescription.push_back ( en );

		cout << (*i).first << " , " << (*i).second << " , " << tot << endl;

		// Filling  blobContentSet
		blobContentSet.insert( *i );

		// Filling  blobContent
		blobContent.push_back ( en );
		// And set the blobItr to the end of the cluster
		blobItr = blobContent.end();

	}

}

/**
 * Make a clone of an existing cluster
 */
blob::blob(const blob & b){

	ClearOneBlobData();

	// structures
	blobContent = b.blobContent;
	blobContentSet = b.blobContentSet;
	clusterDescription = b.clusterDescription;
	clusterDescriptionToA = b.clusterDescriptionToA;
	clusterDescriptionFastToA = b.clusterDescriptionFastToA;
	clusterDescriptionCalibrated = b.clusterDescriptionCalibrated;
	// properties
	bP = b.bP;
	btype = b.btype;

}

/**
 * Make a clone of an existing cluster masking the pixels in "mask"
 */
blob::blob(const blob & b, set<pair<int, int > > mask){

	ClearOneBlobData();

	// blobContent
	list< pair < pair<int, int>, int > >::const_iterator ibc = b.blobContent.begin();
	set<pair<int, int > >::iterator mEnd = mask.end();

	for( ; ibc != b.blobContent.end() ; ibc++) {
		if( mask.find( (*ibc).first ) == mEnd ) {     // value not found --> Ok, copy.
			blobContent.push_back( *ibc );
		} 											  // Otherwise it belongs to the mask and won't be copied
	}

	// Blob content set
	set< pair<int, int> >::const_iterator bcs = b.blobContentSet.begin();
	mEnd = mask.end();

	for( ; bcs != b.blobContentSet.end() ; bcs++){
		if( mask.find( *bcs ) == mEnd ) {             // value not found --> Ok, copy.
			blobContentSet.insert( *bcs );
		} 											  // Otherwise it belongs to the mask and won't be copied
	}

	// Cluster description
	list< pair < pair<int, int>, int > >::const_iterator ibcd = b.clusterDescription.begin();
	mEnd = mask.end();

	for( ; ibcd != b.clusterDescription.end() ; ibcd++) {
		if( mask.find( (*ibcd).first ) == mEnd ) {    // value not found --> Ok, copy.
			clusterDescription.push_back( *ibcd );
		}
	}

	// Cluster description ToA
	list< pair < pair<int, int>, int > >::const_iterator ibcd1 = b.clusterDescriptionToA.begin();
	mEnd = mask.end();

	for( ; ibcd1 != b.clusterDescriptionToA.end() ; ibcd1++) {
		if( mask.find( (*ibcd1).first ) == mEnd ) {    // value not found --> Ok, copy.
			clusterDescriptionToA.push_back( *ibcd1 );
		}
	}

	// Cluster description FastToA
	list< pair < pair<int, int>, int > >::const_iterator ibcd2 = b.clusterDescriptionFastToA.begin();
	mEnd = mask.end();

	for( ; ibcd2 != b.clusterDescriptionFastToA.end() ; ibcd2++) {
		if( mask.find( (*ibcd2).first ) == mEnd ) {    // value not found --> Ok, copy.
			clusterDescriptionFastToA.push_back( *ibcd2 );
		}
	}


	// Cluster description Calibrated
	list< pair < pair<int, int>, int > >::const_iterator ibcdc = b.clusterDescriptionCalibrated.begin();
	mEnd = mask.end();

	for( ; ibcdc != b.clusterDescriptionCalibrated.end() ; ibcdc++) {
		if( mask.find( (*ibcdc).first ) == mEnd ) {    // value not found --> Ok, copy.
			clusterDescriptionCalibrated.push_back( *ibcdc );
		}
	}

	// properties
	CalculateProperties(256, 256, 1);

}

/**
 * Dump all information in the cluster
 */
void blob::AllDump(){

	cout << "--------------------------------------------" << endl;
	// blobContent
	list< pair < pair<int, int>, int > >::iterator ibc = blobContent.begin();
	cout << "[blobContent]                  { list< pair < pair<int, int>, int > > }; size = " << blobContent.size() << " ::: ";
	for( ; ibc != blobContent.end() ; ibc++) {
		cout << "(" << (*ibc).first.first << "," << (*ibc).first.second << ")[" << (*ibc).second << "] , ";
	}
	cout << endl;

	// blobItr
	//cout << "[blobItr] --> (" << (*blobItr).first.first << "," << (*blobItr).first.second << ")[" << (*blobItr).second << "] ( bad entries may mean .end() )" << endl;

	// Blob content set
	set< pair<int, int> >::iterator bcs = blobContentSet.begin();
	cout << "[blobContentSet]               { set< pair<int, int> > };                size = " << blobContentSet.size() << " ::: ";
	for( ; bcs != blobContentSet.end() ; bcs++){
		cout << "(" << (*bcs).first << "," << (*bcs).second << ")    , ";
	}
	cout << endl;

	// Cluster description
	list< pair < pair<int, int>, int > >::iterator ibcd = clusterDescription.begin();
	cout << "[clusterDescription]           { list< pair < pair<int, int>, int > > }; size = " << clusterDescription.size() << " ::: ";
	for( ; ibcd != clusterDescription.end() ; ibcd++) {
		cout << "(" << (*ibcd).first.first << "," << (*ibcd).first.second << ")[" << (*ibcd).second << "] , ";
	}
	cout << endl;

	// Cluster description Calibrated
	list< pair < pair<int, int>, int > >::iterator ibcdc = clusterDescriptionCalibrated.begin();
	cout << "[clusterDescriptionCalibrated] { list< pair < pair<int, int>, int > > }; size = " << clusterDescriptionCalibrated.size() << " ::: ";
	for( ; ibcdc != clusterDescriptionCalibrated.end() ; ibcdc++) {
		cout << "(" << (*ibcdc).first.first << "," << (*ibcdc).first.second << ")[" << (*ibcdc).second << "] , ";
	}
	cout << endl;

	cout << "[Properties]                    clusterTOT  = " << bP.clusterTOT << endl;
	cout << "[ClusterType]                   clusterType = " << btype << endl;

	cout << "--------------------------------------------" << endl;

}

blob & blob::DeletePixelEntry(pair<int, int> p, int padSizex, int padSizey, double chisquare_OverDof_rotation){

	// Blob content set
	list< pair < pair<Int_t, Int_t>, Int_t > > cont = blobContent;
	list< pair < pair<Int_t, Int_t>, Int_t > >::const_iterator itr = cont.begin();
	list< pair < pair<Int_t, Int_t>, Int_t > >::const_iterator iErase = cont.end();
	for ( ; itr != cont.end() ; itr++ ) {
		if ( (*itr).first == p ) {
			iErase = itr;
			break;
		}
	}
	if ( iErase != cont.end() ) blobContent.remove( *iErase );
	/////////////////
	cont = clusterDescription;
	itr = cont.begin();
	iErase = cont.end();
	for ( ; itr != cont.end() ; itr++ ) {
		if ( (*itr).first == p ) {
			iErase = itr;
			break;
		}
	}
	if ( iErase != cont.end() ) clusterDescription.remove( *iErase );
	/////////////////
	cont = clusterDescriptionToA;
	itr = cont.begin();
	iErase = cont.end();
	for ( ; itr != cont.end() ; itr++ ) {
		if ( (*itr).first == p ) {
			iErase = itr;
			break;
		}
	}
	if ( iErase != cont.end() ) clusterDescriptionToA.remove( *iErase );
	/////////////////
	cont = clusterDescriptionFastToA;
	itr = cont.begin();
	iErase = cont.end();
	for ( ; itr != cont.end() ; itr++ ) {
		if ( (*itr).first == p ) {
			iErase = itr;
			break;
		}
	}
	if ( iErase != cont.end() ) clusterDescriptionFastToA.remove( *iErase );

	if ( iErase != cont.end() ) {

		// The first set has a different structure
		set< pair<Int_t, Int_t> >::const_iterator csetItr = blobContentSet.begin();
		set< pair<Int_t, Int_t> >::const_iterator csetItrErase = blobContentSet.end();
		for ( ; csetItr != blobContentSet.end() ; csetItr++ ) {
			if ( *csetItr == p ) {
				csetItrErase = csetItr;
				break;
			}
		}
		blobContentSet.erase( csetItrErase );

	}


	CalculateProperties( padSizex, padSizey, chisquare_OverDof_rotation );

	//this->CalculateProperties(GetFrameWidth(), GetFrameHeight(), 2);

	return *this;
}
void blob::ClearOneBlobData() {

	// clean up
	blobContent.clear();
	blobContentSet.clear();

	bP.nPixels = 0;
	bP.totalCharge = 0;
	bP.width_x = 0;
	bP.width_y = 0;
	bP.geoCenter_x = 0;
	bP.geoCenter_y = 0;
	bP.weightedCenter_x = 0;
	bP.weightedCenter_y = 0;
	bP.boxArea = 0;
	bP.circleArea = 0.;
	bP.minToGeoCenter = 0.;
	bP.maxToGeoCenter = 0.;

	bP.lvl1.clear();

}
#endif
