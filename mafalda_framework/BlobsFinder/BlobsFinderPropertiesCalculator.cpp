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

#include <BlobsFinder.h>
#include <MAFTools/MAFTools.h>

///
/// Calculate properties of one blob once it is complete
///
int blob::CalculateProperties(int xdim, int ydim, double chisquare_OverDof_rotation){

	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator listItr;

	// max an min, x and y in the pad
	Int_t max_x = 0;
	Int_t min_x = xdim-1;
	Int_t max_y = 0;
	Int_t min_y = ydim-1;

	// number of pixels in this blob
	bP.nPixels = clusterDescription.size();
	bP.clusterSize = clusterDescription.size();

	// calculate total charge
	bP.totalCharge = 0; // alias of clusterTOT
	bP.clusterTOT = 0;

	Int_t xcoor = -1;
	Int_t ycoor = -1;
	set< pair<Int_t, Int_t> > blobSet = GetContentSet();
	bP.nInnerPixels = 0;
	Int_t inthits = 0;

	for(listItr = clusterDescription.begin() ; listItr != clusterDescription.end() ; listItr++) {
		// <x, y>
		//bP.totalCharge += GetMatrixElement((*listItr).first);
		bP.totalCharge += (*listItr).second;
		bP.clusterTOT += (*listItr).second; // alias of total charge

		xcoor = (*listItr).first.first;
		ycoor = (*listItr).first.second;

		//std::cout << "x = " << xcoor << " , y = " << ycoor << std::endl;

		if(xcoor  < min_x) min_x = xcoor;
		if(xcoor  > max_x) max_x = xcoor;
		if(ycoor < min_y) min_y  = ycoor;
		if(ycoor > max_y) max_y  = ycoor;

		// Get here the number of inner pixels, cross sape.
		// A pixel is said to be an inner pixel if it is found
		// in the following configuration.
		//         *
		//       * * *
		//         *
		set< pair<Int_t, Int_t> >::iterator itri = blobSet.begin();
		inthits = 0;
		for( ; itri != blobSet.end() ; itri++){
			if((*itri).first == xcoor   && (*itri).second == ycoor+1) inthits++;
			if((*itri).first == xcoor   && (*itri).second == ycoor-1) inthits++;
			if((*itri).first == xcoor-1 && (*itri).second == ycoor  ) inthits++;
			if((*itri).first == xcoor+1 && (*itri).second == ycoor  ) inthits++;
		}

		if(inthits >= MIN_INNERPIXEL){ // internal pixel
			bP.nInnerPixels++;
		}
	}

	// widths of the blob
	bP.width_x = (max_x - min_x) + 1; // it is 1 if the width is 1 pixel
	bP.width_y = (max_y - min_y) + 1; //

	// geometrical center
	bP.geoCenter_x = (max_x + min_x + 1)/2.0;//TMath::CeilNint((max_x + min_x)/2.0);
	bP.geoCenter_y = (max_y + min_y + 1)/2.0;//TMath::CeilNint((max_y + min_y)/2.0);



	// calc min and max distance between all pixels and geoCenter
	CalcMinMaxGeoCenter(bP.geoCenter_x, bP.geoCenter_y);

	// Weighted center
	// Works like the mean of an histogram
	// N = Sum_{i=0}^{M-1} H_{i} // where M is the total number of entries
	//                           // H_{i} is the histogram entry
	// The Mean:  mu = 1/N Sum_{i=0}^{M-1} i H_{i}
	double Nw = 0.;
	double elem = 0.;
	for(listItr = clusterDescription.begin() ; listItr != clusterDescription.end() ; listItr++) {
		// <x, y>
		xcoor = (*listItr).first.first;
		ycoor = (*listItr).first.second;
		// calc N
		//elem = (double)GetMatrixElement((*listItr).first);
		elem = (double) (*listItr).second;
		Nw += elem;
		// calc mu
		bP.weightedCenter_x += (double)xcoor * elem;
		bP.weightedCenter_y += (double)ycoor * elem;

		// the lvl1
		//Int_t lvl1 = GetLVL1((*listItr).first);
		//if(m_lvl1Cut != -1 && lvl1 >= m_lvl1Cut) { // reject this blob
		//	return __REJECTED_BY_LEVEL_1;
		//}
		//bP.lvl1.push_back(GetLVL1((*listItr).first));

	}
	bP.weightedCenter_x /= Nw;
	bP.weightedCenter_y /= Nw;

	// calc box area
	Int_t wx = bP.width_x;
	Int_t wy = bP.width_y;
	bP.boxArea = wx*wy;

	///
	///  Some characteristics are calculated on a rotated cluster and calculated
	///   only if the cluster is not a single, double, triple or quad hit.
	///
	if( SingleHit() || DoubleHit() || TripleHit() || QuadHit() ) {
		bP.circleArea = -1.;
		bP.ellipseArea = -1.;
		bP.ellipseA = -1.;
		bP.ellipseB = -1.;
		bP.rotAngle = -1.;
		bP.chisquare_OverDof = -1;
		bP.fitSlope = -1;
		bP.fitCut = -1;
		bP.balanceToMin = -1;
		bP.balanceToMax = -1;
		return __BLOB_GOOD_TO_STORE;
	}

	// CalculateBalance ... in this algorithm I duplicate the grid to avoid
	//  meaningless information from certain structures.
	//  I also get a guess for the linear regression start parameters.
	double balanceToMin = 0., balanceToMax = 0.,
			slopeGuess = 0., intersectionyGuess = 0.,
			chisquare_OverDof = 0., chisquare = 0., ndf = 0.;
	Int_t minQId = 0, maxQId = 0;
	GetBalance(min_x, max_x, min_y, max_y,
			& balanceToMin, & balanceToMax,
			& minQId, & maxQId,
			& slopeGuess, & intersectionyGuess);

	bP.balanceToMin = balanceToMin;
	bP.balanceToMax = balanceToMax;

	// Calculate linear regression
	double fitSlope = 0.;
	double fitCut = 0.;

	MAFTools::LinearRegression(&fitSlope, &fitCut, &chisquare_OverDof, &chisquare, &ndf,
			slopeGuess, intersectionyGuess,
			GetContentSet());

	bP.chisquare_OverDof = chisquare_OverDof;
	bP.chisquare = chisquare;
	bP.NDF = ndf;
	bP.rotAngle = -1.*TMath::ATan(fitSlope);
	bP.fitSlope = fitSlope;
	bP.fitCut = fitCut;

	// A rotation might be unreliable on a highly pixelized object.  Here's a different approach
	//  to calculate the circle and elipse areas

	// check if the object needs rotation
	if(bP.chisquare_OverDof < chisquare_OverDof_rotation)
	{

		// Find the line passing through the extremes of the blob (using the linear regression is not quite right)
		vector< pair<double, double> > limitPixels = MAFTools::FindLimitPixels(*this);
		double slope = ( limitPixels[1].second - limitPixels[0].second ) / ( limitPixels[1].first - limitPixels[0].first );
		double cut = bP.geoCenter_y - slope*bP.geoCenter_x;

		// Now calculate the perpendicular distance to this line from all points.
		// I will take this as the b of the ellipse
		// Then using the perpline I do the same to find a of the ellipse

		double maxb = 0.;
		double maxa = 0.;
		double dist = 0., distcut = 0.;
		double perpslope = -1./slope;
		double perpcut = bP.geoCenter_y - perpslope*bP.geoCenter_x;
		double cutdistance = MAFTools::CalcDistance(limitPixels[0].first, limitPixels[0].second, limitPixels[1].first, limitPixels[1].second) / bP.width_x;
		double radiusSquared = 0.;

		for (listItr = clusterDescription.begin() ; listItr != clusterDescription.end() ; listItr++) {

			// Tricky here ...
			// To determine b, I scan the points between the line at a perpendicular distance <= 1/5 of max pixels from the perp line
			// the opposite for a

			distcut = MAFTools::CalcPerpDistanceToLine(perpslope, perpcut, (*listItr).first );
			if(distcut <= cutdistance) {
				dist = MAFTools::CalcPerpDistanceToLine(slope, cut, (*listItr).first );
				if(dist > maxb) maxb = dist;
			}

			distcut = MAFTools::CalcPerpDistanceToLine(slope, cut, (*listItr).first );
			if(distcut <= cutdistance) {
				dist = MAFTools::CalcPerpDistanceToLine(perpslope, perpcut, (*listItr).first );
				if(dist > maxa) maxa = dist;
			}

			// If this didn't work the calculation if a and b simplifies to
			// This happens with small cross like (5-hit) clusters.
			if(maxa == 0.) maxa = bP.width_x;
			if(maxb == 0.) maxb = bP.width_y;

		}

		bP.ellipseA = maxa/2.;
		bP.ellipseB = maxb/2.;

		// wx/2 is a, wy/2 is b.  AreaElipse = Pi*a*b
		bP.ellipseArea = TMath::Pi()*bP.ellipseA*bP.ellipseB;

		// fraction --> (sqrt(2)/2)
		if(bP.ellipseA >= bP.ellipseB)
			radiusSquared = bP.ellipseA * bP.ellipseA;
		else
			radiusSquared = bP.ellipseB * bP.ellipseB;

		bP.circleArea = TMath::Pi()*radiusSquared;

	}
	else // no rotation needed
	{
		Double_t radiusSquared = (wx*wx + wy*wy) * 0.5/4.; // fraction --> (sqrt(2)/2)
		bP.circleArea = TMath::Pi()*radiusSquared;
		bP.ellipseA = (Float_t)wx / 2.;
		bP.ellipseB = (Float_t)wy / 2.;
		bP.ellipseArea = TMath::Pi() * bP.ellipseA * bP.ellipseB;
	}

	/*
	// check if the object needs rotation
	if(oneBlob.bP.chisquare_OverDof > m_chisquare_OverDof_rotation)
	{

		// now rotate every point in the blob and find new max and min
		Int_t rot_max_x = -2*GetMatrixXdim()-1;
		Int_t rot_min_x = 2*GetMatrixXdim()-1;
		Int_t rot_max_y = -2*GetMatrixYdim()-1;
		Int_t rot_min_y = 2*GetMatrixYdim()-1;

		pair<Int_t, Int_t> rotPair;

		//Int_t totalItr = 0;
		for(listItr = oneBlob.blobContent.begin() ; listItr != oneBlob.blobContent.end() ; listItr++)
		{
			// <rot_x, rot_y>
			rotPair = RotateVector2D((*listItr).first.first, (*listItr).first.second, oneBlob.bP.rotAngle);
			if(rotPair.first  < rot_min_x) rot_min_x = rotPair.first;
			if(rotPair.first  > rot_max_x) rot_max_x = rotPair.first;
			if(rotPair.second < rot_min_y) rot_min_y = rotPair.second;
			if(rotPair.second > rot_max_y) rot_max_y = rotPair.second;

			//totalItr++;
		}

		Int_t rot_wx = (rot_max_x - rot_min_x) + 1;
		Int_t rot_wy = (rot_max_y - rot_min_y) + 1;

		oneBlob.bP.ellipseA = rot_wx;
		oneBlob.bP.ellipseB = rot_wy;

		// wx/2 is a, wy/2 is b.  AreaElipse = Pi*a*b
		oneBlob.bP.ellipseArea = TMath::Pi()*(rot_wx*rot_wy)/4.0;

		// fraction --> (sqrt(2)/2)
		Double_t radiusSquared = 0.;
		if(rot_wx >= rot_wy)
			radiusSquared = (2.*rot_wx*rot_wx) * (0.5/4.);
		else
			radiusSquared = (2.*rot_wy*rot_wy) * (0.5/4.);

		oneBlob.bP.circleArea = TMath::Pi()*radiusSquared;

	}
	else // no rotation needed
	{
		Double_t radiusSquared = (wx*wx + wy*wy) * 0.5/4.; // fraction --> (sqrt(2)/2)
		oneBlob.bP.circleArea = TMath::Pi()*radiusSquared;
		oneBlob.bP.ellipseArea = TMath::Pi()*(wx * wy)/4.0;
		oneBlob.bP.ellipseA = wx;
		oneBlob.bP.ellipseB = wy;
	}
	 */

	// finally type initialized to _NOTYPE_ASSIGNED
	// PRBasicSpecies will decide later
	btype = _NOTYPE_ASSIGNED;

	//max_x = 0;
	//min_x = GetMatrixXdim()-1;
	//max_y = 0;
	//min_y = GetMatrixYdim()-1;
	return __BLOB_GOOD_TO_STORE;
}

//TString blob::Dump(){
//}

pair<Int_t, Int_t> BlobsFinder::RotateVector2D(Int_t x, Int_t y, double rotAngle){

	pair<Int_t, Int_t> rotPair;

	rotPair.first =
			TMath::FloorNint( (double)x*TMath::Cos(rotAngle) - (double)y*TMath::Sin(rotAngle) );

	rotPair.second =
			TMath::FloorNint( (double)x*TMath::Sin(rotAngle) + (double)y*TMath::Cos(rotAngle) );

	return rotPair;
}


double blob::GetBalance(Int_t min_x, Int_t max_x, Int_t min_y, Int_t max_y,
		double * balanceToMin, double * balanceToMax,
		Int_t * minQId, Int_t * maxQId,
		double * slope, double * intersectiony){

	double mean = 0.;
	Int_t divvalmax = 0, divvalmin = 0;

	mean = CalcMeanForBlob(min_x, max_x, min_y, max_y, &divvalmax, &divvalmin, minQId, maxQId, slope, intersectiony);

	// calculate the width to get the size of one quadrant, extended to
	// the longest side
	Int_t width_x = (max_x - min_x) + 1; // it is 1 if the width is 1 pixel
	Int_t width_y = (max_y - min_y) + 1; //
	if(width_x < width_y)
		width_x = width_y;

	// The balance is the ratio of the minimum occupancy to the area
	// of one quadrant as calculated before.
	*balanceToMin = (double)divvalmin / (double)(width_x*width_x);
	*balanceToMax = (double)divvalmax / mean;

	// return mean, probably not used
	return mean;
}

double blob::CalcMeanForBlob(Int_t min_x, Int_t max_x, Int_t min_y, Int_t max_y,
		Int_t * divvalmax, Int_t * divvalmin,
		Int_t * minQId, Int_t * maxQId,
		double * slope, double * intersectiony){

	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator listItr;
	Int_t posx = -1, posy = -1;
	Int_t newposx = -1, newposy = -1;
	double mean = 0.;

	/////////////////////////////////////////////////////////////
	// Working in a doubled grid !
	min_x *= 2;
	max_x = 2*max_x + 1;
	min_y *= 2;
	max_y = 2*max_y + 1;
	// recalculate centers
	Int_t c_x = TMath::CeilNint((double)(min_x + max_x)/2.0);
	Int_t c_y = TMath::CeilNint((double)(min_y + max_y)/2.0);
	/////////////////////////////////////////////////////////////

#define NQUADRANTS 4

	Int_t firstQ = 0;
	Int_t secondQ = 0;
	Int_t thirdQ = 0;
	Int_t fourthQ = 0;
	Int_t maxQ = 0; // maximum numbers of counts between the four quadrants
	Int_t minQ = 0; // minimum numbers of counts between the four quadrants

	// I need to calculate the mean x and y coordinates in the four cuadrants
	double meansQuadX[] = {0., 0., 0., 0.};
	double meansQuadY[] = {0., 0., 0., 0.};

	for(listItr = clusterDescription.begin() ; listItr != clusterDescription.end() ; listItr++)
	{

		posx = (*listItr).first.first;
		posy = (*listItr).first.second;

		// Now for each pixel in the list I get 4

		for(Int_t i = 0 ; i < 2 ; i++)
		{
			for(Int_t j = 0 ; j < 2 ; j++)
			{
				// 4 positions in the doubled-grid system
				newposx = 2*posx + i;
				newposy = 2*posy + j;

				// now separate quadrants
				if(newposx - c_x < 0) // 3rd or 4th quadrant
				{
					if(newposy - c_y >= 0){ // 4th
						fourthQ++;
						meansQuadX[Q_FOURTH_INDX] += (double)newposx;
						meansQuadY[Q_FOURTH_INDX] += (double)newposy;
					}
					else if(newposy - c_y < 0){ // 3rd
						thirdQ++;
						meansQuadX[Q_THIRD_INDX] += (double)newposx;
						meansQuadY[Q_THIRD_INDX] += (double)newposy;
					}
				}
				else if(newposx - c_x >= 0) // 1st or 2nd quadrant
				{
					if(newposy - c_y >= 0){ // 1st
						firstQ++;
						meansQuadX[Q_FIRST_INDX] += (double)newposx;
						meansQuadY[Q_FIRST_INDX] += (double)newposy;
					}
					else if(newposy - c_y < 0){ // 2nd
						secondQ++;
						meansQuadX[Q_SECOND_INDX] += (double)newposx;
						meansQuadY[Q_SECOND_INDX] += (double)newposy;
					}
				}
			}
		}

	}

	// means per quadrant
	if(fourthQ) meansQuadX[Q_FOURTH_INDX] /= (double)fourthQ;
	else meansQuadX[Q_FOURTH_INDX] = c_x - 1;  // if there is nothing in
	// a cuadrant, take the
	// center +/-1 depending
	// on the cuadrant.  This
	// will at least give a
	// correct guess for the
	// slope and ycut

	if(fourthQ) meansQuadY[Q_FOURTH_INDX] /= (double)fourthQ;
	else meansQuadY[Q_FOURTH_INDX] = c_y + 1;

	if(thirdQ) meansQuadX[Q_THIRD_INDX] /= (double)thirdQ;
	else meansQuadX[Q_THIRD_INDX] = c_x - 1;

	if(thirdQ) meansQuadY[Q_THIRD_INDX] /= (double)thirdQ;
	else meansQuadY[Q_THIRD_INDX] = c_y - 1;

	if(firstQ) meansQuadX[Q_FIRST_INDX] /= (double)firstQ;
	else meansQuadX[Q_FIRST_INDX] = c_x + 1;

	if(firstQ) meansQuadY[Q_FIRST_INDX] /= (double)firstQ;
	else meansQuadY[Q_FIRST_INDX] = c_y + 1;

	if(secondQ) meansQuadX[Q_SECOND_INDX] /= (double)secondQ;
	else meansQuadX[Q_SECOND_INDX] = c_x + 1;

	if(secondQ) meansQuadY[Q_SECOND_INDX] /= (double)secondQ;
	else meansQuadY[Q_SECOND_INDX] = c_y - 1;

	// find max
	if(firstQ > maxQ){
		maxQ = firstQ;
		*maxQId = Q_FIRST; // Id of the quadrant
	}
	if(secondQ > maxQ){
		maxQ = secondQ;
		*maxQId = Q_SECOND;
	}
	if(thirdQ > maxQ){
		maxQ = thirdQ;
		*maxQId = Q_THIRD;
	}
	if(fourthQ > maxQ){
		maxQ = fourthQ;
		*maxQId = Q_FOURTH;
	}
	// This could be the only option in a perfectly symmetric blob
	minQ = maxQ;
	*minQId = *maxQId;

	// find min
	if(firstQ < minQ){
		minQ = firstQ;
		*minQId = Q_FIRST;
	}
	if(secondQ < minQ){
		minQ = secondQ;
		*minQId = Q_SECOND;
	}
	if(thirdQ < minQ){
		minQ = thirdQ;
		*minQId = Q_THIRD;
	}
	if(fourthQ < minQ){
		minQ = fourthQ;
		*minQId = Q_FOURTH;
	}

	// calc mean
	mean = (firstQ + secondQ + thirdQ + fourthQ)/4.0;

	// decide the sign of the slope and calculate
	if(*maxQId == Q_FOURTH || *maxQId == Q_SECOND) { // negative slope
		*slope = (meansQuadY[Q_FOURTH_INDX] - meansQuadY[Q_SECOND_INDX])/
				(meansQuadX[Q_FOURTH_INDX] - meansQuadX[Q_SECOND_INDX]);

		*intersectiony = meansQuadY[Q_SECOND_INDX] - ( (*slope) * meansQuadX[Q_SECOND_INDX] );
		// convert to the simple grid (this working-grid was doubled !)
		*intersectiony /= 2.;
	}
	else if (*maxQId == Q_FIRST || *maxQId == Q_THIRD) { // positive slope
		*slope = (meansQuadY[Q_FIRST_INDX] - meansQuadY[Q_THIRD_INDX])/
				(meansQuadX[Q_FIRST_INDX] - meansQuadX[Q_THIRD_INDX]);

		*intersectiony = meansQuadY[Q_THIRD_INDX] - ( (*slope) * meansQuadX[Q_THIRD_INDX] );
		// convert to the simple grid (this working-grid was doubled !)
		*intersectiony /= 2.;

	}

	//Log << MSG::DEBUG << "guess slope = " << *slope << endreq;
	//Log << MSG::DEBUG << "guess ycut  = " << *intersectiony << endreq;

	// other passed-by-reference values
	*divvalmax = maxQ;
	*divvalmin = minQ;

	return mean;
}

/**
 * Copy properties of all blobs in one frame for ntuple dump
 *
 */
void BlobsFinder::CopyPropertiesToNtupleData(){

	if(aB->allBlobs.size() > MAX_BLOBS_PER_FRAME){
		Log << MSG::ERROR << "Too many blobs found, and can't save to the ntuple.  Currently = " << aB->allBlobs.size() << endreq;
		Log               << "You need to change MAX_BLOBS_PER_FRAME and recompile this algorithm." << endreq;
		Log               << "Or, decide not to save to ntuple at all." << endreq;
		exit(1);
	}

	vector<blob>::iterator blobsItr = aB->allBlobs.begin();

	if(m_bP_nBlobs != (Int_t) aB->allBlobs.size()){ // castint  uint --> int.  Safe.
		Log << MSG::ERROR << "Something went really bad here, check BlobsFinder::CalculateProperties()" << endreq;
		exit(1);
	}

	Int_t itr = 0;
	for( ; blobsItr != aB->allBlobs.end() ; blobsItr++)
	{

		m_bP_nPixels[itr] = (*blobsItr).bP.nPixels;
		m_bP_clusterSize[itr] = (*blobsItr).bP.nPixels;
		m_bP_totalCharge[itr] = (*blobsItr).bP.totalCharge;

		m_bP_width_x[itr] = (*blobsItr).bP.width_x;
		m_bP_width_y[itr] = (*blobsItr).bP.width_y;

		m_bP_geoCenter_x[itr] = (*blobsItr).bP.geoCenter_x;
		m_bP_geoCenter_y[itr] = (*blobsItr).bP.geoCenter_y;

		m_bP_weightedCenter_x[itr] = (*blobsItr).bP.weightedCenter_x;
		m_bP_weightedCenter_y[itr] = (*blobsItr).bP.weightedCenter_y;

		m_bP_boxArea[itr] = (*blobsItr).bP.boxArea;
		m_bP_circleArea[itr] = (*blobsItr).bP.circleArea;
		m_bP_ellipseArea[itr] = (*blobsItr).bP.ellipseArea;

		m_bP_balanceToMin[itr] = (*blobsItr).bP.balanceToMin;
		m_bP_balanceToMax[itr] = (*blobsItr).bP.balanceToMax;

		itr++;
	}

}

TString blob::GetTypeAsString(){

	if(btype == _NOTYPE_ASSIGNED)
		return TString("Not_Assigned");
	else if(btype == _SINGLE_HIT)
		return TString("Single_Hit");
	else if(btype == _DOUBLE_HIT)
		return TString("Double_Hit");
	else if(btype == _TRIPLE_HIT)
		return TString("Triple_Hit");
	else if(btype == _QUAD_HIT)
		return TString("Quad_Hit");
	else if(btype == _LONG_GAMMA)
		return TString("Long_Gamma");
	else if(btype == _MIP)
		return TString("Mip");
	else if(btype == _HEAVY_BLOB)
		return TString("Heavy_Blob");
	else if(btype == _CURLY)
		return TString("Curly");
	else if(btype == _UNKNOWN)
		return TString("Unknown");

	return TString("Unknown");
}

///
/// Most basic cluster shapes
///
/// Single :  *
bool blob::SingleHit(){

	if(bP.nPixels == 1) {
		btype = _SINGLE_HIT;
		return true;
	}

	return false;
}

/// Double :  ** || *  || *  ||  *
///                 *      *    *
bool blob::DoubleHit(){

	if(
			(bP.nPixels == 2 && bP.width_x == 2 && bP.width_y == 2)
			||
			(bP.nPixels == 2 && bP.width_x == 2 && bP.width_y == 1)
			||
			(bP.nPixels == 2 && bP.width_x == 1 && bP.width_y == 2)
	) {
		btype = _DOUBLE_HIT;
		return true;
	}

	return false;
}

/// Triple :  ** || ** || *  ||  *
///            *    *     **    **
bool blob::TripleHit(){

	if(bP.nPixels == 3 && bP.width_x == 2 && bP.width_y == 2){
		btype = _TRIPLE_HIT;
		return true;
	}
	return false;
}

/// Quad   :  **
///           **
bool blob::QuadHit(){

	if(bP.nPixels == 4 && bP.width_x == 2 && bP.width_y == 2) {
		btype = _QUAD_HIT;
		return true;
	}

	return false;
}

void BlobsFinder::FillClusterDescription(){

	list< pair < pair<Int_t, Int_t>, Int_t > > bc = oneBlob.GetBlobContent();
	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator i = bc.begin();

	for ( ; i != bc.end() ; i++) {

		oneBlob.clusterDescription.push_back( pair< pair<Int_t, Int_t>, Int_t > ( (*i).first, GetMatrixElement( (*i).first) ) );
		//Log << MSG::ALWAYS << "+ " << GetToA( (*i).first ) << endreq;
		oneBlob.clusterDescriptionToA.push_back( pair< pair<Int_t, Int_t>, Int_t > ( (*i).first, GetToA( (*i).first ) ) );
		oneBlob.clusterDescriptionFastToA.push_back( pair< pair<Int_t, Int_t>, Int_t > ( (*i).first, GetFastToA( (*i).first ) ) );

	}

}

void blob::CalcMinMaxGeoCenter(double geox, double geoy){

	Int_t xcoor = -1;
	Int_t ycoor = -1;
	double dist = 0.;
	// start point
	xcoor = (*(clusterDescription.begin())).first.first;
	ycoor = (*(clusterDescription.begin())).first.second;
	dist = TMath::Sqrt((xcoor - geox)*(xcoor - geox)
			+
			(ycoor - geoy)*(ycoor - geoy));
	bP.minToGeoCenter = dist;
	bP.maxToGeoCenter = dist;
	// loop over all elements
	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator listItr;
	for(listItr = clusterDescription.begin() ; listItr != clusterDescription.end() ; listItr++)
	{
		// <x, y>
		xcoor = (*listItr).first.first;
		ycoor = (*listItr).first.second;

		dist = TMath::Sqrt((xcoor - geox)*(xcoor - geox)
				+
				(ycoor - geoy)*(ycoor - geoy));
		if(dist < bP.minToGeoCenter)
			bP.minToGeoCenter = dist;
		else if(dist > bP.maxToGeoCenter)
			bP.maxToGeoCenter = dist;
	}

}
