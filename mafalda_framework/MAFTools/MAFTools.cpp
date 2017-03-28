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


#ifndef MAFTools_cpp
#define MAFTools_cpp

#include "MAFTools.h"
#include "MPXAlgo/Highlighter.h"
#include "MPXAlgo/MediPixAlgo.h"

#include <set>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>

#include <stdlib.h>

#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>
#include <TGraph2D.h>

using namespace std;

namespace MAFTools {

Bool_t LinearRegression(double * slope, double * cut,
		double * chisquare_OverDof, double * chisquare, double * ndf,
		double initslope, double initycross,
		vector< pair<Int_t, Int_t> > blobContent) {

	// copy the vector into a set and do the normal thing
	set<pair<int, int> > bc;
	vector<pair<int, int> >::iterator i = blobContent.begin();
	for( ; i != blobContent.end() ; i++){
		bc.insert( *i );
	}

	return LinearRegression(slope, cut, chisquare_OverDof, chisquare, ndf,
			initslope, initycross, bc);

}

Bool_t LinearRegression(double * slope, double * cut,
		double * chisquare_OverDof, double * chisquare, double * ndf,
		double initslope, double initycross,
		set< pair<Int_t, Int_t> > blobContent,
		pair<Int_t, Int_t> pairAvoid) {

	// Create a copy of blobContent avoiding the point "pairAvoid"
	set<pair<int, int> > truncatedBlogContent;
	set<pair<int, int> >::iterator i = blobContent.begin();
	int newnentries = 0;
	for ( ; i != blobContent.end() ; i++ ) {
		if ( *i != pairAvoid ) {
			truncatedBlogContent.insert( *i );
			newnentries++;
		}
	}

	cout << "Attempting to take out one element.  Original size = " << blobContent.size()
  													<< " | new size = " << truncatedBlogContent.size() << endl;

	//cout << "New set : " << endl;
	//set<pair<int, int> >::iterator i2;
	//for(i2 = truncatedBlogContent.begin() ; i2 != truncatedBlogContent.end() ; i2++) {
	//	cout << "(" << (*i2).first << "," << (*i2).second << " | ";
	//}
	//cout << endl;

	// And proceed
	return LinearRegression(slope, cut, chisquare_OverDof, chisquare, ndf,
			initslope, initycross, truncatedBlogContent);
}


Bool_t LinearRegression(double * slope, double * cut,
		double * chisquare_OverDof, double * chisquare, double * ndf,
		double initslope, double initycross,
		set< pair<Int_t, Int_t> > blobContent) {


	// data -->N points
	Int_t N = (Int_t) blobContent.size();

	// Fit to ax + b.  N has to be greater or equal to the number of
	// linearly indepentend functions, in this case 2.
	if(N < 2) {
		// fit failed
		*slope = 0;
		*chisquare_OverDof = 0.;
		return false;
	}

	double * xd = new double[N];
	double * yd = new double[N];
	Int_t cntr = 0;

	set< pair<Int_t, Int_t> >::iterator setItr = blobContent.begin();
	for( ; setItr != blobContent.end() ; setItr++) {
		xd[cntr] = (double)(*setItr).first;
		yd[cntr] = (double)(*setItr).second;
		cntr++;
	}

	// get the fit done, linear fit
	TF1 * f1 = new TF1("f1","[0]*x + [1]",0. , 1.);
	// enter the guess for the parameters
	f1->SetParameters(initslope, initycross);
	TGraph * g1 = new TGraph(N, xd, yd);
	g1->Fit(f1,"NQ");

	*chisquare_OverDof = f1->GetChisquare()/(double)f1->GetNDF();
	*chisquare = f1->GetChisquare();
	*ndf = (double)f1->GetNDF();
	*slope = f1->GetParameter(0);
	*cut = f1->GetParameter(1);

	delete f1;
	delete g1;
	delete [] xd;
	delete [] yd;

	return true;
}

TGraph * GetTGraphZeroSuppression(TGraph2D * g2D){

	int N = g2D->GetN();
	double * x = g2D->GetX();
	double * y = g2D->GetY();
	double * z = g2D->GetZ();
	// get the number of non zero entries
	int NnonZero = 0;
	for ( int i = 0 ; i < N ; i++ ) {
		if ( z[i] > 0. ) NnonZero++;
	}
	// Copy them to a subset
	double * x_nz = new double[NnonZero];
	double * y_nz = new double[NnonZero];
	int nz_cntr = 0;
	for ( int i = 0 ; i < N ; i++ ) {
		if ( z[i] > 0. ) {
			x_nz[nz_cntr] = x[i];
			y_nz[nz_cntr] = y[i];
			nz_cntr++;
		}
	}
	// create the TGraph
	TGraph * g = new TGraph(NnonZero, x_nz, y_nz);

	return g;
}

double TriangleArea(double * x, double * y) {

	double x0 = x[0];
	double y0 = y[0];
	double x1 = x[1];
	double y1 = y[1];
	double x2 = x[2];
	double y2 = y[2];

	double base = CalcDistance(x0, y0, x1, y1);
	double h = 0.;
	// check the slope of the parallel won't be infinite
	if ( y1-y0 == 0. ) {
		h = TMath::Abs( y2 - y1 );
	} else {
		// h is the perp distance between the line (x0,y0)-->(x1,y1) and the point (x2,y2)
		double slope = (y1 - y0) / (x1 - x0);
		double cut = y1 - slope*x1;
		h = CalcPerpDistanceToLine(slope, cut, x2, y2);
	}

	return h*base/2.;
}

Int_t XYtoC(pair<Int_t, Int_t> position, Int_t dimX) {
	return XYtoC(position.first, position.second, dimX);
}

Int_t XYtoC(int x, int y, Int_t dimX){
	return y * dimX + x;
}

Int_t XYtoX(int x, int y, Int_t dimX){
	return y * dimX + x;
}

Int_t XYtoX(pair<int, int> pix, Int_t dimX){
	return pix.second * dimX + pix.first;
}

Int_t ConvertXYPositionToInteger(pair<Int_t, Int_t> position, Int_t dimX, Int_t dimY){

	return (position.first)%dimX + (position.second)*dimY;

}

Int_t ConvertXYPositionToInteger(Int_t positionX, Int_t positionY, Int_t dimX, Int_t dimY){

	pair<Int_t, Int_t> coordinates (positionX, positionY);

	return ConvertXYPositionToInteger(coordinates, dimX, dimY);

}

pair<int, int> XtoXY(int X, int dimX){
	return make_pair(X % dimX, X/dimX);
}

Double_t CalcDistance(Double_t x1, Double_t y1,
		Double_t x2, Double_t y2){

	return TMath::Sqrt((y2 - y1)*(y2 - y1)
			+
			(x2 - x1)*(x2 - x1));

}

Double_t CalcDistance(std::pair<Double_t, Double_t> p1, std::pair<Double_t, Double_t> p2) {
	return CalcDistance(p1.first, p1.second, p2.first, p2.second);
}

Double_t CalcDistance(std::pair<int, int> p1, std::pair<int, int> p2) {
	return CalcDistance((double)p1.first, (double)p1.second, (double)p2.first, (double)p2.second);
}

///
/// Distance to line can be positive or negative in this form
///
Double_t CalcPerpDistancePointToLineWithSign(Float_t slope, Float_t cut, pair<Int_t, Int_t> point){
	return CalcPerpDistancePointToLineWithSign(slope, cut, point.first, point.second);
}
Double_t CalcPerpDistancePointToLineWithSign(Float_t slope, Float_t cut, Int_t px, Int_t py){
	//calculate dictance from point to line
	Float_t distance = (slope*px - py + cut)/sqrt(slope*slope +1);
	return distance;
}

Double_t CalcPerpDistanceToLine(double slope, double cut, pair<Int_t, Int_t> point){
	return CalcPerpDistanceToLine(slope, cut, point.first, point.second);
}

Double_t CalcPerpDistanceToLine(double slope, double cut, Int_t px, Int_t py){

	double cutPrim = 0.;
	double slopePrim = -1./slope;

	// Cut of the perpendicular line
	cutPrim  = py;
	cutPrim -= slopePrim * px;

	// intersection of the two lines
	double interx = (cutPrim - cut)/(slope - slopePrim);
	double intery = interx * slope + cut;

	return CalcDistance(px, py, interx, intery);
}

void DrawLocalCoordSystem(MediPixAlgo * algo, blob b){

	vector< pair<double, double> > limitPixels = MAFTools::FindLimitPixels(b);
	limitPixels[0].first -= 10; // little offset for graphics
	limitPixels[1].first += 10; // little offset for graphics

	double cut = b.bP.fitCut;
	double slope = b.bP.fitSlope;

	MAFTools::DrawLine(algo,
			limitPixels[1].first,
			limitPixels[1].first*slope + cut,
			limitPixels[0].first,
			limitPixels[0].first*slope + cut,
			2, 1, kRed);

	double pslope = -1./slope;
	double pcut = b.bP.geoCenter_y - (b.bP.geoCenter_x * pslope);

	MAFTools::DrawLine(algo,
			(b.bP.geoCenter_x - 1),
			(b.bP.geoCenter_x - 1)*pslope + pcut,
			(b.bP.geoCenter_x + 1),
			(b.bP.geoCenter_x + 1)*pslope + pcut,
			2, 1, kGreen);

}

void DrawLine(MediPixAlgo * algo, double x1, double y1, double x2, double y2, Int_t wd = 2, Int_t st = 2, EColor color = kBlack) {

	Highlighter * lineD = new Highlighter(x1, y1, x2, y2, "line", algo);
	lineD->SetLineColor(color);
	lineD->SetLineWidth(wd);
	lineD->SetLineStyle(st);
	algo->PullToStoreGateAccess(lineD, MPXDefs::DO_NOT_SERIALIZE_ME);

}

void DrawHull(MediPixAlgo * algo, vector< pair<int, int> > hc, int lcolor, float lshift){

	vector < pair<int, int> >::iterator setItr = hc.begin();
	vector < pair<int, int> >::iterator setItr2 = hc.begin();

	// keep first point
	pair<int, int> firstp = *setItr;
	double x1, y1, x2, y2;
	for ( ; setItr != hc.end() ; setItr++) {

		setItr2++;

		x1 = (*setItr).first;
		y1 = (*setItr).second;

		if (setItr2 != hc.end()) {
			x2 = (*setItr2).first;
			y2 = (*setItr2).second;
		} else {
			x2 = firstp.first;
			y2 = firstp.second;
		}

		Highlighter * lineD = new Highlighter(x1 + lshift, y1, x2 + lshift, y2, "line", algo);
		lineD->SetLineColor(lcolor);
		//lineD->SetLineWidth(wd);
		//lineD->SetLineStyle(st);
		algo->PullToStoreGateAccess(lineD, MPXDefs::DO_NOT_SERIALIZE_ME);

	}

}

std::string DumpWordInBinary(UInt_t word) {

	string wordS;

	UInt_t oneBit = 0x0;
	Int_t sizeOfWord = sizeof(word)*8; // number of bits

	for (int i = 0 ; i < sizeOfWord ; i++) { // number of bits

		oneBit = word & 0x1;

		if(oneBit == 0x1) { wordS.append(1,'1'); }
		else { wordS.append(1,'0'); }

		word = word >> 1;
	}

	reverse(wordS.begin(),wordS.end());

	return wordS;
}

vector<pair<double, double> > FindLimitPixels(blob theBlob) {

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
			}

		}

	}

	return limits;
}


Bool_t BlobAtADistanceFromEdges(blob, vector<pair<double, double> > l, Int_t width, Int_t height,
		Int_t cutx, Int_t cuty) {

	Int_t dist = 0; // in units of pixels

	// lower border
	for (int i = 0 ; i < (int)l.size() ; i++) {
		dist = TMath::FloorNint(l[i].second);
		if ( dist <= cuty ) return true;
	}
	// left border
	for (int i = 0 ; i < (int)l.size() ; i++) {
		dist =  TMath::FloorNint(l[i].first);
		if ( dist <= cutx ) return true;
	}
	// top border
	for (int i = 0 ; i < (int)l.size() ; i++) {
		dist = height - TMath::FloorNint(l[i].second);
		if ( dist <= cuty ) return true;
	}
	// right border
	for (int i = 0 ; i < (int)l.size() ; i++) {
		dist = width - TMath::FloorNint(l[i].first);
		if ( dist <= cutx ) return true;
	}

	return false;
}

int FindDiscontinuity(blob b){

	// get the content set x,y
	set< pair<Int_t, Int_t> > cs = b.GetContentSet();
	set< pair<Int_t, Int_t> >::iterator csItr = cs.begin();
	set< pair<Int_t, Int_t> >::iterator csItr2 = cs.begin();

	double dist = 0.;
	int nClose = 0;
	for( ; csItr != cs.end() ; csItr++) {

		nClose = 0;
		for(csItr2 = cs.begin() ; csItr2 != cs.end() ; csItr2++) {

			if(csItr != csItr2) {
				dist = CalcDistance(*csItr, *csItr2);
			}

			if(dist < 1.5) nClose++;

		}
		if (nClose == 1) return 1; // discontinuity
		//std::cout << endl;
		//std::cout << "nClose = " << nClose << " ---------" << std::endl;

	}

	return 0;
}


/**
 *  All clustering stuff
 *
 */
bool StartBlob(int xi, int yi, Int_t disctTolerance, AllBlobsContainer * blobsContainer, blob * oneBlob, int maxx, int maxy, std::map<int,int> totmap){

	// I will recluster here on a cluster already found (if discontinuity > 0).
	// I will have to rewind the internal iterator so the search is possible again.
	//(*oneBlob.blobItr).second = 0;

	pair<Int_t, Int_t> startPoint = make_pair(xi, yi);

	//cout << "--> " << xi << " , " << yi << endl;

	/////////////////////////////////////////////////////////////////////////////////////
	// Check if a blob containing this point has already been saved
	//  I only need to find the 'startPoint' in 'bloblContentSet' in each blob
	vector<blob>::iterator blobsItr = blobsContainer->allBlobs.begin();
	set< pair<Int_t, Int_t> >::iterator itr;
	for( ; blobsItr != blobsContainer->allBlobs.end() ; blobsItr++)
	{
		itr = (*blobsItr).blobContentSet.find(startPoint);
		if(itr != (*blobsItr).blobContentSet.end()) // point found in previous blob !!!
			return false;
	}

	// If we hit this point we are ready to start a new blob
	// Insert in <set>, push_back in <list>, first time
	oneBlob->blobContentSet.insert( startPoint );
	oneBlob->blobContent.push_back( pair< pair<int, int>, int > (startPoint, __CENTER));
	oneBlob->blobItr = oneBlob->blobContent.begin();

	// Follow
	steerSpiral spiral;
	FollowBlob(oneBlob, &spiral, disctTolerance, maxx, maxy, totmap);

	return true;
}

/**
 *  Using recursion to follow a blob up to its end.
 *
 */
bool FollowBlob(blob * oneBlob, steerSpiral * spiral, int disctTolerance, int maxx, int maxy, std::map<int,int> totmap){

	//cout << " +++ looking at = " << (*(oneBlob->blobItr)).first.first << ", " << (*(oneBlob->blobItr)).first.second << endl;

	Int_t spiralEnds = GetLengthOfSpiral(disctTolerance);

	// test insertion in the set
	pair< set < pair<Int_t, Int_t> >::iterator , bool> rtest;

	// spiral starts, i.e. search through the neighbors
	RewindSpiral(spiral);

	// spiral starts at
	pair<Int_t, Int_t> nextNeighbor = (*oneBlob->blobItr).first;

	(*oneBlob->blobItr).second = 0;
	//cout << " +++ iterator = " << (*oneBlob->blobItr).second << ", " << spiralEnds << endl;

	while((*oneBlob->blobItr).second < spiralEnds) {

		// try current position, next neighbor
		//nextNeighbor = GetNextPosition( (*oneBlob.blobItr).first );
		nextNeighbor = GetNextPosition( nextNeighbor, spiral );

		//cout << "next: " <<  nextNeighbor.first << " , " << nextNeighbor.second << " : " << totmap[ nextNeighbor.second * maxx + nextNeighbor.first ] << endl;

		// see if the next position does not fall in the pad
		/*
		if(IsUnsafePosition(nextNeighbor, maxx, maxy)) {
			// increment neighbor counter, and skip the rest withing the while loop
			(*oneBlob.blobItr).second++;
			continue;
		}
		 */
		//Log << MSG::LOOP_DEBUG << "next: " <<  nextNeighbor.first << " , " << nextNeighbor.second << " ---> "
		//<< GetMatrixElement(nextNeighbor) << endreq;

		// if there is something
		if(totmap[ nextNeighbor.second * maxx + nextNeighbor.first ] != 0)
		{
			//Log << MSG::LOOP_DEBUG << "next: " <<  nextNeighbor.first << " , " << nextNeighbor.second << " ---> "
			//		<< GetMatrixElement(nextNeighbor) << endreq;

			//cout << "next: " <<  nextNeighbor.first << " , " << nextNeighbor.second << " ---> "
			//<< totmap[ nextNeighbor.second * maxx + nextNeighbor.first ] << endl;

			// check the Set see if it existed
			rtest = oneBlob->blobContentSet.insert( nextNeighbor );
			// if insertion is possible it means this point didn't exist.  So, push_back new entry in the list
			if(rtest.second == true) {
				oneBlob->blobContent.push_back( pair< pair<int, int>, int > (nextNeighbor, __CENTER) );
			}

		}

		// increment neighbor counter
		(*oneBlob->blobItr).second++;
	}

	// go to the next in the list, if there are any left
	if(++oneBlob->blobItr != oneBlob->blobContent.end()) {
		FollowBlob(oneBlob, spiral, disctTolerance, maxx, maxy, totmap);
	}

	return true;
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
Int_t GetLengthOfSpiral(Int_t tolerancePixels){

	if(tolerancePixels == -1)
		return 0;
	else
		return (2*tolerancePixels + 2)*4 +
				GetLengthOfSpiral(--tolerancePixels);

}

/**
 *  Follows a spiral shape around a pixel
 *
 */
pair<Int_t, Int_t> GetNextPosition(pair<Int_t, Int_t> current, steerSpiral * spiral) {

	pair<Int_t, Int_t> newPos;

	if(spiral->xySwitch == X_STEP)
	{
		newPos.first = current.first - spiral->dir;
		newPos.second = current.second;
		spiral->xCntr++;

		if(spiral->xCntr > spiral->localMax)
			spiral->xySwitch = Y_STEP;
	}
	else if(spiral->xySwitch == Y_STEP)
	{
		newPos.second = current.second - spiral->dir;
		newPos.first = current.first;
		spiral->yCntr++;

		if(spiral->yCntr > spiral->localMax)
		{
			spiral->xySwitch = X_STEP;
			spiral->localMax++;
			spiral->xCntr = 1;
			spiral->yCntr = 1;
			spiral->dir *= -1; // change direction
		}
	}

	return newPos;
}

void RewindSpiral(steerSpiral * spiral){

	spiral->xySwitch = X_STEP;
	spiral->localMax = 1;
	spiral->xCntr = 1;
	spiral->yCntr = 1;
	spiral->dir = 1;

}

void FillValuesForDisplay(Highlighter * hl, blob bl, int nextras, vector<TString> names, vector<double> extras){

	// First fill the properties in the blob and return if no extras
	FillValuesForDisplay(hl, bl);
	if(nextras == 0) return;

	for (int i = 0 ; i < nextras ; i++) {
		hl->UploadDisplayValue(names[i].Data(), extras[i]);
	}

}


void FillValuesForDisplay(Highlighter * hl, blob bl){

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

}

/**
 *  Check if the position we are about to step in is safe, i.e. see if
 *  it falls inside the pad.
 *
 */
bool IsUnsafePosition(pair<Int_t, Int_t> pos, int maxx, int maxy){

	if(pos.first < 0 || pos.second < 0)
		return true;
	if(pos.first > maxx-1 || pos.second > maxy-1)
		return true;

	return false;
}

void ClearOneBlobData(blob * oneBlob){

	// clean up
	oneBlob->blobContent.clear();
	oneBlob->blobContentSet.clear();

	oneBlob->bP.nPixels = 0;
	oneBlob->bP.totalCharge = 0;
	oneBlob->bP.width_x = 0;
	oneBlob->bP.width_y = 0;
	oneBlob->bP.geoCenter_x = 0;
	oneBlob->bP.geoCenter_y = 0;
	oneBlob->bP.weightedCenter_x = 0;
	oneBlob->bP.weightedCenter_y = 0;
	oneBlob->bP.boxArea = 0;
	oneBlob->bP.circleArea = 0.;
	oneBlob->bP.minToGeoCenter = 0.;
	oneBlob->bP.maxToGeoCenter = 0.;

	oneBlob->bP.lvl1.clear();

}

int DigitalBlobReclustering(int discontinuityTolerance, AllBlobsContainer * bContainer, blob bl,
		int maxx, int maxy, std::map<int,int> totmap){

	// Start point of a blob
	// This clustering runs on the contents of an existing blob
	set< pair<Int_t, Int_t> > bCSet = bl.GetContentSet();
	set< pair<Int_t, Int_t> >::iterator bCSetItr = bCSet.begin();

	// Blob structure to find
	blob * oneBlob = new blob;
	ClearOneBlobData(oneBlob);

	for ( ; bCSetItr != bCSet.end() ; bCSetItr++) {

		int xi = (*bCSetItr).first;
		int yi = (*bCSetItr).second;

		if( StartBlob(xi, yi, discontinuityTolerance, bContainer, oneBlob, maxx, maxy, totmap) ) {
			bContainer->allBlobs.push_back(*oneBlob);
			ClearOneBlobData(oneBlob);
		}

	}

	delete oneBlob;

	return 0;
}

// 3D representation of a blob with z being the TOT

bool GoodCandidateGraph2D(blob b){

	// If the cluster is 1D line this operation is not
	//  well defined.
	if(b.bP.width_x == 1 || b.bP.width_y == 1) {
		return false;
	}

	return true;
}

blob BuildBlob(TGraph2D * g){

	// The new blob
	blob b;

	// Get information from TGraph2D
	int N = g->GetN();
	double * x = g->GetX();
	double * y = g->GetY();
	double * z = g->GetZ();

	set<pair<int,int> > cs; // BlobContentSet
	list< pair < pair<int, int>, int > > c; // BlobContent
	list< pair < pair<int, int>, int > > cd; // Cluster description
	for(int i = 0 ; i < N ; i++){
		// The cluster in TGraph2D contains pixels with zeros
		// Get rid of those here
		if(z[i] == 0) continue;
		cs.insert( make_pair(x[i], y[i]) );
		c.push_back(  make_pair( make_pair(x[i], y[i]),    0 )  ) ; // the second part is a clustering ordering, doesn't mean anything here
		cd.push_back(  make_pair( make_pair(x[i], y[i]), z[i] )  ) ; // This one really contains the TOT info
	}

	// insert the content set
	b.SetContentSet(cs);
	b.SetBlobContent(c);
	b.SetClusterDescription(cd);

	return b;
}

///
/// Get a slice following the Linear regression of the blob --> dir = "X"
/// or in the perpendicular direction --> dir = "Y"

TH1 * ExtractSlice(blob b, int eid, TString dir){

#define __sl_x  100
#define __sl_y  200

	TString htitle = "slice_";
	htitle += eid;
	htitle += "_";
	htitle += dir;
	htitle += "_";
	htitle += b.bP.geoCenter_x;
	htitle += "_";
	htitle += b.bP.geoCenter_y;

	int dirv = __sl_x; // default
	if(dir == "Y") dirv = __sl_y;

	// Needed cluster information
	double slope = b.bP.fitSlope;
	double cut = b.bP.fitCut;
	double perpslope = -1./slope;
	double perpcut = b.bP.geoCenter_y - (perpslope * b.bP.geoCenter_x);

	list< pair < pair<int, int>, int > > cdes = b.GetClusterDescription();
	list< pair < pair<int, int>, int > >::iterator i = cdes.begin();

	// Get the closest points to the Linear Regression and
	// calculate the min and max in X
	double dist = 0.;
	list< pair < pair<int, int>, int > > closepts;
	int min = 0x0fffffff;
	int max = 0;
	for ( ; i != cdes.end() ; i++){

		// Depending on the shape this histos may have more points in X that in Y and viceversa
		// For instance take a comet-shape ... X(along the linear reggresion) more points
		// satisfying the condition will be found
		if(dirv == __sl_x) dist = CalcPerpDistanceToLine(slope, cut, (*i).first );
		if(dirv == __sl_y) dist = CalcPerpDistanceToLine(perpslope, perpcut, (*i).first );

		if( dist < 1.0 ) { // inside one pixel
			closepts.push_back( make_pair( (*i).first , (*i).second ) );
			if( (*i).first.first < min ) min = (*i).first.first;
			if( (*i).first.first > max ) max = (*i).first.first;
		}
	}

	// I need to make a map to pile up the averages first
	// X-projected , AverageCounts
	map<int, double> averages;
	map<int, int> entries;

	for (i = closepts.begin() ; i != closepts.end() ; i++){
		// just add up for the average
		averages[ (*i).first.first ] += (*i).second;
		// keep track of the number of entries
		if( entries.find( (*i).first.first ) == entries.end() ) { // if it's not in the map then push
			entries[ (*i).first.first ] = 1;
		} else {
			entries[ (*i).first.first ]++; // add if there is a new entry
		}
	}
	// Calculate the averages
	map<int, double>::iterator itra;
	for(itra = averages.begin() ; itra != averages.end() ; itra++){
		//cout << (*itra).second << " | " << entries[(*itra).first] << endl;
		(*itra).second /= (double)entries[(*itra).first];
	}
	// Finally fill the histogram
	//cout << "----------" << endl;
	TH1F * h = new TH1F(htitle, htitle, max-min, min, max);
	for(itra = averages.begin() ; itra != averages.end() ; itra++){
		//cout << (*itra).second << endl;
		h->Fill((*itra).first, (*itra).second);
	}

	return static_cast<TH1 *>(h);
}

int IdendifyPeaks_A1 (TH1 * h, vector<double> * criticalPoints) {

	queue<int> ent;
	int qsize = 3;
	int npeaksguess = 0;

	int Nbins = h->GetNbinsX();
	if(Nbins < 5) return 0; // Nothing to look for here

	double slope = 0.;
	int sign, prevsign, nPos = 0, nNeg = 0;
	int postimer = 0, negtimer = 0;

	// init queue --> qsize - 1
	for (int i = 1 ; i <= qsize ; i++) {
		ent.push( h->GetBinContent( i ) );
	}
	SimpleLinearReggresion(ent, &slope);
	prevsign = GetSign(slope);
	//cout << "init sign = " << prevsign << " | ";

	if(prevsign == 1) postimer++;
	if(prevsign == -1) negtimer++;

	// continue
	for (int i = qsize+1 ; i <= Nbins ; i++) {

		// Get rid of oldest, get ready for this iteration
		ent.pop();

		// push a new
		ent.push( h->GetBinContent( i ) );
		SimpleLinearReggresion(ent, &slope);

		// decide on sign and flip
		sign = GetSign(slope);

		if (sign >= 0 and prevsign == -1) { // critical point

			// first check previous condition
			if(negtimer > qsize-1) {
				nNeg++;
				criticalPoints->push_back( i ); // saving bin
			}

			postimer++;   // start pos timer
			negtimer=0;   // stop neg timer

		}
		if (sign <= 0 and prevsign == 1) { // critical point

			// first check previous condition
			if(postimer > qsize-1) {
				nPos++;
				criticalPoints->push_back( i ); // saving bin
			}

			postimer=0; // stop pos timer
			negtimer++; // start neg timer

		}

		// Keep timers runing
		if (sign ==  1 and postimer > 0) postimer++;
		if (sign == -1 and negtimer > 0) negtimer++;

		// keep sign
		prevsign = sign;

		//cout << sign << "[" << postimer << "," << negtimer << "] , ";

	}

	// Check final timers
	if(postimer > qsize-1) nPos++;
	if(negtimer > qsize-1) nNeg++;


	if(nPos == nNeg) npeaksguess = nPos;
	if(nPos == nNeg-1 || nPos-1 == nNeg) npeaksguess = SmallMin(nPos, nNeg);

	//cout << " ||| nPos = " << nPos << " , nNeg = " << nNeg << endl;

	return npeaksguess;
}

int SmallMin(int a, int b){
	if(a < b) return a;
	if(b < a) return b;
	return a;
}

int GetSign(double slope){
	if(slope > 0) return 1;
	if(slope < -1) return -1;
	if(slope == 0) return 0;
	return 0;
}

///
/// This is only good for the peaks identification.
/// The cut information here is meaningless
///
void SimpleLinearReggresion(queue<int> q, double * slope){

	int N = (int) q.size();
	double x = 0;
	vector<pair<double, double> > q2;
	// build the x,y pairs
	while(!q.empty()) {
		// retrieve oldest element
		q2.push_back( make_pair( x , (double) q.front() ) );
		// and remove it
		q.pop();
		x += 1; // next x value ... arbitrary
	}

	vector<pair<double, double> >::iterator i = q2.begin();
	vector<double> xy, x2;
	double sumx =  0.;
	double sumy =  0.;
	double sumxy = 0.;
	double sumx2 =  0.;
	for( ; i != q2.end() ; i++) {
		xy.push_back( (*i).first * (*i).second );
		x2.push_back( (*i).first * (*i).first );
		sumx  += (*i).first;
		sumy  += (*i).second;
		sumxy += (*i).first * (*i).second;
		sumx2 += (*i).first * (*i).first;
	}

	// slope
	*slope = ( N*sumxy - sumx*sumy ) / ( N*sumx2 - sumx2*sumx2);

}

TGraph * ClusterProfile(blob b, long frameId, TString dir, TString type) {

	// first get the histo
	TH2F * h = CropCluster(b, frameId, type);

	int nx = h->GetNbinsX();
	int ny = h->GetNbinsY();

	vector<double> x;
	vector<double> y;
	for(int i = 1 ; i <= nx ; i++) {
		// take a bin in x
		double columnadd = 0.;
		for(int j = 1 ; j <= ny ; j++) {
			// and add up all y
			columnadd += h->GetBinContent(i,j);
		}
		x.push_back(i);
		y.push_back(columnadd);
	}

	TGraph * g = new TGraph((int)x.size(), &x[0], &y[0]);
	TString hname = "g_profile_";
	hname += b.bP.geoCenter_x; hname += "_"; hname += b.bP.geoCenter_y;
	hname += "_"; hname += frameId;
	g->SetName(hname);

	return g;
}

vector<TGraph *> GetAllProfiles(blob b, long frameId, TString dir, TString type) {

	vector<TGraph *> graphs;

	// first get the histo
	TH2F * h = CropCluster(b, frameId, type, "extra");

	int nx = h->GetNbinsX();
	int ny = h->GetNbinsY();

	vector<double> x;
	vector<double> y;
	for(int i = 1 ; i <= nx ; i++) {
		// take a bin in x
		for(int j = 1 ; j <= ny ; j++) {
			// build each profile per each column
			y.push_back( h->GetBinContent(i,j) );
			x.push_back(j);
		}

		TGraph * g = new TGraph((int)x.size(), &x[0], &y[0]);
		TString hname = "g_profile_";
		hname += b.bP.geoCenter_x; hname += "_"; hname += b.bP.geoCenter_y;
		hname += "_"; hname += frameId; hname += "_sliceY_"; hname += i;
		g->SetName(hname);
		graphs.push_back( g );

		// get ready for next slice
		x.clear();
		y.clear();
	}

	return graphs;
}

TH2F * CropCluster (blob b, long frameId, TString type, TString extraname) {

	int wx = b.bP.width_x;
	int wy = b.bP.width_y;

	TString hname = "h_crop_";
	hname += b.bP.geoCenter_x; hname += "_"; hname += b.bP.geoCenter_y;
	hname += "_"; hname += frameId; hname += "_"; hname += extraname;

	TH2F * h = new TH2F(hname, hname, wx, 0, wx, wy, 0, wy);
	// fill non empty spaces

	list< pair < pair<Int_t, Int_t>, Int_t > > des;
	if(type == TString("E")) {
		des = b.GetClusterDescriptionCalibrated();
	} else {
		des = b.GetClusterDescription();
	}

	pair<int, int> corner = FindBottomLeftCorner(b);
	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator i = des.begin();

	int x,y,c;
	for( ; i != des.end() ; i++){

		// translation to the origin
		x = (*i).first.first - corner.first;
		y = (*i).first.second - corner.second;
		c = (*i).second;

		h->Fill(x,y,c);

	}

	return h;
}

pair<int, int> FindBottomLeftCorner (cluster b) {

	pair<int, int> corner;

	float wx = b.bP.width_x;
	float wy = b.bP.width_y;
	float cx = b.bP.geoCenter_x;
	float cy = b.bP.geoCenter_y;

	corner.first  = TMath::Floor( cx - (wx/2.) );
	corner.second = TMath::Floor( cy - (wy/2.) );

	//cout << "Corner : " << corner.first << ", " << corner.second << endl;

	return corner;
}

TGraph2D * ConvertClusterToGraph2D (blob b, int div, long frameId, TString ctype, TString appendname) {

	TString funcname = "TGraph2D_x";
	funcname += (int)TMath::Floor(b.bP.geoCenter_x);
	funcname += "_y";
	funcname += (int)TMath::Floor(b.bP.geoCenter_y);
	funcname += "_frame";
	funcname += frameId;
	return ConvertClusterToGraph2D(b, div, funcname, ctype, appendname);
}

list< pair < pair<int, int>, int > >::iterator GetLowestToA( list< pair < pair<int, int>, int > > c ) {

	list< pair < pair<int, int>, int > >::iterator itr = c.begin();
	list< pair < pair<int, int>, int > >::iterator lowItr = itr;

	int lowest = 0xEFFFFFFF;
	for ( ; itr != c.end() ; itr++ ) {
		if ( (*itr).second < lowest && (*itr).second > 0 ) {
			lowest = (*itr).second;
			lowItr = itr;
		}
	}

	return lowItr;
}
TGraph2D * ConvertClusterToGraph2D(blob b, int div, TString gname, TString ctype, TString appendname){

	// Some definitions
	vector<double> x;
	vector<double> y;
	vector<double> z;
	list< pair < pair<int, int>, int > >::iterator itr;
	list< pair < pair<int, int>, int > > cld;
	if(ctype == "E") {
		cld = b.GetClusterDescriptionCalibrated();
	} else if( ctype == "ToA" ) {
		cld = b.GetClusterDescriptionToA();
	} else {
		cld = b.GetClusterDescription();
	}

	// Check first if this is a good candidate
	if( !GoodCandidateGraph2D(b) ) {
		gname += "_ERROR_1D_OBJECT";
		// In this case foce div = 0.
		// For a 1Dim cluster the bigger mesh
		// is not well defined
		div = 0;
	}

	// append name if requested
	if(appendname.Length() != 0) {
		gname += "_";
		gname += appendname;
	}

	// First I need to copy the cluster in a
	// matrix of size (width_x)*(width_y)

	int i, j;

	int ** c;
	int Nx = b.bP.width_x; // size of the original grid x
	int Ny = b.bP.width_y; // size of the original grid y
	c = new int * [Nx];
	for(int i = 0 ; i < Nx ; i++) {
		c[i] = new int[Ny];
	}
	// Fill with zeros
	for(int i = 0 ; i < Nx ; i++) {
		for(int j = 0 ; j < Ny ; j++) {
			c[i][j] = 0;
		}
	}
	// Find the left-bottom extreme
	int ext_x = TMath::FloorNint(b.bP.geoCenter_x) - (Nx/2);
	int ext_y = TMath::FloorNint(b.bP.geoCenter_y) - (Ny/2);
	list< pair < pair<int, int>, int > >::iterator lowItr;
	if(ctype == "ToA") {
		lowItr = GetLowestToA(cld);
	}
	// Now copy the old cluster.  If the value is not available we will have 0
	for ( itr = cld.begin() ; itr != cld.end() ; itr++ ) {
		// bring it back to i,j indexes in the matrix
		i = (*itr).first.first - ext_x;
		j = (*itr).first.second - ext_y;

		if(ctype == "ToA") {
			c[i][j] = (*itr).second - (*lowItr).second;
			cout << c[i][j] << endl;
		} else {
			c[i][j] = (*itr).second; // TOT or E or ToA of FastToA
		}
		//cout << "i = " << i << ", j = " << j << " : " << c[i][j] << endl;
	}

	// Print the cluster in the original matrix size width_x * width_y with extra zeros
	//PrintMatrix(c, Nx, Ny);

	if(div == 0) {

		// Fill the vectors for TGraph2D from the original copy
		// of the matrix.
		for(int i = 0; i < Nx ; i++) {
			for(int j = 0 ; j < Ny ; j++) {
				x.push_back( i );
				y.push_back( j );
				z.push_back( c[i][j] ); // TOT
				//z.push_back( i*j ); // TOT
			}
		}

	} else {

		int Sx = Nx;
		int Sy = Ny;

		int ** cm = c; // initialize with the old grid
		int ** c_erase; int sox;
		for (int gridItr = 0 ; gridItr < div ; gridItr++) {
			// store point to the old grid to erase after
			c_erase = cm; sox = Sx;
			// initially Sx is the size of the old grid
			cm = RefineGridOnce(cm, &Sx, &Sy);
			// now Sx, Sy is the size of the new grid
			//PrintMatrix(cm, Sx, Sy);
			// Erase old grid
			for(int i = 0 ; i < sox ; i++){delete [] c_erase[i];}
			delete [] c_erase;
		}

		// Fill the vectors for TGraph2D
		for(int i = 0; i < Sx ; i++) {
			for(int j = 0 ; j < Sy ; j++) {
				x.push_back( i );
				y.push_back( j );
				z.push_back( cm[i][j] ); // TOT
			}
		}
		// erase last matrix
		for(int i = 0 ; i < Sx ; i++){delete [] cm[i];}
		delete [] cm;

	}

	TGraph2D * g2 = new TGraph2D((int)x.size(), &x[0], &y[0], &z[0]);
	g2->SetName(gname);


	return g2;
}

TH2I * ConvertClusterToTH2I(cluster cl, long frameId, TString appendname) {

	TString funcname = "TGraph2D_x";
	funcname += (int)TMath::Floor(cl.bP.geoCenter_x);
	funcname += "_y";
	funcname += (int)TMath::Floor(cl.bP.geoCenter_y);
	funcname += "_frame";
	funcname += frameId;
	if(appendname.Length() != 0) {
		funcname += "_";
		funcname += appendname;
	}

	TH2I * h2 = new TH2I(funcname, funcname, cl.bP.width_x, 0, cl.bP.width_x, cl.bP.width_y, 0, cl.bP.width_y);

	list< pair < pair<int, int>, int > > des = cl.GetClusterDescription();
	list< pair < pair<int, int>, int > >::iterator desItr = des.begin();
	int x,y;
	pair<int, int> bl_extreme = GetBottomLeftExtreme(cl);
	for ( ; desItr != des.end() ; desItr++) {
		x = (*desItr).first.first - bl_extreme.first;
		y = (*desItr).first.second - bl_extreme.second;
		h2->Fill(x, y, (*desItr).second); //
	}

	return h2;
}

pair<int, int> GetBottomLeftExtreme(cluster cl) {

	pair<int, int> ext;
	double xwidth = cl.bP.width_x;
	double ywidth = cl.bP.width_y;
	double cx = cl.bP.geoCenter_x;
	double cy = cl.bP.geoCenter_y;

	ext.first = TMath::Floor( cx - xwidth/2. );
	ext.second = TMath::Floor( cy - ywidth/2. );

	return ext;
}

int ** RefineGridOnce(int ** c, int * prev_x, int * prev_y){

	// newNx and newNy bring the previous size of the grid Nx * Ny and are
	// also used to return the new

	// Bring the cluster to a matrix of size (2*width_x - 1)*(2*width_y - 1)
	// prepare the new matrix

	int ** cm;
	int Sx = 2*(*prev_x) - 1; // size of the new grid x // new size !
	int Sy = 2*(*prev_y) - 1; // size of the new grid y

	//cout << "New grid size " << Sx << ", " << Sy << endl;

	cm = new int * [Sx];
	for(int i = 0 ; i < Sx ; i++) {
		cm[i] = new int[Sy];
	}
	// Fill with zeros
	for(int i = 0 ; i < Sx ; i++) {
		for(int j = 0 ; j < Sy ; j++) {
			cm[i][j] = 0;
		}
	}
	// 1) Fill the normal nodes
	for(int i = 0 ; i < Sx ; i+=2) {
		for(int j = 0 ; j < Sy ; j+=2) {
			cm[i][j] = c[i/2][j/2];
		}
	}
	// 2) Fill the inner with only 2 neighbors
	for(int i = 1; i < Sx ; i+=2) {
		for(int j = 0 ; j < Sy ; j+=2) {
			cm[i][j] = (cm[i-1][j] + cm[i+1][j]) / 2;
		}
	}
	for(int i = 0; i < Sx ; i+=2) {
		for(int j = 1 ; j < Sy ; j+=2) {
			cm[i][j] = (cm[i][j-1] + cm[i][j+1]) / 2;
		}
	}
	// 3) Crossed values with 4 neighbors
	for(int i = 1; i < Sx ; i+=2) {
		for(int j = 1 ; j < Sy ; j+=2) {
			cm[i][j] = (cm[i-1][j-1] + cm[i+1][j-1] + cm[i-1][j+1] + cm[i+1][j+1]) / 4;
		}
	}

	// new grid size
	*prev_x = Sx;
	*prev_y = Sy;

	return cm;
}

void PrintMatrix(int ** c, int Nx, int Ny){

	int val;
	int nchars;

	// find max value
	int max = 0;
	for(int i = 0 ; i < Nx ; i++) {
		for(int j = 0 ; j < Ny ; j++) {
			val = c[i][j];
			if(val > max) max = val;
		}
	}
	// according to max value this is the number of
	// caracters to print in each entry to make the
	// matrix look good
	nchars = TMath::FloorNint( TMath::Log10(max) );

	// extract each char and print properly
	char valC[256];
	char uc = 'a';
	int uc_cntr = 0;

	for(int j = Ny-1 ; j >= 0 ; j--) {  // Watch the direction here !
		for(int i = 0 ; i < Nx ; i++) { // trying to match what you observe in the screen :)

			val = c[i][j];
			sprintf(valC,"%d",val);
			uc = 'a'; uc_cntr = 0;
			while(uc != '\0'){
				uc = valC[uc_cntr++];
			}
			uc_cntr--; // <-- number of characters in valC
			// Now print nchars - uc_cntr blank spaces + valC + space
			for(int k = 0 ; k <= nchars-uc_cntr ; k++) { printf(" "); }
			printf("%s",valC);
			printf(" ");
		}
		printf("\n");
	}

}

#define __inner_occ 4

bool PixelIsInner(pair<int,int> pix, blob b) {

	int innercntr = 0;
	int x = pix.first;
	int y = pix.second;
	set<pair<int,int> > cs = b.GetContentSet();
	// test up down right and left
	if( cs.find( make_pair( x+1 , y   ) ) != cs.end() ) innercntr++;
	if( cs.find( make_pair( x-1 , y   ) ) != cs.end() ) innercntr++;
	if( cs.find( make_pair( x   , y+1 ) ) != cs.end() ) innercntr++;
	if( cs.find( make_pair( x   , y-1 ) ) != cs.end() ) innercntr++;

	if(innercntr == __inner_occ) return true;

	return false;
}

vector<pair<int,int> > GetVectorOfInnerPixels(blob b) {

	vector<pair<int,int > > inners;  // inner pixels

	set<pair<int,int> > cs = b.GetContentSet();
	set<pair<int,int> >::iterator itr = cs.begin();

	for( ; itr != cs.end() ; itr++) {
		if( PixelIsInner(*itr, b) ) {
			inners.push_back(*itr);
		}
	}

	return inners;
}

set<pair<int,int> > GetSetOfInnerPixels(blob b) {

	set<pair<int,int > > inners;  // inner pixels

	set<pair<int,int> > cs = b.GetContentSet();
	set<pair<int,int> >::iterator itr = cs.begin();

	for( ; itr != cs.end() ; itr++) {
		if( PixelIsInner(*itr, b) ) {
			inners.insert(*itr);
		}
	}

	return inners;
}

vector<pair<int,int> > GetClusterSkeletton_InnerNegative(blob b) {

	// The skeleton will be cs - inners
	vector<pair<int,int> > inners = GetVectorOfInnerPixels( b );
	vector<pair<int,int> >::iterator itr = inners.begin();
	// Cluster content
	set<pair<int,int> > cs = b.GetContentSet();

	for( ; itr != inners.end() ; itr++){
		cs.erase( *itr ); // erase the inners
	}

	// The remaining in cs is the skeleton
	// Copying it to a std::vector
	vector<pair<int,int> > sk;
	set<pair<int,int> >::iterator itr2 = cs.begin();
	for( ; itr2 != cs.end() ; itr2++){
		sk.push_back( *itr2 );
	}

	return sk;
}

///
/// Same as GetClusterSkeletton_InnerNegative but the skeleton points ought to
/// have a neighboor which is inner
///
vector<pair<int,int> > GetClusterSkeletton_InnerNegativeAndAdjacent(blob b) {

	// The skeleton will be cs - inners
	set<pair<int,int> > inners = GetSetOfInnerPixels( b );
	set<pair<int,int> >::iterator itr = inners.begin();
	// Cluster content
	set<pair<int,int> > cs = b.GetContentSet();

	// Erase the inners
	for( ; itr != inners.end() ; itr++){
		cs.erase( *itr ); // erase the inners
	}
	// Also erase those which are not close to an inner
	set<pair<int,int> >::iterator itrI = cs.begin();
	int x = 0, y = 0;
	int closeToInner = 0;
	vector<pair<int,int > > deltaSubSet;
	for( ; itrI != cs.end() ; itrI++){
		closeToInner = 0;
		x = (*itrI).first;
		y = (*itrI).second;
		if( inners.find( make_pair( x+1 , y   ) ) != inners.end() ) closeToInner++;
		if( inners.find( make_pair( x-1 , y   ) ) != inners.end() ) closeToInner++;
		if( inners.find( make_pair( x   , y+1 ) ) != inners.end() ) closeToInner++;
		if( inners.find( make_pair( x   , y-1 ) ) != inners.end() ) closeToInner++;
		if( inners.find( make_pair( x+1 , y+1 ) ) != inners.end() ) closeToInner++;
		if( inners.find( make_pair( x-1 , y-1 ) ) != inners.end() ) closeToInner++;
		if( inners.find( make_pair( x-1 , y+1 ) ) != inners.end() ) closeToInner++;
		if( inners.find( make_pair( x+1 , y-1 ) ) != inners.end() ) closeToInner++;

		// If there are no inner adjacent to this pixel, schedule to erase
		// !!! This Small subset is interesting because most of the time it corresponds to a delta ray !!!
		if(closeToInner == 0) {
			deltaSubSet.push_back( *itrI );
		}

	}
	// Erase the subset
	vector<pair<int,int > >::iterator itrSub = deltaSubSet.begin();
	for( ; itrSub != deltaSubSet.end() ; itrSub++){
		cs.erase( *itrSub );
	}

	// The remaining in cs is the skeleton, stripped from external legs
	// Copying it to a std::vector
	vector<pair<int,int> > sk;
	set<pair<int,int> >::iterator itr2 = cs.begin();
	for( ; itr2 != cs.end() ; itr2++){
		sk.push_back( *itr2 );
	}

	return sk;
}

list< pair < pair<int, int>, int > >::iterator FindHotestPixel(list< pair < pair<int, int>, int > > * cdes){

	list< pair < pair<int, int>, int > >::iterator i = cdes->begin();
	list< pair < pair<int, int>, int > >::iterator hot;

	int max = 0.;
	for ( ; i != cdes->end() ; i++){
		if ( (*i).second > max ) {
			max = (*i).second;
			hot = i;
		}
	}

	// return the iterator to hottest pixel
	return hot;
}

list< pair < pair<int, int>, int > >::iterator FindColdestPixel(list< pair < pair<int, int>, int > > * cdes, int hotval){

	list< pair < pair<int, int>, int > >::iterator i = cdes->begin();
	list< pair < pair<int, int>, int > >::iterator cold;

	int min = hotval;
	for ( ; i != cdes->end() ; i++) {
		if ( (*i).second < min ) {
			min = (*i).second;
			cold = i;
		}
	}

	// return the iterator to hottest pixel
	return cold;
}

/**
 * Find the hotest pixel and return the TOT value
 */
double GetTOTOfHotestPixel(cluster cl){

	list< pair < pair<int, int>, int > > cd = cl.GetClusterDescription();
	list< pair < pair<int, int>, int > >::iterator i = cd.begin();

	double max = 0.;
	for ( ; i != cd.end() ; i++){
		if ( (*i).second > max ) {
			max = (*i).second;
		}
	}

	// return the value of the hotest pixel
	return max;
}

double GetEnergyOfHotestPixel(cluster cl){

	list< pair < pair<int, int>, int > > cd = cl.GetClusterDescriptionCalibrated();
	list< pair < pair<int, int>, int > >::iterator i = cd.begin();

	double max = 0.;
	for ( ; i != cd.end() ; i++){
		if ( (*i).second > max ) {
			max = (*i).second;
		}
	}

	// return the value of the hotest pixel
	return max;
}


set<pair<int, int> > FindPixelsAboveAverageTOT(cluster cl) {

	set<pair<int, int> > hotpixels;

	list< pair < pair<int, int>, int > > cdes = cl.GetClusterDescription();
	list< pair < pair<int, int>, int > >::iterator i = cdes.begin();

	double averageTOT = (double) cl.bP.clusterTOT / (double) cl.bP.clusterSize;

	for ( ; i != cdes.end() ; i++ ) {
		if( (double) ( (*i).second ) > averageTOT ) {
			hotpixels.insert( (*i).first );
		}
	}

	return hotpixels;
}

set<pair<int, int> > FindPixelsBelowOrEqualAverageTOT(cluster cl) {

	set<pair<int, int> > coldpixels;

	list< pair < pair<int, int>, int > > cdes = cl.GetClusterDescription();
	list< pair < pair<int, int>, int > >::iterator i = cdes.begin();

	double averageTOT = (double) cl.bP.clusterTOT / (double) cl.bP.clusterSize;

	for ( ; i != cdes.end() ; i++ ) {
		if( (double) ( (*i).second ) <= averageTOT ) {
			coldpixels.insert( (*i).first );
		}
	}

	return coldpixels;
}

set<pair<int, int> > FindPixelsBelow(cluster cl, int th) {

	set<pair<int, int> > coldpixels;

	list< pair < pair<int, int>, int > > cdes = cl.GetClusterDescription();
	list< pair < pair<int, int>, int > >::iterator i = cdes.begin();

	for ( ; i != cdes.end() ; i++ ) {
		if( (double) ( (*i).second ) <= th ) {
			coldpixels.insert( (*i).first );
		}
	}

	return coldpixels;
}


/**
 *   Mask (set to 0) the pixels coming in the vector "mask"
 */
int MaskInMap(map<int, int> & totmap, vector<pair<int, int> > mask, int dimX){

	int nMasked = 0;
	vector<pair<int, int> >::iterator i = mask.begin();

	for ( ; i != mask.end() ; i++ ) {
		totmap [ XYtoC(*i, dimX) ] = 0; // masking
		nMasked++;
	}

	return nMasked;
}

/**
 *   Mask (set to 0) the pixels coming in the vector "mask"
 */
int MaskInMap(map<int, int> & totmap, set<pair<int, int> > mask, int dimX){

	int nMasked = 0;
	set<pair<int, int> >::iterator i = mask.begin();

	for ( ; i != mask.end() ; i++ ) {
		totmap [ XYtoC(*i, dimX) ] = 0; // masking
		nMasked++;
	}

	return nMasked;
}


/**
 *   Mask (set to 0) the pixels coming in the vector "mask"
 */
int MaskInCluster(cluster & cl, vector<pair<int, int> > mask, int dimX){

	int nMasked = 0;
	vector<pair<int, int> >::iterator i = mask.begin();

	for ( ; i != mask.end() ; i++ ) {

		cl.DeletePixelEntry( *i , dimX, dimX );
		nMasked++;
	}

	return nMasked;
}

vector<pair<double, double> > FindLimitPixels_MaxMinXY(blob theBlob) {

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


bool ClusterTouchesBorder(cluster cl, int dimx, int dimy, int borderExclusion){

	set<pair<int, int> > cs = cl.GetContentSet();
	set<pair<int, int> >::iterator i = cs.begin();
	int x,y;

	for ( ; i != cs.end() ; i++ ) {
		x = (*i).first;
		y = (*i).second;
		if( x >= dimx - borderExclusion ) return true;
		if( x <= borderExclusion)         return true;
		if( y >= dimy - borderExclusion ) return true;
		if( y <= borderExclusion )        return true;
	}

	return false;
}

set< pair<int, int> > VectorToSet(vector<pair<int, int> > v) {

	set< pair<int, int> > s;
	vector<pair<int, int> >::iterator i = v.begin();
	for( ; i != v.end() ; i++ ) {
		s.insert( *i );
	}

	return s;
}

bool PixelInListIsAClusterConstituent(blob cl, int * mask , int sizex, int sizey) {

	list< pair < pair<int, int>, int > > cdes = cl.GetClusterDescription();

	list< pair < pair<int, int>, int > >::iterator i = cdes.begin();

	int occ = sizex*sizey;
	pair<int, int> pix;
	int j;
	for ( ; i != cdes.end() ; i++ ) {
		pix = (*i).first;
		if( mask[ XYtoX( pix, sizex ) ] != 0) { // constituent found in the mask
			return true;
		}
	}

	return false;
}

bool PixelInListIsAClusterConstituent(blob cl, vector<int>  mask , int sizex, int sizey) {

	list< pair < pair<int, int>, int > > cdes = cl.GetClusterDescription();

	list< pair < pair<int, int>, int > >::iterator i = cdes.begin();

	int occ = sizex*sizey;
	pair<int, int> pix;
	int j;
	for ( ; i != cdes.end() ; i++ ) {
		pix = (*i).first;
		if( mask[ XYtoX( pix, sizex ) ] != 0) { // constituent found in the mask
			return true;
		}
	}

	return false;
}


double GausFuncAdd(double * x, double * par) {

	double xx = x[0];

	// very first parameter is the nbins which is the
	// same range of the kernel density function
	int N = par[0];

	// C-array "par" contains 3N + 1 parameters
	// N = par[               0 ] : Number of entries per parameter
	// a = par[       1 -->   N ] : N constants
	// b = par[   N + 1 --> 2*N ] : N mean vaues
	// c = par[ 2*N + 1 --> 3*N ] : N sigma values

	// sum func
	double func = 0;
	int Nd = 2*N;
	for ( int i = 0 ; i < N ; i++ ) {
		func += par[ i + 1 ] * exp( - ( ( xx - par[ N + i + 1 ] )*( xx - par[ N + i + 1 ] ) ) / ( 2.*par[ Nd + i + 1 ]*par[ Nd + i + 1 ] ) );
	}

	return func;
}

TF1 * CreateKernelDensityFunction(int id, vector<double> hist, double bandwidth) {

	int nbins = (int) hist.size();

	// Set of parameters for the full kernel density function
	// Npars = m_nbins*( constants + sigmas + means ) + size
	// ex: for 100 bins we need 301 parameters
	double * par = new double[nbins * 3 + 1];
	par[0] = nbins;
	//double distmax = *max_element(hist.begin(), hist.end());

	int i = 1;
	for ( ; i <= nbins ; i++ ) {

		par[i] = hist[i-1]; // constant

		par[i + nbins] = i; //i-1; // mean

		par[i + 2*nbins] = bandwidth; // sigma

	}
	//cout << "-------" << endl;

	// Final kernel density
	TString kernelname = "kernel_";
	//kernelname += m_calhandler->GetSourcename();
	kernelname += "_";
	kernelname += id;
	//cout << "Creating kernel function : " << kernelname << endl;
	TF1 * fker = new TF1(kernelname, GausFuncAdd, 1, nbins, (3 * nbins) + 1);
	fker->SetParameters( par );

	return fker;
}

double DerivativeFivePointsStencil(TF1 * f, double x, double h) {

	double der = 0.;

	der -=      f->Eval( x + 2*h );
	der += 8. * f->Eval( x +   h );
	der -= 8. * f->Eval( x -   h );
	der +=      f->Eval( x - 2*h );
	der /= 12.*h;

	return der;
}

TGraph * GetDerivative(TGraph * g) {

	int N = g->GetN();
	double * x = g->GetX();
	double * y = g->GetY();

	queue<double> qd;
	vector<double> der;
	vector<double> xp;

	//for( int i = 0 ; i < N ; i++ ) {
	for( int i = N-1 ; i >= 0 ; i-- ) {

		//if( i < 5 ) {
		if( i > N-6 ) {

			qd.push( y[i] );

		} else {
			//der.push_back(  DerivativeFivePointsStencil( qd, x[i-2] - x[i-3] )  );
			der.push_back(  DerivativeFivePointsStencil( qd, x[i+3] - x[i+2] )  );

			//xp.push_back( x[i-3] );
			xp.push_back( x[i+3] );

			// remove oldest and push a new one
			qd.pop();
			qd.push( y[i] );
		}

	}

	TGraph * gder = new TGraph( (int)xp.size(), &xp[0], &der[0] );

	return gder;
}

// Applying triangular smoothing.
// http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html

TGraphErrors * TriangularSmooth(TGraphErrors * g) {

	int ns = g->GetN();
	if(ns < 5) return 0x0; // at least 5 points are needed

	double * x = g->GetX();
	double * y = g->GetY();

	vector<double> x_smooth;
	vector<double> y_smooth;
	vector<double> x_err_smooth;
	vector<double> y_err_smooth;

	for ( int j = 2 ; j < ns - 2 ; j++ ) {

		x_smooth.push_back( x[j] );
		y_smooth.push_back( ( y[j-2] + 2*y[j-1] + 3*y[j] + 2*y[j+1] + y[j+2] ) / 9. );
		x_err_smooth.push_back( 0. );
		y_err_smooth.push_back( 0. );

	}

	TGraphErrors * gs = new TGraphErrors(
			(int)x_smooth.size(),
			&x_smooth[0],
			&y_smooth[0],
			&x_err_smooth[0],
			&y_err_smooth[0]
	);

	return gs;
}

/**
 * Rough 5-points stencil derivative on a queue
 */
double DerivativeFivePointsStencil(queue<double> q, double h) {

	// The queue needs to have at least 5 points.
	if(q.size() < 5) return 0.;

	double der = 0.;

	// Calculate the derivative considering only the 5 "oldest" entries
	//double h = 1;

	der -=      q.front(); // reference to oldest element
	q.pop();               // and pop out

	der += 8. * q.front();
	q.pop();

	q.pop();

	der -= 8. * q.front();
	q.pop();

	der +=      q.front();
	q.pop();

	der /= 12.*h;

	return der;
}

enum {
	__s_pos = 0,
	__s_neg,
};

int GetCriticalPoints(TF1 * f, int nbins, vector<double> & min, vector<double> & max) {

	int ncrit = 0;
	short sign = __s_pos;

	double der = 0.;
	double step = 1.;
	for (double x = 1. ; x <= (double)nbins ; x+= step) {

		//der = f->Derivative(x, 0x0, 0.1);
		der = DerivativeFivePointsStencil(f, x, 0.1);

		if(der < 0 && sign == __s_pos) { // Flip from pos to neg --> maximum
			//max.push_back( x - step );  // Critical point in the previous step
			max.push_back( x ); // Critical point here
			ncrit++;
		}
		if(der > 0 && sign == __s_neg) { // Flip from neg to pos --> minimum
			//min.push_back( x - step );  // Critical point in the previous step
			min.push_back( x );  // Critical point here
			ncrit++;
		}

		if(der > 0) sign = __s_pos;
		if(der < 0) sign = __s_neg;

		//cout << " x = " << x << " --> der = " << der << endl;
	}

	return ncrit;
}

/**
 *  Follows a spiral shape around a pixel
 */
void RewindSpiral(steerSpiral & spiral){

	spiral.xySwitch = X_STEP;
	spiral.localMax = 1;
	spiral.xCntr = 1;
	spiral.yCntr = 1;
	spiral.dir = 1;

}

pair<int, int> GetNextPosition(pair<int, int> current, steerSpiral & spiral) {

	pair<Int_t, Int_t> newPos;

	if(spiral.xySwitch == X_STEP) {

		newPos.first = current.first - spiral.dir;
		newPos.second = current.second;
		spiral.xCntr++;

		if(spiral.xCntr > spiral.localMax)
			spiral.xySwitch = Y_STEP;

	} else if(spiral.xySwitch == Y_STEP) {

		newPos.second = current.second - spiral.dir;
		newPos.first = current.first;
		spiral.yCntr++;

		if(spiral.yCntr > spiral.localMax)
		{
			spiral.xySwitch = X_STEP;
			spiral.localMax++;
			spiral.xCntr = 1;
			spiral.yCntr = 1;
			spiral.dir *= -1; // change direction
		}

	}

	return newPos;
}

int CheckDiscontinuityWithModifiedMatrix(blob bl, map<int,int> totmap, int discontinuityTolerance, int width, int height, MediPixAlgo * algo){

	// New object containing all the blobs
	AllBlobsContainer * reclusteringContainer = new AllBlobsContainer(algo);
	//Log << MSG::DEBUG << "Re-clustering container --> " << reclusteringContainer << endreq;

	MAFTools::DigitalBlobReclustering(discontinuityTolerance, reclusteringContainer, bl, width, height, totmap);

	//Log << MSG::DEBUG << "Re-clustering ... found " << reclusteringContainer->GetBlobsVector().size() << " clusters in the blob" << endreq;

	int nSubClusters = reclusteringContainer->GetBlobsVector().size();
	delete reclusteringContainer;

	return nSubClusters;
}

int FindSubclusters (blob b, map<int, int> totmap, double multipeak_divider, int width, int height, MediPixAlgo * algo) {

	// Find the pixels below the average threshold, i.e. cold pixels
	double meanTOT = (double) b.bP.clusterTOT / (double) b.bP.clusterSize;
	double hottest = MAFTools::GetTOTOfHotestPixel( b );
	// cut-off half way between the mean TOT and the max
	double cutoff_cold = ( meanTOT + hottest ) / multipeak_divider;
	set<pair<int, int> > coldpixels = MAFTools::FindPixelsBelow(b, cutoff_cold);

	// Make a copy of the cluster but masking the cold pixels
	cluster b2(b, coldpixels);
	//b2.AllDump();

	// Re-cluster and see if the hot spot is fragmented --> If that is the case this is probably an overlap
	//map<int, int> totmap = GetTOTMap();
	MAFTools::MaskInMap(totmap, coldpixels, width);
	int subclusters = CheckDiscontinuityWithModifiedMatrix(b2, totmap, 0, width, height, algo);
	//Log << MSG::DEBUG << "number of sub-clusters : " << subclusters << endreq;

	return subclusters;
}


double CalcStdDev(vector<double> v) {

	double sd = 0.;
	double mean = CalcMean(v);
	if( v.size() == 0 ) return sd;
	double N = (double) v.size();
	vector<double>::iterator i;

	// Calculate Sample Standard Deviation
	double fac;
	for ( i = v.begin() ; i != v.end() ; i++ ) {
		fac = (*i) - mean;
		sd += fac*fac;
	}
	sd /= (N - 1.);
	sd = TMath::Sqrt( sd );

	return sd;
};

double CalcStdDev(vector<int> v) {

	double sd = 0.;
	double mean = CalcMean(v);
	if( v.size() == 0 ) return sd;
	double N = (double) v.size();
	vector<int>::iterator i;

	// Calculate Sample Standard Deviation
	double fac;
	for ( i = v.begin() ; i != v.end() ; i++ ) {
		fac = (double)(*i) - mean ;
		sd += fac*fac;
	}
	sd /= (N - 1.);
	sd = TMath::Sqrt( sd );

	return sd;
};

double CalcMean(vector<int> v) {

	double mean = 0.;
	int N = (int) v.size();
	vector<int>::iterator i;

	// Calculate mean
	for ( i = v.begin() ; i != v.end() ; i++ ) {
		mean += *i;
	}
	mean /= (double)N;

	return mean;
}

double CalcMean(vector<double> v) {

	double mean = 0.;
	int N = (int) v.size();
	vector<double>::iterator i;

	// Calculate mean
	for ( i = v.begin() ; i != v.end() ; i++ ) {
		mean += *i;
	}
	mean /= (double) N;

	return mean;
}

void DumpVector(vector<double> v ) {

	vector<double>::iterator i;

	// Calculate mean
	cout << "< ";
	for ( i = v.begin() ; i != v.end() ; i++ ) {
		cout << *i;
		if( i != v.end() ) cout << ", ";
		else cout << " >";
	}
	cout << endl;
}

void DumpVector(vector<int> v ) {

	vector<int>::iterator i;

	// Calculate mean
	cout << "< ";
	int indx = 0;
	for ( i = v.begin() ; i != v.end() ; i++ ) {
		cout << "[" << indx++ << "]" << *i;
		if( i != v.end() ) cout << ", ";
		else cout << " >";
	}
	cout << endl;
}

TH1F * GetHistogramFromTGraph(TGraph * g, data_direction dd ){

	TString hname = g->GetName();
	hname += "_extractHisto";
	double * x = g->GetX();
	double * y = g->GetY();
	int N = g->GetN();

	// Conver x to a vector just o identify min and max
	vector<double> xv;
	for(int i = 0 ; i < N ; i++){
		xv.push_back( x[i] );
	}
	int min = *( min_element( xv.begin(), xv.end()) );
	int max = *( max_element( xv.begin(), xv.end()) );

	cout << "min = " << min << " | max = " << max << endl;
	TH1F * gh = new TH1F(hname, hname, N, min, max);

	if ( dd == __righthanded ) {
		for(int i = 0 ; i < N ; i++) {
			gh->SetBinContent(i+1, y[i]);
		}
	} else {
		int cntr = 1;
		for(int i = N ; i > 0 ; i--) { // Reverse, still in bins !
			gh->SetBinContent(cntr, y[i-1]);
			cntr++;
		}
	}
	return gh;
}

void FindMinMaxForFit(TH1 * h, double mean, double & minf, double & maxf) {

	// start point
	minf = mean - 1;
	maxf = mean + 1;

	int centerBin = h->FindBin(mean);
	double heightAtKernelHint = h->GetBinContent(centerBin);

	unsigned int nBins = h->GetXaxis()->GetNbins();
	unsigned int rightBin = h->FindBin( maxf );
	unsigned int leftBin = h->FindBin( minf );

	int underThresholdCntr = 0;
	for(unsigned int i=centerBin; i<=nBins; ++i ) {
		if( h->GetBinContent(i) < heightAtKernelHint * __fraction_of_height_range_id ) {
			underThresholdCntr++;
			rightBin = i;
			if(underThresholdCntr > 1) break;
		}
	}
	underThresholdCntr = 0;
	for(unsigned int i=centerBin; i>1; --i ) {
		if( h->GetBinContent(i) < heightAtKernelHint * __fraction_of_height_range_id ) {
			underThresholdCntr++;
			leftBin = i;
			if(underThresholdCntr > 1) break;
		}
	}
	// take one more
	if(rightBin < nBins) rightBin++;
	if(leftBin > 1) leftBin--;

	// define maxf and minf
	maxf = h->GetBinCenter(rightBin);
	minf = h->GetBinCenter(leftBin);

}

/*
template<typename T> double CalcStdDev(vector<T> v) {

	double sd = 0., mean = 0.;
	if( v.size() == 0 ) return sd;
	double N = (double) v.size();

	vector<T>::iterator i;
	// calc mean for Sample Standard Deviation
	for ( i = v.begin() ; i != v.end() ; i++ ) {
		mean += (double)(*i);
	}
	mean /= (N - 1.);

	// calc Sample Standard Deviation
	double fac;
	for ( i = v.begin() ; i != v.end() ; i++ ) {
		fac = ( (double)(*i) - mean );
		sd += fac*fac;
	}
	sd /= (N - 1.);
	sd = TMath::Sqrt( sd );

	return sd;
};
 */

}
#endif
