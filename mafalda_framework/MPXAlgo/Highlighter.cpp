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

#ifndef Highlighter_cxx
#define Highlighter_cxx

#include <typeinfo>

#include "Highlighter.h"

ClassImp(Highlighter)
ClassImp(MAFTArrow)
ClassImp(MAFTEllipse)
ClassImp(MAFTLine)

/**
 * Class for drawing highlighters in the Viewer.
 * I will may be move this somewhere else later. 
 *
 */

Highlighter::Highlighter() : 
CandidateContainer () {

}

// Exclusive for user defined ellipse
Highlighter::Highlighter(Double_t x, Double_t y, Double_t a, Double_t b, MediPixAlgo * algo) :
						  CandidateContainer(algo, MPXDefs::VIS) {

	Log.setAlgoName("Highlighter");

	m_arrow = 0;
	m_ellipse = 0;
	m_line = 0;

	// make an marker
	m_ellipse = new MAFTEllipse(this, x, y, a, b);
	m_ellipse->SetLineWidth(2);
	m_ellipse->SetFillStyle(1);
	m_type = ELLIPSE;

	m_characteristicsCntr = 0;
}

// Exclusive for user defined circle
Highlighter::Highlighter(Double_t x, Double_t y, Double_t r, MediPixAlgo * algo) :
						  CandidateContainer(algo, MPXDefs::VIS) {

	Log.setAlgoName("Highlighter");

	m_arrow = 0;
	m_ellipse = 0;
	m_line = 0;

	// make an marker
	m_ellipse = new MAFTEllipse(this, x, y, r, r);
	m_ellipse->SetLineWidth(2);
	m_ellipse->SetFillStyle(1);
	m_type = ELLIPSE;

	m_characteristicsCntr = 0;
}


Highlighter::Highlighter(Double_t x, Double_t y, TString opt, MediPixAlgo * algo) :
						  CandidateContainer(algo, MPXDefs::VIS) {

	Log.setAlgoName("Highlighter");

	m_arrow = 0;
	m_ellipse = 0;
	m_line = 0;

	// make an marker
	if(opt == "circle")
	{
		m_ellipse = new MAFTEllipse(this, x, y, ELLIPSE_R1, ELLIPSE_R2);
		m_ellipse->SetLineWidth(2);
		m_ellipse->SetFillStyle(1);
		m_type = ELLIPSE;
	}
	else if(opt == "solidcircle")
	{
		m_ellipse = new MAFTEllipse(this, x, y, ELLIPSE_SOLID1_R1, ELLIPSE_SOLID1_R2);
		m_ellipse->SetLineWidth(2);
		m_ellipse->SetFillStyle(1);
		m_type = ELLIPSE;
	}
	else if (opt == "arrow")
	{
		Double_t arrowcoords[] = {x, y, 0, 0};
		PrepareArrowCoordinates(arrowcoords);
		m_arrow = new MAFTArrow(this, *arrowcoords, *(arrowcoords+Y1),
				*(arrowcoords+X2), *(arrowcoords+Y2),
				0.02, "<");
		m_arrow->SetLineWidth(2);
		m_type = ARROW;
	}
	else if(opt == "line")
	{
		Double_t arrowcoords[] = {x, y, 0, 0};
		PrepareArrowCoordinates(arrowcoords);
		m_line = new MAFTLine(this, *arrowcoords, *(arrowcoords+Y1),
				*(arrowcoords+X2), *(arrowcoords+Y2));
		m_line->SetLineWidth(2);
		m_type = LINE;
	}
	else {
		Log << MSG::WARNING << "When defining a marker, it ought to be one of the following :" << endreq;
		Log << MSG::WARNING << "circle, arrow or line" << endreq;
		Log << MSG::WARNING << "Building and arrow by default !" << endreq;

		Double_t arrowcoords[] = {x, y, 0, 0};
		PrepareArrowCoordinates(arrowcoords);
		m_arrow = new MAFTArrow(this, *arrowcoords, *(arrowcoords+Y1),
				*(arrowcoords+X2), *(arrowcoords+Y2),
				0.02, "<");
		m_arrow->SetLineWidth(2);
		m_type = ARROW;

	}

	m_characteristicsCntr = 0;

}

Highlighter::Highlighter(std::pair<Double_t, Double_t> xy, TString opt, MediPixAlgo * algo) :
						  CandidateContainer(algo, MPXDefs::VIS) {

	Highlighter(xy.first, xy.second, opt, algo);

}

Highlighter::Highlighter(Double_t x, Double_t y, Double_t x2, Double_t y2, TString, MediPixAlgo * algo) :
						  CandidateContainer(algo, MPXDefs::VIS) {

	Log.setAlgoName("Highlighter");

	m_arrow = 0;
	m_ellipse = 0;
	m_line = 0;

	m_line = new MAFTLine(this, x, y, x2, y2);
	m_line->SetLineWidth(2);
	m_type = LINE;


	m_characteristicsCntr = 0;
}

void Highlighter::SetLineColor(Color_t color){
	if(m_arrow) m_arrow->SetLineColor(color);
	if(m_ellipse) m_ellipse->SetLineColor(color);
	if(m_line) m_line->SetLineColor(color);
}

void Highlighter::SetLineWidth(Int_t wd){
	if(m_arrow) m_arrow->SetLineWidth(wd);
	if(m_ellipse) m_ellipse->SetLineWidth(wd);
	if(m_line) m_line->SetLineWidth(wd);
}

void Highlighter::SetLineStyle(Int_t wd){
	if(m_arrow) m_arrow->SetLineStyle(wd);
	if(m_ellipse) m_ellipse->SetLineStyle(wd);
	if(m_line) m_line->SetLineStyle(wd);
}

void Highlighter::PrepareArrowCoordinates(Double_t * coords){

	Double_t * x1 = coords;
	Double_t * y1 = coords+Y1;
	Double_t * x2 = coords+X2;
	Double_t * y2 = coords+Y2;

	Double_t xS = GetFrameXSize();
	Double_t yS = GetFrameYSize();

	// North or south hemisphere, calculate the position of the arrow.  Something that
	//  looks easy in the screen.  This equation produces arrows pointing radially inwards.
	if(*y1 <= yS/2 && *y1 > 0) // south hemisphere
	{
		*x2 = TMath::FloorNint( ARROW_LENGTH*
				TMath::Sqrt(
						1.0/(
								1.0 + TMath::Tan( TMath::Pi()/2.0 - (TMath::Pi()*(Double_t)(*y1)/(Double_t)yS) )
						)

				) );
		*y2 = TMath::Sqrt(ARROW_LENGTH*ARROW_LENGTH - (Double_t)(*x2)*(Double_t)(*x2));
		*y2 = (*y1 - *y2);
	}
	else if(*y1 > yS/2 && *y1 < yS-1)  // north hemisphere
	{
		*x2 = TMath::FloorNint( ARROW_LENGTH*
				TMath::Sqrt(
						1.0/(
								1.0 + TMath::Tan( (TMath::Pi()*(Double_t)(*y1)/(Double_t)yS) - TMath::Pi()/2.0 )
						)
				) );
		*y2 = TMath::Sqrt(ARROW_LENGTH*ARROW_LENGTH - (Double_t)(*x2)*(Double_t)(*x2));
		*y2 = (*y1 + *y2);
	}
	else if(*y1 == yS-1) // north bound
	{
		*x2 = 0;
		*y2 = *y1 + ARROW_LENGTH;
	}
	else if(*y1 == 0)  // south bound
	{
		*x2 = 0;
		*y2 = *y1 - ARROW_LENGTH;
	}
	else
	{} // no other cases.

	// west or east emisphere
	if(*x1 <= xS/2)
	{
		*x2 = *x1 - (*x2);
	}
	else
	{
		*x2 = *x1 + (*x2);
	}

}

TObject * Highlighter::GetDrawable(){

	if(m_type == ELLIPSE){
		return m_ellipse;
	}else if(m_type == ARROW){
		return m_arrow;
	}else if(m_type == LINE){
		return m_line;
	}

	// default
	return m_arrow;
}

//template function void Highlighter::UploadDisplayValue(Int_t);

void Highlighter::DumpDisplayValues(){

	if(m_characteristicsMap.empty()){
		Log << MSG::ALWAYS << "No characteristics loaded" << endreq;
		return;
	}

	map<Int_t, TString>::iterator mapItr = m_characteristicsMap.begin();

	for( ; mapItr != m_characteristicsMap.end() ; mapItr++){
		Log << MSG::ALWAYS << (*mapItr).second << endreq;
	}
	Log << MSG::ALWAYS << endreq;

}

/**
 * Upload a value to the highlighter which can be dump
 *   by the DumpMAFaldaValues() member registered in the menu.
 * 
 */
template <class T> void Highlighter::UploadDisplayValue(const Char_t * info, T val){

	double d1;
	float f1;
	int i1;

	m_characteristicsMap[m_characteristicsCntr].Append( TString(info) );
	m_characteristicsMap[m_characteristicsCntr] += " : ";
	TString valS;
	if(typeid(val).name() == typeid(i1).name() ) {
		valS += val;
	} else if(typeid(val).name() == typeid(d1).name() || typeid(val).name() == typeid(f1).name()){
		valS = TString::Format("%.2f", val);
	}
	m_characteristicsMap[m_characteristicsCntr] += valS;
	m_characteristicsCntr++;

}

// In a shared lib scenario (or if you define the template in a different file, like here) 
//  you have to tell what are you going to instantiate
template void Highlighter::UploadDisplayValue(const Char_t *, Int_t);
template void Highlighter::UploadDisplayValue(const Char_t *, Float_t);
template void Highlighter::UploadDisplayValue(const Char_t *, Double_t);

// Case with vectors
void Highlighter::UploadDisplayValue(const Char_t * info, vector<Int_t> val){

	m_characteristicsMap[m_characteristicsCntr].Append( TString(info) );
	m_characteristicsMap[m_characteristicsCntr] += " : ";

	vector<Int_t>::iterator itr = val.begin();
	for( ; itr != val.end() ; itr++){
		m_characteristicsMap[m_characteristicsCntr] += (*itr);
		if(itr+1 != val.end()) m_characteristicsMap[m_characteristicsCntr] += ", ";
	}

	m_characteristicsCntr++;

}
//template void Highlighter::UploadDisplayValue(const Char_t *, vector<Int_t>);


/** 
 *  MAFTArrow
 *
 *
 */

MAFTArrow::MAFTArrow() : 
						  TArrow() {

}

MAFTArrow::MAFTArrow(Highlighter * hl, Double_t x1, Double_t y1, Double_t x2, Double_t y2, 
		Float_t arrowsize = 0.05, Option_t* option = ">") :
		TArrow(x1, y1, x2, y2, arrowsize, option) {

	m_HighligterOwner = hl;

}

void MAFTArrow::DumpMAFaldaValues() {

	// display values are in the Highlighter
	m_HighligterOwner->DumpDisplayValues();

}

/** 
 *  MAFTEllipse
 *
 *
 */

MAFTEllipse::MAFTEllipse() : 
						  TEllipse() {

}

MAFTEllipse::MAFTEllipse(Highlighter * hl, Double_t x1, Double_t y1, Double_t r1, Double_t r2) : 
						  TEllipse(x1, y1, r1, r2) {

	m_HighligterOwner = hl;

}

void MAFTEllipse::DumpMAFaldaValues() {

	// display values are in the Highlighter
	m_HighligterOwner->DumpDisplayValues();

}

/** 
 *  MAFTLine
 *
 *
 */

MAFTLine::MAFTLine() : 
						  TLine() {

}

MAFTLine::MAFTLine(Highlighter * hl, Double_t x1, Double_t y1, Double_t r1, Double_t r2) : 
						  TLine(x1, y1, r1, r2) {

	m_HighligterOwner = hl;

}

void MAFTLine::DumpMAFaldaValues() {

	// display values are in the Highlighter
	m_HighligterOwner->DumpDisplayValues();

}


#endif
