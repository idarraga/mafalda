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

#ifndef Highlighter_h
#define Highlighter_h

#include <vector>

#include <TROOT.h>
#include <TArrow.h>
#include <TEllipse.h>
#include <TBrowser.h>
#include <TEllipse.h>
#include <TString.h>
#include <TMath.h>
#include <TLine.h>

#include "MPXAlgo/MediPixAlgo.h"
#include "MPXStoreGate/MPXStoreGate.h"
#include "MPXStoreGate/CandidateContainer.h"

#include "OutputMng.h"

#define ARROW_LENGTH     20
#define ARROW_LENGTH_X   20
#define ARROW_LENGTH_Y   20
#define ELLIPSE_R1       10
#define ELLIPSE_R2       10
#define ELLIPSE_SOLID1_R1   5
#define ELLIPSE_SOLID1_R2   5
#define ARROW           100
#define ELLIPSE         200
#define LINE            300
#define USER_DEF_ELLIPSE 400
#define USER_DEF_CIRCLE  500

#define X1 0
#define Y1 1
#define X2 2
#define Y2 3

using namespace std;

/**
 * Class for drawing highlighters in the Viewer.
 * I will may be move this somewhere else later. 
 *
 */

class MAFTArrow;
class MAFTEllipse;
class MAFTLine;

class Highlighter : public CandidateContainer {

public:
	Highlighter();
	Highlighter(Double_t, Double_t, TString, MediPixAlgo *);
	Highlighter(Double_t, Double_t, Double_t, MediPixAlgo *);
	Highlighter(Double_t, Double_t, Double_t, Double_t, MediPixAlgo *);
	Highlighter(pair<Double_t, Double_t>, TString, MediPixAlgo *);

	// coordinates completely defined by the user
	Highlighter(Double_t, Double_t, Double_t, Double_t, TString, MediPixAlgo *);

	~Highlighter(){};
	OutputMng Log;
	MSG::Endreq endreq;

	Int_t GetType(){return m_type;};
	TObject * GetDrawable();

#ifndef __ROOTCINT__ // hiding this from rootcint
	template <class T> void UploadDisplayValue(const Char_t *, T);
	void UploadDisplayValue(const Char_t *, vector<Int_t>);
#endif

	void SetLineColor(Color_t);
	void SetLineWidth(Int_t);
	void SetLineStyle(Int_t);
	void DumpDisplayValues();

private:

	void PrepareArrowCoordinates(Double_t *);
	// vars
	MAFTArrow * m_arrow;
	MAFTEllipse * m_ellipse;
	MAFTLine * m_line;

	Int_t m_type;

	// values to display through
	map<Int_t, TString>  m_characteristicsMap;
	Int_t m_characteristicsCntr;

	ClassDef(Highlighter, 0)
};

/**
 * A MAFalda wrap for TArrow
 *
 */
class MAFTArrow : public TArrow {

public:
	MAFTArrow();
	MAFTArrow(Highlighter *, Double_t, Double_t,
			Double_t, Double_t, Float_t, Option_t *);
	~MAFTArrow(){};
	void DumpMAFaldaValues(); // *MENU*

private:
	Highlighter * m_HighligterOwner;

	ClassDef(MAFTArrow, 0)
};

/**
 * A MAFalda wrap for TEllipse
 *
 */
class MAFTEllipse : public TEllipse {

public:
	MAFTEllipse();
	MAFTEllipse(Highlighter *, Double_t, Double_t, Double_t, Double_t);
	~MAFTEllipse(){};
	void DumpMAFaldaValues(); // *MENU*

private:
	Highlighter * m_HighligterOwner;

	ClassDef(MAFTEllipse, 0)
};

/**
 * A MAFalda wrap for TLine
 *
 */

class MAFTLine : public TLine {

public:
	MAFTLine();
	MAFTLine(Highlighter *, Double_t, Double_t, Double_t, Double_t);
	~MAFTLine(){};
	void DumpMAFaldaValues(); // *MENU*

private:
	Highlighter * m_HighligterOwner;

	ClassDef(MAFTLine, 0)
};

#endif
