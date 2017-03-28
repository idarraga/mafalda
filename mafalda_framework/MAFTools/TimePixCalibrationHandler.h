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

#ifndef TIMEPIXCALIBRATIONHANDLER_H_
#define TIMEPIXCALIBRATIONHANDLER_H_

#ifndef __MAFALDA_EMBEDDED
#include <TF1.h>
#endif

#include <vector>
#include <set>
#include <map>

#include "MPXAlgo/OutputMng.h"

using namespace std;

namespace MAFTools {

///////////////////////////////////////////////
// TimePix specific
// Handler of the calibration Matrix

class TimePixCalibrationHandler {

public:
	// Provide the 4 calibration files, a,b,c,t
	TimePixCalibrationHandler(){};
	TimePixCalibrationHandler(const char *, const char *, const char *, const char *, int, int);
	~TimePixCalibrationHandler();

	double SurrogateCalibFunction(double x, double *);
	int GetTOT(std::pair<int, int>, double);
	void GetCuadraticSolutions (pair<int,int> pix, int tot, double & sol1, double & sol2, double & sol);
	void GetSolutionHighEnergy (pair<int,int> pix, int tot, double & sol1, double & sol2, double & sol);
	double GetE(std::pair<int, int>, int);
	int XYtoC(std::pair<int, int>, int);
	int XYtoC(int x, int y, int);

	void WordWithIncompleteMatrix(){ m_workWithIncompleteMatrix = true; };

	// log service
	void SetOutputLevel(MSG::Level l){Log.OutputLevel = l;};
	OutputMng Log;
	MSG::Endreq endreq;

#ifndef __MAFALDA_EMBEDDED
	TF1 * GetSurrogateTF1(std::pair<int, int>);
#endif

private:

	double * m_a;
	double * m_b;
	double * m_c;
	double * m_t;

	int m_width;
	int m_height;

	bool m_workWithIncompleteMatrix;

	int m_badCalibCntr;

#ifndef __MAFALDA_EMBEDDED
	std::map<std::pair<int, int>, TF1* > m_surrogateMap;
#endif
	std::set<std::pair<int, int> > m_surrogateSet; // verify content
	std::set<std::pair<int, int> >::iterator m_surrogateItr;

};

}

#endif /* TIMEPIXCALIBRATIONHANDLER_H_ */
