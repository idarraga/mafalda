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

#include "TimePixCalibrationHandler.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>


namespace MAFTools {

///
/// Timepix calibration handler.
/// Reads ascci files for parameters a,b,c,t of the surrogate
/// function of the form: TOT = aE + b - c/(E-t), for all pixels.
/// The class can be queried to obtain E(TOT) for any pixel.
///

TimePixCalibrationHandler::TimePixCalibrationHandler(const char * af,
		const char * bf,
		const char * cf,
		const char * tf,
		int width, int height){

	// Set the Log service
	Log.setAlgoName("CalibrationHandler");
	Log.OutputLevel = MSG::INFO;

	m_width = width;
	m_height = height;

	ifstream afs(af, fstream::in);
	ifstream bfs(bf, fstream::in);
	ifstream cfs(cf, fstream::in);
	ifstream tfs(tf, fstream::in);

	////////////////////////////////////
	// a
	int badcalibCntr = 0;

	m_a = new double[m_width*m_height];
	double temp = 0.;
	int i = 0;

	while ( afs.good() ) {

		afs >> temp;
		if(afs.eof()) break;
		m_a[i++] = temp;

		// If the parameter is 0 this pixel doesn't have a proper calibration
		// Reasons: 1) the pixel was masked for some reason, 2) the calibration
		//  algorithm didn't succeed.  Probably a noisy pixel not masked.
		if(temp <= 0.) {
			badcalibCntr++;
		}

	}
	afs.close();

	if (i != m_width*m_height) {
		cout << "Calibration does not seem to be complete for paremeter a." << endl;
		cout << "Got " << i << " items.  " << m_width << "*" << m_height
				<< " were requested.  Giving up." << endl;
		exit(1);
	}

	////////////////////////////////////
	// b
	m_b = new double[m_width*m_height];
	i = 0;

	while ( bfs.good() ) {
		bfs >> temp;
		if(bfs.eof()) break;
		m_b[i++] = temp;
	}
	bfs.close();
	if(i != m_width*m_height){
		cout << "Calibration does not seem to be complete for paremeter b." << endl;
		cout << "Got " << i << " items.  " << m_width << "*" << m_height
				<< " were requested.  Giving up." << endl;
		exit(1);
	}

	////////////////////////////////////
	// c
	m_c = new double[m_width*m_height];
	i = 0;

	while ( cfs.good() ) {
		cfs >> temp;
		if(cfs.eof()) break;
		m_c[i++] = temp;
	}
	cfs.close();
	if(i != m_width*m_height){
		cout << "Calibration does not seem to be complete for paremeter c." << endl;
		cout << "Got " << i << " items.  " << m_width << "*" << m_height
				<< " were requested.  Giving up." << endl;
		exit(1);
	}

	////////////////////////////////////
	// t
	m_t = new double[m_width*m_height];
	i = 0;

	while ( tfs.good() ) {
		tfs >> temp;
		if(tfs.eof()) break;
		m_t[i++] = temp;
	}
	tfs.close();
	if(i != m_width*m_height){
		cout << "Calibration does not seem to be complete for paremeter t." << endl;
		cout << "Got " << i << " items.  " << m_width << "*" << m_height
				<< " were requested.  Giving up." << endl;
		exit(1);
	}

	m_badCalibCntr = badcalibCntr;

	double percentage = ((double)m_badCalibCntr/(double)(m_width*m_height)) * 100.;
	cout.precision(1);
	Log << MSG::INFO << "Number of pixels with incorrect calibration : " << m_badCalibCntr
			<< " (" << percentage << "%) --> Calibration OK." << endreq;

}

TimePixCalibrationHandler::~TimePixCalibrationHandler(){
}

///
/// E(TOT)
///

void TimePixCalibrationHandler::GetSolutionHighEnergy(pair<int,int> pix, int tot, double & sol1, double & sol2, double & sol) {

	int index = XYtoC(pix, m_width);

	// Quadratic solution
	double totval = double(tot);
	double a = m_a[index];
	//a *= 2.1943; // @ high energy !
	//a *= 0.95; // @ high energy !

	double b = m_b[index] - m_a[index]*m_t[index] - totval;
	double c = -m_c[index] - m_t[index]*m_b[index] + totval*m_t[index];

	sol1 = -b + sqrt(b*b - 4*a*c);
	sol1 /= 2*a;
	sol2 = -b - sqrt(b*b - 4*a*c);
	sol2 /= 2*a;

	//cout << "__________________H___________________" << endl;
	//cout << "sol1 : " << sol1 << " | sol2 : " << sol2 << endl;

	if (sol1 > 0. && sol2 > 0.) { // If both solution are positive
		double maxsol = sol1;
		if(sol2 > maxsol) maxsol = sol2;
		sol = maxsol;
	} else if(sol2 <= 0 && sol1 > 0.) {
		sol = sol1; // Otherwise use the positive solution
	} else sol = sol2;

}


void TimePixCalibrationHandler::GetCuadraticSolutions(pair<int,int> pix, int tot, double & sol1, double & sol2, double & sol) {

	int index = XYtoC(pix, m_width);

	if( m_a[index] == 0. ) { // this happens when the calibration is bad for this particular pixel
		sol1 = sol2 = sol = -1.;
		return;
	}

	// Quadratic solution
	double totval = double(tot);
	double a = m_a[index];
	double b = m_b[index] - m_a[index]*m_t[index] - totval;
	double c = -m_c[index] - m_t[index]*m_b[index] + totval*m_t[index];

	sol1 = -b + sqrt(b*b - 4*a*c);
	sol1 /= 2*a;
	sol2 = -b - sqrt(b*b - 4*a*c);
	sol2 /= 2*a;

	//cout << "______________________________________" << endl;
	//cout << "sol1 : " << sol1 << " | sol2 : " << sol2 << endl;

	if (sol1 > 0. && sol2 > 0.) { // If both solution are positive
		double maxsol = sol1;
		if(sol2 > maxsol) maxsol = sol2;
		sol = maxsol;
	} else if(sol2 <= 0 && sol1 > 0.) {
		sol = sol1; // Otherwise use the positive solution
	} else sol = sol2;

}

double TimePixCalibrationHandler::GetE(pair<int,int> pix, int tot){

	double sol1, sol2, sol;
	GetCuadraticSolutions(pix, tot, sol1, sol2, sol);

	//cout << " sol1, sol2, sol = " << sol1 << ", " << sol2 << ", " << sol << endl;

	if(sol != sol) return 0.; // avoid nan

	return sol;

	/*
	// Effective approach based on C12 studies
	if(sol > 1500) { // high energy regime > 2 MeV / pixel
		// Get an increase due to cross-talk
		// tot = tot * 1.15;
		GetSolutionHighEnergy(pix, tot, sol1, sol2, sol);
		//cout << " High charge --> " << sol << endl;
	}
	//cout << "--------------------------------------" << endl;

	return sol;
	*/
}

int TimePixCalibrationHandler::GetTOT (pair<int, int> pix, double E) {

	int index = XYtoC(pix, m_width);
	double par[] = { m_a[index], m_b[index], m_c[index], m_t[index] };

	return SurrogateCalibFunction(E, par);
}

double TimePixCalibrationHandler::SurrogateCalibFunction(double x, double * par){
	return par[0]*x + par[1] - (par[2]/(x-par[3]));
}

#ifndef __MAFALDA_EMBEDDED
///
/// Only if needed one can fetch a few surrogate functions as TF1 objects
/// Not included for embedded applications.
///
TF1 * TimePixCalibrationHandler::GetSurrogateTF1(pair<int, int> pix){

	// check first if the function has been defined.  Use the Set.
	m_surrogateItr = m_surrogateSet.find(pix);

	// If surrogate TF1 does not exist yet, build the function first
	if( m_surrogateItr == m_surrogateSet.end() ) {
		TString ftitle = "surrogate_";
		ftitle += pix.first; ftitle += "_"; ftitle += pix.second;
		// Usually the parameters come given in keV !!! WARNING !!!
		m_surrogateMap[pix] = new TF1(ftitle,"[0]*x + [1] - ([2]/(x-[3]))",0,1000); // 1 MeV
		// Populate parameters.  Look for right index first
		int index = XYtoC(pix, m_width);
		m_surrogateMap[pix]->SetParameters(m_a[index], m_b[index], m_c[index], m_t[index]);
	}

	return m_surrogateMap[pix];
}
#endif

///
/// These two functions are already defined in MAFTools but for certain tests
/// I want this class to be self-contained.
///
int TimePixCalibrationHandler::XYtoC(pair<int, int> position, int dimX) {
	return XYtoC(position.first, position.second, dimX);
}
int TimePixCalibrationHandler::XYtoC(int x, int y, int dimX) {
	return y * dimX + x;
}

} // end namespace MAFTools
