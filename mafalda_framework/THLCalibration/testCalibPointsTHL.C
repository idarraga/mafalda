 /*
 * 	Copyright 2013 John Idarraga
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


#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TH1.h>
#include <TF1.h>

#include <iostream>
#include <vector>

using namespace std;

TF1 * GetCalibFunc(pair<double, double>, pair<double, double>);
void DumpInfo(TF1 *, int, int);

void testCalibPointsTHL(){

	// 13.9  keV
	// 26.3  keV
	// 59.54 keV

	pair<double, double> thl_E_1 = make_pair (276, 13.9);

	pair<double, double> thl_E_2 = make_pair (337, 13.9);
	pair<double, double> thl_E_3 = make_pair (295, 26.3);
	pair<double, double> thl_E_4 = make_pair (257, 26.3);
	pair<double, double> thl_E_5 = make_pair (225, 26.3);
	pair<double, double> thl_E_6 = make_pair (166, 26.3);

	TF1 * fcal13 = GetCalibFunc(thl_E_1, thl_E_3);
	DumpInfo(fcal13, 438, 478-438); // operation point
	DumpInfo(fcal13, 478, 0);
	DumpInfo(fcal13, 50, 0);

	//TF1 * fcal36 = GetCalibFunc(thl_E_3, thl_E_6);
	//DumpInfo(fcal36, 438);
	//DumpInfo(fcal36, 478);


}

void DumpInfo(TF1 * f, int thl, int thl_overnoise) {

	double EnergyAtTHL = f->Eval(thl) * 1000;
	double nelectrons = EnergyAtTHL/3.6;
	double oneTHLUnit = nelectrons/thl_overnoise;

	cout << "--------------------------------------------" << endl;
	cout << "Energy at thl : " << thl << " = " << EnergyAtTHL << " eV" << endl;
	cout << " n electrons = " << nelectrons << endl;
	if(thl_overnoise > 0) {
		cout << " One THL unit = " << oneTHLUnit << " electrons" << endl;
	}
}

TF1 * GetCalibFunc(pair<double, double> p1, pair<double, double> p2) {

	double eta = (p2.second - p1.second) / (p2.first - p1.first);
	double b = p1.second - eta*p1.first;

	TF1 * f = new TF1("calibf","[0]*x + [1]",0,500);
	f->SetParameter(0, eta);
	f->SetParameter(1, b);

	return f;
}
