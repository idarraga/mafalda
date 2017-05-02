/**
 *  Author: Jean-Samuel Roux <roux@lps.umontreal.ca>
 *  Load calibration root file for analysis and validation
 */

#include <TH2I.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>

#include <vector>

#include <map>
#include <iostream>
#include <TROOT.h>
#include <TString.h>

using namespace std;

class TOTCalib;

// Prototypes and global
void LoadFile(TString);
TOTCalib* calib;

R__LOAD_LIBRARY(libTOTCalib.so);

void runCalibAnalysis (  ) {

	// Load calibration library
	gSystem->Load("libTOTCalib.so");

	int pix = 14000; // Work on this set of pixel
	
	TString file = "/export/home/zp/roux/github/mafalda/TOTCalib/macro_test_2nd.root";//
	LoadFile(file);

	calib->DrawFullPixelCalib(pix);

	//TO DO
	/*
	calib->DrawSurrogate(pix);
	calib->DrawSpectrum(pix, "Am241");
	calib->DrawStatusMap();
	calib->DrawParameterMap("a");
	calib->DrawChiSquaredMap();
	*/

}


void LoadFile(TString s){

	// We create an empty TOTCalib (source) object to inherit every variable and function from this class (eg DrawFullPixelCalib)
	calib = new TOTCalib();

	// Load ROOT file and link branches to variable
	TFile * f = new TFile(s.Data());
	TTree * tsurr = (TTree*) f->Get("surrogateFunction");
	TTree * tparam =(TTree*) f->Get("parameters");

	map<int, vector<double> > * surr_p;
	map<int, int> *surr_status;
	tsurr->SetBranchAddress("parameters", &surr_p);
	tsurr->SetBranchAddress("status", &surr_status);
	tsurr->GetEntry(0); // this TTree has a single entry

	map<int, vector< pair<double, double> > > * points; //energy and mean
	map<int, vector<double> > * sigmas;
	map<int, vector<double> > * constants;
	pair<double, double> * thres; //threshold + error
	vector<TString> * sn ; //source name

	// low energy fit parameter
	map<int, vector<double> > *calibPoints_ia;
	map<int, vector<double> > *calibPoints_ib;
	map<int, vector<double> > *calibPoints_ic;
	map<int, vector<double> > *calibPoints_it;

	tparam->SetBranchAddress("energyAndGausMean", &points);
	tparam->SetBranchAddress("gausSigma", &sigmas);
	tparam->SetBranchAddress("gausConst" , &constants);
	tparam->SetBranchAddress("threshold", &thres);
	tparam->SetBranchAddress("sources", &sn);
	tparam->SetBranchAddress("low_ia", &calibPoints_ia);
	tparam->SetBranchAddress("low_ib", &calibPoints_ib);
	tparam->SetBranchAddress("low_ic", &calibPoints_ic);
	tparam->SetBranchAddress("low_it", &calibPoints_it);

	// Dump variables into calib
	tparam->GetEntry(0); // this TTree has a single entry
	calib-> DumpCalibParametersFromSavedFile(*points, *sigmas, *constants, *calibPoints_ia, *calibPoints_ib, *calibPoints_ic, *calibPoints_it, *surr_p, *surr_status);
	calib-> SetGlobalThresholdEnergyAndErr((*thres).first, (*thres).second);
	
	// Create sub-TOTCalib objects from TTrees
	vector<TString>::iterator it;
	TOTCalib * tempSource;
	vector< vector<double> > * spectrum;
	map<int, double> *CHpoints; // points defined in calib handler object
	map<int, int> *CHregions; // regions defined in calib handler object
	map<int, vector<double> > *max; // identified peaks (required by DrawFullPixelCalib)
	double bandwidth;

	for (it = (*sn).begin(); it != (*sn).end(); it++){
		TTree * tsour = (TTree*) f->Get( (*it).Data() );

		tsour->SetBranchAddress("spectrumVec", &spectrum);
		tsour->SetBranchAddress("bandwidth", &bandwidth);
		tsour->SetBranchAddress("regions", &CHregions);
		tsour->SetBranchAddress("points", &CHpoints);
		tsour->SetBranchAddress("maximums", &max);
		tsour -> GetEntry(0); // this TTree has a single entry

		tempSource = new TOTCalib();
		tempSource -> SetKernelBandWidth(bandwidth);
		tempSource -> DumbSpectrumVectorFromSavedFile(*spectrum);
		tempSource -> CreateCalibHandlerFromSavedFile( (*it) );
		tempSource -> DumpSourceInfoFromSavedFile(*CHpoints, *CHregions, *max);

		calib->AddSingleSourceFromSavedFile(tempSource);

	}
	calib->CreateGausAndLowEnergFitFunction();

}
