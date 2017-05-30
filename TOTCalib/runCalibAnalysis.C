/**
 *  Author: Jean-Samuel Roux <roux@lps.umontreal.ca>
 *  Load calibration root file for analysis and validation
 */

#include <TH2I.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>

#include <vector>
#include <string>

#include <map>
#include <iostream>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include "TExec.h"

using namespace std;

class TOTCalib;

// global
TOTCalib* calib;
Int_t palette[4] = {kWhite, kRed, kBlue, kGreen}; // status color palette


R__LOAD_LIBRARY(libTOTCalib.so);

void runCalibAnalysis (  ) {

	// Load calibration library
	gSystem->Load("libTOTCalib.so");

	int pix = 1000; // Work on this set of pixel
	int x = 20; int y = 35;
	
	TString file = "/export/home/zp/roux/github/mafalda/TOTCalib/macro_test.root";//
	LoadFile(file);


	//DrawFullPixelCalib(pix);	// pixel ID
	////DrawFullPixelCalib(x,y); 	// same, but (x,y)
//
	//DrawSpectrum(pix, "Am241"); // source can be modified
	////DrawSpectrum(x,y, "Am241");
//
	//DrawSurrogate(pix);
	////DrawSurrogate(x,y);
//
	DrawParameterMap("a");// use surrogate parameters a,b,c or t
	DrawStatusMap();
	DrawChiSquareMap();

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
		tempSource -> DumpSpectrumVectorFromSavedFile(*spectrum);
		tempSource -> CreateCalibHandlerFromSavedFile( (*it) );
		tempSource -> DumpSourceInfoFromSavedFile(*CHpoints, *CHregions, *max);

		calib->AddSingleSourceFromSavedFile(tempSource);

	}
	calib->CreateGausAndLowEnergFitFunction();

}

void DrawFullPixelCalib(int pix){
	calib->DrawFullPixelCalib(pix);
}

void DrawFullPixelCalib(int x, int y){
	calib->DrawFullPixelCalib(x, y);
}

void DrawSurrogate(int x, int y){
	DrawSurrogate(calib->XYtoX(make_pair(x,y), calib->GetMatrixWidth()));
}

void DrawSurrogate(int pix){ // from DrawFullPixelCalib

	TString cname = "Surr_pix_";
	cname += pix;
	pair<int, int> pix_xy = calib->XtoXY(pix, calib->GetMatrixWidth());

	TString ctitle = "Surrogate - pixel ";
	ctitle += pix;
	ctitle += " (";
	ctitle += pix_xy.first;
	ctitle += ",";
	ctitle += pix_xy.second;
	ctitle += ")";

	cout << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "-- Requesting surrogate function for pixel " << pix << "(" << pix_xy.first << ","
			<< pix_xy.second << ")" <<  endl;

	TCanvas * c1 = new TCanvas(cname, ctitle);

	// Graph
	TGraphErrors * g = calib->GetCalibGraph(pix);
	g->Draw("A*");

	TF1 * s = calib->GetSurrogateFunction(pix);
	s->Draw("same");
	s->GetXaxis()->SetTitle("Energy (keV)");
	s->GetYaxis()->SetTitle("TOT");
	s->SetNpx(1000);

	// Set the maximum in Y
	double * ym = g->GetY();
	vector<double> ym_v;
	for(int i = 0 ; i < g->GetN() ; i++) ym_v.push_back(ym[i]);
	double maxel_y = *max_element( ym_v.begin(), ym_v.end() );
	s->GetYaxis()->SetRangeUser( 0., maxel_y*1.1 );
	g->GetYaxis()->SetRangeUser( 0., maxel_y*1.1 );


	// Set the maximum in X
	double * xm = g->GetX();
	vector<double> xm_v;
	for(int i = 0 ; i < g->GetN() ; i++) xm_v.push_back(xm[i]);
	double maxel_x = *max_element( xm_v.begin(), xm_v.end() );
	//s->GetXaxis()->SetRangeUser( 0., maxel_x*1.1 );
	TLatex * l2 = new TLatex();

	TString parS;
	int npar = s->GetNpar();
	for(int i = 0 ; i < npar ; i++) {
		parS = "";
		if ( i == 0 ) parS = "a = ";
		if ( i == 1 ) parS = "b = ";
		if ( i == 2 ) parS = "c = ";
		if ( i == 3 ) parS = "t = ";
		parS += TString::Format("%.2f"/* +/- %.2f*/, s->GetParameter(i), s->GetParError(i) ); 	// error on fit parameters is not computed
		l2->DrawLatex(maxel_x/2, maxel_y * (1 - (i/10.)), parS); 
	}

	c1->Update();
}

void DrawStatusMap(){
	TString cname = "status_map";
	TString ctitle = "Calibration status";
	TCanvas * c1 = new TCanvas(cname, ctitle);

	cout << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "-- Requesting calibration status map --" <<  endl;
	cout << "STATUS LEGEND: " << endl;
	cout << "> kWhite: No calibration was performed on this pixel." << endl;
	cout << "> kRed  : No data, pixel is probably dead. No surrogate parameters." << endl;
	cout << "> kBlue : Missing info from one or many sources. Calibration may be incorrect." << endl;
	cout << "> kGreen: Calibration was succesful." << endl;

	int width = calib->GetMatrixWidth();
	int height = calib->GetMatrixHeight();
	int size = calib->GetMatrixSize();

	TH2I * h = new TH2I(cname, ctitle, width, 0, width, height, 0, height);
	map<int, int> status = calib->GetSurrogateStatusMap();

	int it;
	for (it = 0; it < size; it++){
		pair<int, int> pix_xy = calib->XtoXY(it, calib->GetMatrixWidth());
		if (status.find(it) == status.end()) {
			h->SetBinContent(pix_xy.first+1,pix_xy.second+1, -2);  // no calibration made on this pixel
		} else {
			h->SetBinContent(pix_xy.first+1, pix_xy.second+1, status[it]);	//binID has to have +1 because there is no bin/row/column 0

		}
	}

	c1->cd();
	h->SetStats(kFALSE);

	TExec *ex1 = new TExec("ex1","gStyle->SetPalette(4, palette);");
	h->Draw("col");
	ex1->Draw();
	h->Draw("col same"); // must draw twice to set the right color palette
}

void DrawParameterMap(string p){
	int index =-1;
	if (p == "a") {index = 0;}
	else if (p == "b"){index = 1;}
	else if (p == "c"){index = 2;}
	else if (p == "t"){index = 3;}
	else {
		cout << "[ERROR] Parameter " << p << " unrecognized."<< endl;
		cout << "        Please use one of the following parameters: a, b, c, t." << endl;
		return;
	}

	TString cname = p+"_param_map";
	TString ctitle = "Surrogate parameter "+p;
	TCanvas * c1 = new TCanvas(cname, ctitle);

	cout << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "-- Requesting surrogate parameter " << p << " map --" <<  endl;


	int width = calib->GetMatrixWidth();
	int height = calib->GetMatrixHeight();
	int size = calib->GetMatrixSize();

	TH2D * h = new TH2D(cname, ctitle, width, 0, width, height, 0, height);
	map<int, vector<double> > map = calib->GetSurrogateParamMap();

	int it;
	vector<double> v;
	for (it = 0; it < size; it++){
		pair<int, int> pix_xy = calib->XtoXY(it, calib->GetMatrixWidth());
		v = map[it];
		if (v.empty()){
			h->SetBinContent(pix_xy.first+1, pix_xy.second+1, 0.);	//binID has to have +1 because there is no bin/row/column 0
		} else{
			h->SetBinContent(pix_xy.first+1, pix_xy.second+1, v[index]);
		}
	}
	c1->cd();
	h->SetStats(kFALSE);
	TExec *ex2 = new TExec("ex2","gStyle->SetPalette(1);");
	h->Draw("colz");
	ex2->Draw();
	h->Draw("colz same"); // must draw twice to set the right color palette
}

void DrawSpectrum(int x, int y, TString s){
	DrawSpectrum( calib->XYtoX( make_pair(x,y), calib->GetMatrixWidth()) , s);
}


void DrawSpectrum(int pix, TString s){
	bool source_identified = false;
	TOTCalib * source;
	TString name;

	vector<TOTCalib*> sourVec = calib->GetSourcesVector();
	int it; // iterator

	int orderCntr = 0;
	for( it = 0; it < sourVec.size(); it++){
		name = (sourVec[it])->GetCalibHandler()->GetSourcename();
		if (name == s){
			source = sourVec[it];
			source_identified = true;
			break;
		} else {
			orderCntr += ( sourVec[it]->GetCalibHandler()->GetCalibPoints() ).size();
		}
	}

	if(!source_identified){
		cout << "[ERROR] Source " << s.Data() << " unrecognized. Sources in your file are:" << endl;
		for(it = 0; it < sourVec.size(); it++){
			cout << sourVec[it]->GetCalibHandler()->GetSourcename() << "   " ; }
		cout << endl;
		return;
	}

	vector<double> fit_const = calib->GetMapCalibPointsConstants()[pix];
	vector< pair<double, double> > fit_mean = calib->GetMapCalibPoints()[pix];
	vector<double> fit_sigmas = calib->GetMapCalibSigmas()[pix];
	vector<double> fit_ia = calib->GetMapCalibIA()[pix];
	vector<double> fit_ib = calib->GetMapCalibIB()[pix];
	vector<double> fit_ic = calib->GetMapCalibIC()[pix];
	vector<double> fit_it = calib->GetMapCalibIT()[pix];

	// this part is from DrawFullPixelCalib
	TString cname = "pix_";
	cname += pix;
	cname += "_";
	cname += name;
	pair<int, int> pix_xy = calib->XtoXY(pix, calib->GetMatrixWidth());

	TString ctitle = "Data pixel ";
	ctitle += pix;
	ctitle += " (";
	ctitle += pix_xy.first;
	ctitle += ",";
	ctitle += pix_xy.second;
	ctitle += ") - ";
	ctitle += name;

	cout << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "-- Requesting " << name <<" spectrum for pixel " << pix << "(" << pix_xy.first << ","
			<< pix_xy.second << ")" <<  endl;

	TCanvas * c1 = new TCanvas(cname, ctitle);

	TH1 * h; TF1 * kf;

	TLegend * leg1 = new TLegend(0.6, 0.6, 0.9, 0.9);
	leg1->SetFillColor(kWhite);
	leg1->SetBorderSize(1);

	TF1 * gf; TF1 * gf_clone;

	// Data
	h = source->GetHisto(pix, "summary");
	h->Draw("HIST");
	h->GetXaxis()->SetTitle("TOT");
	h->GetYaxis()->SetTitle("entries");

	// Kernel density
	kf = source->GetKernelDensityFunction(pix);
	kf->Draw("same");
	kf->SetLineColor(kBlack);
	kf->SetLineStyle(2);
	kf->SetLineWidth(1);

	kf->GetXaxis()->SetTitle("TOT");
	kf->GetYaxis()->SetTitle("entries");

	leg1->AddEntry(kf, "k.d.f. (envelope)", "L");
	leg1->AddEntry(h, "Calibration data", "L");

	TLatex * l1 = new TLatex();
	l1->DrawLatex(source->GetNBins()/2, h->GetMaximum() * 0.5 , name);

	// Get list of peaks identified
	map<int, double> calibPoints = source->GetCalibHandler()->GetCalibPoints();
	vector<double> peaks = (source->GetMaxPeaksIdentified())[pix];

	// Points used in the fit
	vector< pair<double, double> > calibFitPoints = calib->GetCalibPoints(pix);

	int nCalibPoints = (int)calibPoints.size();
	double pos_offset = 1.;
	// Check if it's possible to draw
	if ( (int) peaks.size() < nCalibPoints ) {
		cout << "Not enough peaks were identified for this pixel !" << endl;
		return;
	}

	// loop over calip points per source
	for(int p = 0 ; p < nCalibPoints ; p++) {
		if(calibPoints[p]<0){ // skips energies < 0 (ignored points)
			orderCntr++;
			continue;
		}

		// Fit on data
		TString funcName = "f_histofit_";
		funcName += pix;
		funcName += "_";
		funcName += p;
		gf = calib->FittingFunctionSelector( calibPoints[p], source, p );
         
        gf->FixParameter( 0, fit_const[orderCntr] );
        gf->FixParameter( 2, fit_sigmas[orderCntr] );
        if( TString(gf->GetName()).Contains("gf_lowe") ) {
            gf->FixParameter( 1, fit_mean[orderCntr].first ); // fixed energy
            gf->FixParameter( 3, fit_ia[orderCntr] );
            gf->FixParameter( 4, fit_ib[orderCntr] );
            gf->FixParameter( 5, fit_ic[orderCntr] );
            gf->FixParameter( 6, fit_it[orderCntr] );

        }else{
            gf->SetParameter( 1, fit_mean[orderCntr].second ); // mean tot from fit
        }
        
         
		cout<<"--------- Source: " << name << " ---------" <<endl << "The fitting function for this source is : " << gf->GetName() << endl;
		//gf->Dump();
		gf_clone = static_cast<TF1 * > ( gf->Clone( funcName ) );
		//gf_clone->Dump();

		cout << "Point : " << calibPoints[p] << " | order : " << orderCntr << " [ ";
         
         if( TString(gf->GetName()).Contains("gf_lowe") ) {

             cout << "Constant: "<< fit_const[orderCntr] << ", " << "Mean: " <<fit_mean[orderCntr].first << ", " <<"Sigma: "<< fit_sigmas[orderCntr] << ", ";
             cout <<"a: "<< fit_ia[orderCntr] << ", " <<"b: "<< fit_ib[orderCntr] << ", " <<"c: "<< fit_ic[orderCntr] << ", " <<"t: "<< fit_it[orderCntr];

         }else{

             cout << "Constant: "<< fit_const[orderCntr] << ", " << "Mean: " <<fit_mean[orderCntr].second << ", " <<"Sigma: "<< fit_sigmas[orderCntr] << ", ";

         }
		cout << " ]" << endl;
         
		gf_clone->SetLineColor(kRed);
		gf_clone->Draw("same");              

		leg1->AddEntry(gf_clone, "Fit on spectrum", "L");


		//TString peak = TString::Format( "%.1f  -> %.3f keV", calibFitPoints[orderCntr].second, calibPoints[p] );
		TString peak = TString::Format( "%.1f  -> %.3f keV", fit_mean[orderCntr].second, calibPoints[p] );
		l1->DrawLatex( peaks[p], h->GetMaximum() * pos_offset , peak );

		pos_offset -= 0.2;

		orderCntr++;
	}

	leg1->Draw();
	c1->Update();

}


void DrawChiSquareMap(){
	TString cname = "ChiSquare_map";
	TString ctitle = "Surrogate fit - Chi2";
	TCanvas * c1 = new TCanvas(cname, ctitle);

	cout << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "-- Requesting surrogate chi square map --" <<  endl;


	int width = calib->GetMatrixWidth();
	int height = calib->GetMatrixHeight();
	int size = calib->GetMatrixSize();

	TH2D * h = new TH2D(cname, ctitle, width, 0, width, height, 0, height);
	map<int, vector<double> > map = calib->GetSurrogateParamMap();

	TF1 * f = new TF1("surrogate", surrogatefunc_calib, calib->GetGlobalThresholdEnergyAndErr().first, 60.0, 4); // surrogate to get chi2


	int pix;
	vector<double> v;
	for (pix = 0; pix < size; pix++){
		pair<int, int> pix_xy = calib->XtoXY(pix, calib->GetMatrixWidth());
		v = map[pix];
		if (v.empty()){
			h->SetBinContent(pix_xy.first+1, pix_xy.second+1, 0.);
		} else {
			TGraphErrors * g = calib->GetCalibGraph(pix);
			for(int i = 0; i < 4; i++) f->FixParameter(i, v[i]);
			TFitResultPtr fitr = g->Fit("surrogate", "NRSQ", ""); // not an actual fit since every parameter is fixed
			h->SetBinContent(pix_xy.first+1, pix_xy.second+1, fitr.Get()->Chi2());
		}

	}
	c1->cd();

	h->SetStats(kFALSE);
	TExec *ex2 = new TExec("ex2","gStyle->SetPalette(1);");
	h->Draw("colz");
	ex2->Draw();
	h->Draw("colz same"); // must draw twice to set the right color palette
}


double surrogatefunc_calib(double * x, double * par) {

	// independent var
	double xx = x[0];

	// pars
	double a = par[0];
	double b = par[1];
	double c = par[2];
	double t = par[3];

	double func = a * xx + b;
	func -= ( c / ( xx - t) );

	return func;
}