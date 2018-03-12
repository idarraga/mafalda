#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>

using namespace std;

TH1 * GetHisto(TString fn, TString tn, TString bn, TString cut, int, int, int, int);

void runCalibrationValidation(){

	TString file1 = "MAFOutput_MPXNtuple_Am241_dataset01_20MHz.root";
	//TString file1 = "MAFOutput_MPXNtuple_Fe55_dataset01_20MHz.root";
	int rebin_factor = 2;

	TString theCut_s = "";
	TString theCut_d = ""; //"doubleHitE > 20.0 && doubleHitE < 36.0";
	TString theCut_t = ""; //"doubleHitE > 20.0 && doubleHitE < 36.0";

	// TOT
	TH1 * h_tot_s = GetHisto(file1, "CalibrationValidation", "singleHitTOT", theCut_s, 0, 150, kRed, rebin_factor);
	TH1 * h_tot_d = GetHisto(file1, "CalibrationValidation", "doubleHitTOT", theCut_d, 0, 150, kOrange, rebin_factor);
	//TH1 * h_tot_d1 = GetHisto(file1, "CalibrationValidation", "doubleHitTOT_1", theCut_d, 0, 150, kBlue, rebin_factor);
	//TH1 * h_tot_d2 = GetHisto(file1, "CalibrationValidation", "doubleHitTOT_2", theCut_d, 0, 150, kGreen, rebin_factor);
	TH1 * h_tot_t = GetHisto(file1, "CalibrationValidation", "tripleHitTOT", theCut_t, 0, 150, kBlack, rebin_factor);
	TH1 * h_tot_q = GetHisto(file1, "CalibrationValidation", "quadHitTOT", "", 0, 150, kBlue, rebin_factor);
	TH1 * h_tot_all = GetHisto(file1, "CalibrationValidation", "allTOT", "", 0, 150, kGreen, rebin_factor);

	// E
	TH1 * h_E_s = GetHisto(file1, "CalibrationValidation", "singleHitE", theCut_s, 0, 150, kRed, rebin_factor);
	TH1 * h_E_d = GetHisto(file1, "CalibrationValidation", "doubleHitE", theCut_d, 0, 150, kOrange, rebin_factor);
	//TH1 * h_E_d1 = GetHisto(file1, "CalibrationValidation", "doubleHitE_1", theCut_d, 0, 150, kBlue, rebin_factor);
	//TH1 * h_E_d2 = GetHisto(file1, "CalibrationValidation", "doubleHitE_2", theCut_d, 0, 150, kGreen, rebin_factor);
	TH1 * h_E_t = GetHisto(file1, "CalibrationValidation", "tripleHitE", theCut_t, 0, 150, kBlack, rebin_factor);
	TH1 * h_E_q = GetHisto(file1, "CalibrationValidation", "quadHitE", "", 0, 150, kBlue, rebin_factor);
	TH1 * h_E_all = GetHisto(file1, "CalibrationValidation", "allE", "", 0, 150, kGreen, rebin_factor);

	TCanvas * c1 = new TCanvas("TOT_Spectrum", "TOT_Spectrum");
	c1->cd();

	h_tot_s->SetTitle("Am241");

	h_tot_s->Draw();
	h_tot_d->Draw("same");
	//h_tot_d1->Draw("same");
	//h_tot_d2->Draw("same");
	h_tot_t->Draw("same");
	h_tot_q->Draw("same");
	//h_tot_all->Draw("same");

	TLegend * leg1 = new TLegend(0.6, 0.6, 0.9, 0.9);
	leg1->AddEntry(h_tot_s, "single", "L");
	leg1->AddEntry(h_tot_d, "double", "L");
	//leg->AddEntry(h_E_d1, "double 1", "L");
	//leg->AddEntry(h_E_d2, "double 2", "L");
	leg1->AddEntry(h_tot_t, "triple", "L");
	leg1->AddEntry(h_tot_q, "quad", "L");
	//leg->AddEntry(h_E_q, "quad", "L");
	//leg1->AddEntry(h_tot_all, "all", "L");
	leg1->Draw();

	TCanvas * c2 = new TCanvas("E_Spectrum", "E_Spectrum");
	c2->cd();

	h_E_s->Draw();
	h_E_d->Draw("same");
	//h_E_d1->Draw("same");
	//h_E_d2->Draw("same");
	h_E_t->Draw("same");
	//h_E_q->Draw("same");
	h_E_all->Draw("same");

	TLegend * leg = new TLegend(0.6, 0.6, 0.9, 0.9);
	leg->AddEntry(h_E_s, "single", "L");
	leg->AddEntry(h_E_d, "double", "L");
	//leg->AddEntry(h_E_d1, "double 1", "L");
	//leg->AddEntry(h_E_d2, "double 2", "L");
	leg->AddEntry(h_E_t, "triple", "L");
	//leg->AddEntry(h_E_q, "quad", "L");
	leg->AddEntry(h_E_all, "all", "L");

	leg->Draw();

	//h_E_t->Draw("same");
	//h_E_q->Draw("same");

}


TH1 * GetHisto(TString fn, TString tn, TString bn, TString cut, int min, int max, int color, int rebin) {

	TFile * f = new TFile(fn);
	TTree * T = static_cast<TTree * >( f->Get(tn) );

	TString histoName = "h_";
	histoName += bn;

	TString drawS = bn;
	drawS += ">>";
	drawS += histoName;
	drawS += "(";
	drawS += max-min; drawS += ",";
	drawS += min; drawS += ",";
	drawS += max; drawS += ")";

	cout << "Draw : " << drawS << endl;

	T->Draw(drawS, cut, "", 10000, 0);

	TH1 * h = static_cast<TH1 *> ( gROOT->FindObject(histoName) );
	h->SetLineColor(color);

	if(rebin > 1) h->Rebin(rebin);

	return h;
}
