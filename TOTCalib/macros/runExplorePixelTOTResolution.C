/*
 * Author: Thomas Billoud <billoud@lps.umontreal.ca>
 * 
 * File to display pixel spectra when mouse is pointing at a particuler in
 * a map from output_pixelResolution root file (e.g. "SingleHitFitMeans").
 *
 * Use it in CINT as (you need to be in TOTCalib folder):
 * root [] TFile f("output_pixelResolution.root");
 * root [] SingleHitFitMeans->Draw("colz");
 * root [] gStyle->SetOptStat(0);
 * root [] c1->AddExec("",".x runExplorePixelTOTResolution.C");
 * root [] c1->ToggleEventStatus(); // to display pixel coordinate in the status bar (bottom line of canvas)
 * root [] gStyle->SetPalette(52); // Easier to distinguish abnormal fits by eye with black & white picture.
 *
 * Inspired from exec2.C in root examples
 * Warning: no cout here !
*/

int XYtoX(int pixX, int pixY, int dimX){
    return pixY * dimX + pixX;
}

double fit_func_lowen(double * x, double * par) {

    // Function proposed in J. Jakubek / Nuclear Instruments and Methods in Physics Research A 633 (2011) S262–S266

    // independent var
    Double_t xx = x[0]; // TOT

    // parameters for gaussian
    Double_t gconst = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];

    // surrogate ^ -1
    Double_t a = par[3];
    Double_t b = par[4];
    Double_t c = par[5];
    Double_t t = par[6];

    // Inverse of the surrogate ( = energy in keV)
    Double_t surrogate_inverse, delta, num1, num2, denum;
    delta = TMath::Power((b-a*t-xx),2) - 4 * a * (xx*t-c-b*t);
    num1 = a*t + xx - b;
    num2 = TMath::Sqrt(delta) ;
    denum = 2 * a;
    surrogate_inverse = (num1 +  num2) / denum;

    // Gaussian of the inversed surrogate
    Double_t arg = (surrogate_inverse - mean)/sigma;
    Double_t func = gconst * TMath::Exp(-0.5*arg*arg);

    return func;
}



void runExplorePixelTOTResolution()
{

    //----------------------- Look for trees and branches in root file --------------------------------    
    
    TTree *T = (TTree*)f.Get("SavePixelResolution");
    Int_t pixelID;
    TH1I *hpx = 0;
    TF1 *kernel_func = 0;
    vector<double> *br_double_sigmafit=0;
    vector<double> *br_double_totmeanfit=0;
    vector<double> *br_double_totmeanfitError=0;
    vector<double> *br_double_constantfit=0;
    vector<double> *br_double_Chi2fit=0;
    vector<double> *br_double_NDFfit=0;    
    vector<double> *br_double_afit=0;
    vector<double> *br_double_bfit=0;
    vector<double> *br_double_cfit=0;
    vector<double> *br_double_tfit=0;
    
    Int_t selectedPeakID;
    T->SetBranchAddress("pixelID",&pixelID);
    T->SetBranchAddress("selectedPeakID",&selectedPeakID);
    T->SetBranchAddress("Histo_Spectrum",&hpx);
    T->SetBranchAddress("Kernel_Function",&kernel_func);
    T->SetBranchAddress("FitSigma",&br_double_sigmafit);
    T->SetBranchAddress("FitMean",&br_double_totmeanfit);
    T->SetBranchAddress("FitMeanError",&br_double_totmeanfitError);
    T->SetBranchAddress("FitConstant",&br_double_constantfit);
    T->SetBranchAddress("FitChi2",&br_double_Chi2fit);
    T->SetBranchAddress("FitNDF",&br_double_NDFfit);    
    T->SetBranchAddress("Fita",&br_double_afit);
    T->SetBranchAddress("Fitb",&br_double_bfit);
    T->SetBranchAddress("Fitc",&br_double_cfit);
    T->SetBranchAddress("Fitt",&br_double_tfit);
        
    // In case calibrated data is present    
    vector<double> *br_double_sigmafit_calibrated=0; 
    vector<double> *br_double_totmeanfit_calibrated=0;
    vector<double> *br_double_constantfit_calibrated=0;
    TH1D *hpx_calibrated = 0;
    static TString calibBranch("FitConstant_calibrated");
    TObject* calib_branch = T->GetListOfBranches()->FindObject(calibBranch);
    bool calib_present = false;
    if (calib_branch!=nullptr){ calib_present = true;}
    if (calib_present){
        T->SetBranchAddress("Histo_Spectrum_calibrated",&hpx_calibrated);        
        T->SetBranchAddress("FitConstant_calibrated",&br_double_constantfit_calibrated);
        T->SetBranchAddress("FitMean_calibrated",&br_double_totmeanfit_calibrated);
        T->SetBranchAddress("FitSigma_calibrated",&br_double_sigmafit_calibrated);         
    }  
    

    //------------------------------ Plot interactive histograms ----------------------------------------
    
    
    T->GetEntry(0);
    Int_t FirstPixelInTree = pixelID;

    TObject *select = gPad->GetSelected();
    if(!select) return;
    if (!select->InheritsFrom(TH2::Class())) {gPad->SetUniqueID(0); return;}
    gPad->GetCanvas()->FeedbackMode(kTRUE);

    //erase old position and draw a line at current position
    int pyold = gPad->GetUniqueID();
    int px = gPad->GetEventX();
    int py = gPad->GetEventY();
    float uxmin = gPad->GetUxmin();
    float uxmax = gPad->GetUxmax();
    int pxmin = gPad->XtoAbsPixel(uxmin);
    int pxmax = gPad->XtoAbsPixel(uxmax);
    if(pyold) gVirtualX->DrawLine(pxmin,pyold,pxmax,pyold);
    gVirtualX->DrawLine(pxmin,py,pxmax,py);
    gPad->SetUniqueID(py);
    Float_t upx = gPad->AbsPixeltoX(px);
    Float_t x = gPad->PadtoX(upx);
    Float_t upy = gPad->AbsPixeltoY(py);
    Float_t y = gPad->PadtoY(upy);

    //create or set the new canvas c2
    TVirtualPad *padsav = gPad;
    TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
    if(c2) {delete c2->GetPrimitive("Projection");} else {c2 = new TCanvas("c2");}
    c2->cd();

    //draw slice corresponding to mouse position
    TH2 *h = (TH2*)select;
    Int_t binx = h->GetXaxis()->FindBin(x);
    Int_t biny = h->GetYaxis()->FindBin(y);
    Int_t PixCoorX = binx -1 ;
    Int_t PixCoorY = biny -1 ;

    // Load root file entry of pointed pixel
    int PixelToPlot = XYtoX(PixCoorX,PixCoorY,256);
    Int_t nevent = T->GetEntries();
    T->GetEntry(PixelToPlot-FirstPixelInTree); // I soustract firstPixelInTree for cases where root file starts with a pixelID different than 0.

    // First, draw histogram and kernel function  (stored in root file)
    hpx->GetXaxis()->SetTitle("TOT");
    hpx->GetYaxis()->SetTitle("Counts");    
    hpx->Draw();
    kernel_func->Draw("same");
    kernel_func->SetLineColor(kBlack);
    kernel_func->SetLineStyle(2);
    kernel_func->SetLineWidth(1);

    // Secondly, draw fit functions (from parameters stored in root file)
    TF1 *fit_func = new TF1("gf_linear", "gaus(0)", 0.,hpx->GetNbinsX());
    TF1 *fit_func_lowenergy = new TF1("gf_lowen",fit_func_lowen, 0.,hpx->GetNbinsX(),7);    
    Int_t number_of_fitted_peaks = br_double_constantfit->size();
    for(Int_t i=0 ; i<number_of_fitted_peaks ; i++ ) {

       if (br_double_afit->at(i) == 0.0){
           fit_func->SetParameter( 0, br_double_constantfit->at(i) ); // I don't get the use of orderCntr in DrawFullPixelCalib
           fit_func->SetParameter( 1, br_double_totmeanfit->at(i) );
           fit_func->SetParameter( 2, br_double_sigmafit->at(i));
           if (i==selectedPeakID) {
                fit_func->SetLineColor(kGreen);
                double chi2 = br_double_Chi2fit->at(i);
                TString chi2_str = to_string(chi2);
                double NDF = br_double_NDFfit->at(i);
                TString NDF_str = to_string(NDF);
                double fitmeanError = br_double_totmeanfitError->at(i);
                TString MeanError_str = to_string(fitmeanError);                
                TString str = "Chi2/NDF: "+chi2_str+" / "+NDF_str+" --- Mean error: "+MeanError_str;                
                TText *t = new TText(.5,.5,str);
                //t->SetTextAlign(22);
                t->SetTextFont(43);
                t->SetTextSize(20);
                t->SetTextColor(kGreen);
                t->SetBBoxCenterX(200);
                t->SetBBoxCenterY(100);                
                t->Draw();
           }else{
               fit_func->SetLineColor(kRed);
           }
           fit_func->DrawCopy("same");  
       }else{
           fit_func_lowenergy->SetParameter( 0, br_double_constantfit->at(i) ); // I don't get the use of orderCntr in DrawFullPixelCalib
           fit_func_lowenergy->SetParameter( 1, br_double_totmeanfit->at(i) );
           fit_func_lowenergy->SetParameter( 2, br_double_sigmafit->at(i));
           fit_func_lowenergy->SetParameter( 3, br_double_afit->at(i));
           fit_func_lowenergy->SetParameter( 4, br_double_bfit->at(i));
           fit_func_lowenergy->SetParameter( 5, br_double_cfit->at(i));
           fit_func_lowenergy->SetParameter( 6, br_double_tfit->at(i));
                      
           if (i==selectedPeakID) {
               fit_func_lowenergy->SetLineColor(kGreen);
           }else{
               fit_func_lowenergy->SetLineColor(kRed);
           }
           fit_func_lowenergy->DrawCopy("same");
       }

    }
    c2->Update();
    padsav->cd();
    
    
    //------------------- Add another canvas if calibrated data is present in root file ---------------------------------
    
    if (calib_present){
        
        //create or set the new canvas c2
        TVirtualPad *padsav2 = gPad;
        TCanvas *c3 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c3");
        if(c3) delete c3->GetPrimitive("Projection");
        else c3 = new TCanvas("c3");
        c3->cd();
    
        // First, draw histogram and kernel function  (stored in root file)
        hpx_calibrated->GetXaxis()->SetTitle("Energy (keV)");
        hpx_calibrated->GetYaxis()->SetTitle("Counts"); 
        hpx_calibrated->Draw();
        
        // Secondly, draw fit functions (from parameters stored in root file)
        TF1 *fit_func_calibrated = new TF1("gf_linear", "gaus(0)", 0.,hpx->GetNbinsX());
        Int_t number_of_fitted_peaks = br_double_constantfit_calibrated->size();
        for(Int_t i=0 ; i<number_of_fitted_peaks ; i++ ) {
    
           fit_func_calibrated->SetParameter( 0, br_double_constantfit_calibrated->at(i) ); // I don't get the use of orderCntr in DrawFullPixelCalib
           fit_func_calibrated->SetParameter( 1, br_double_totmeanfit_calibrated->at(i) );
           fit_func_calibrated->SetParameter( 2, br_double_sigmafit_calibrated->at(i));
           if (i==selectedPeakID) {
               fit_func_calibrated->SetLineColor(kGreen);
           }else{
               fit_func_calibrated->SetLineColor(kRed);
           }
           fit_func_calibrated->DrawCopy("same");  
        } 
        
        c3->Update();
        padsav2->cd();
    }

}

