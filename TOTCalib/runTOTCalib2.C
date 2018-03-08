/**
 *  Author: Thomas Billoud <billoud@lps.umontreal.ca>
 *  TOT calibration for a Timepix device
 *  March 2018
 */

#include <TString.h>
using namespace std;
class TOTCalib;

int nTotalFrames;
int minpix; int maxpix;

R__LOAD_LIBRARY(libTOTCalib) // ROOT6 (cling)

void runTOTCalib2() {

    nTotalFrames = -1; // -1 to run all input frames, or a specfic number of frames for fast tests
    //minpix = 0; maxpix = 256*256-1;  // Calibrate the whole detector
    minpix = 0; maxpix = 5000;  // Work on this set of pixels for fast tests

    TString folder = "/Users/thomasbilloud/workspace/DATA/MAFOutput/2017-02-GaAs500-E06W0203/";
    int verb = TOTCalib::__VER_QUIET; // __VER_INFO for info pixel by pixel
    
    // Process one source in the linear region ////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pAm = new TOTCalib(folder+"MAFOutput_MPXNtuple_Am_-300V.root","Am241", minpix, maxpix, 400, nTotalFrames); // 400 is highest TOT value in the spectra
    pAm->SetVerboseLevel(verb);
    pAm->SetHistoRebinning(2); //rebinning is important in case of noisy spectra (must be a divider of the highest defined TOT value above)
    pAm->SetKernelBandWidth(15); // must be close to the sigma of the peak 
    //pAm->OptimizeOnePeak(0.3); // optimizes peak search, to use with caution (see README) 
    pAm->Loop(); // Performs peak searching in spectra

    // Process the 2nd source in the linear region /////////////////////////////////////////////////////////////////////////////
    TOTCalib * pCd = new TOTCalib(folder+"MAFOutput_MPXNtuple_Cd_-300V.root","Cd_fluo_TPXGaAs", minpix, maxpix, 250, nTotalFrames);
    pCd->SetVerboseLevel(verb);
    pCd->SetKernelBandWidth(10);
    pCd->SetHistoRebinning(2); 
    //pCd->OptimizeOnePeak(0.3);    
    pCd->Loop();    
    
    // Process the low energy source ///////////////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pCu = new TOTCalib(folder+"MAFOutput_MPXNtuple_Cu_-300V.root","Cu_fluo", minpix, maxpix, 120, nTotalFrames);
    pCu->SetVerboseLevel(verb);
    pCu->SetHistoRebinning(2);     
    pCu->SetKernelBandWidth(10); 
    pCu->SetLowEnFit_Params(80.,10.,200.,1.); // First guess of low energy paramters: constant, sigma, c, t
    pCu->Loop(); 
    
    // Calibrate and create a, b, c, t ASCII files //////////////////////////////////////////////////////////////////////////////
    pAm->Choose2LinearPeaks(59.5409,23.1); // must be equal to peaks defined in TOTCalib.cpp->CalibHandler (in keV)
    pAm->SetLowEnergyPeak(8.05); // must be equal to the peak defined in TOTCalib.cpp->CalibHandler (in keV)
    pAm->Blender2(pCd,pCu,"test"); // last argument is the base name for ASCII files

    // Plots (uncomment for tests before calibrating the whole detector) ////////////////////////////////////////////////
    // Draw full info for a few pixels. 
    // int reportCntr = 0;
    // for(int i = minpix; i <= maxpix ; i++){
    //     pAm->DrawFullPixelCalib2( i );
    //     if(reportCntr++ > 2) break; // make only 5 reports
    // }

}
