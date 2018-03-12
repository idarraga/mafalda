/*
 *  Author: Thomas Billoud <billoud@lps.umontreal.ca>
 *  Macro to save spectra maps for one source (output root file intended to be used with runExplorePixelResolution.C)
 */

#include <TH2I.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>

#include <vector>

using namespace std;

class TOTCalib;

// Prototypes and a few harmless globals
int nTotalFrames;
int minpix; int maxpix;

R__LOAD_LIBRARY(../build/libTOTCalib) // ROOT6 (cling)

void runSavePixelResolution_Thomas() {

    nTotalFrames = -1;            	 // -1: run all frames in input dataset
    minpix = 0; maxpix = 256*256-1;  // Work on this set of pixels
    //minpix = 0; maxpix = 5000;  // Work on this set of pixels

    TString folder_local = "/Users/thomasbilloud/workspace/DATA/MAFOutput/2017-02-GaAs500-E06W0203/";
    
    // Process one source ///////////////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pCu = new TOTCalib(folder_local+"MAFOutput_MPXNtuple_Cd_-30V.root","Cd_fluo_TPXGaAs", minpix, maxpix, 200, nTotalFrames);
    pCu->SetVerboseLevel(TOTCalib::__VER_QUIET);
    pCu->SetHistoRebinning(2); //rebinning is important in case of noisy spectra (must be a divider of the highest defined TOT value above)    
    pCu->SetKernelBandWidth(10.); // Kernel Bandwidth
    pCu->OptimizeOnePeak(0.3);        
    pCu->Loop(); // Run this source
    
    // Without calibration
    //pCu->SavePixelResolution("SavePixelResolutionOutput_Cu_-50V"); // the argument is the name of the output root file.  
    
    // With calibration
    TString calib_folder = "/Users/thomasbilloud/workspace/DETECTEURS/2017_02_GaAs500_E06W0203/calib/2018_03_07_AmCdCu/";
    TString a = calib_folder+"test_a.txt";
    TString b = calib_folder+"test_b.txt";
    TString c = calib_folder+"test_c.txt";
    TString t = calib_folder+"test_t.txt";
    pCu->Choose2LinearPeaks(59.5409,23.1); // must be equal to peaks defined in TOTCalib.cpp->CalibHandler (in keV)    
    pCu->SavePixelResolution2("SavePixelResolutionOutput_Cd_-30V",a,b,c,t); // the first argument is the name of the output root file.   


}
