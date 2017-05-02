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

R__LOAD_LIBRARY(libTOTCalib) // ROOT6 (cling)

void runSavePixelResolution() {

    nTotalFrames = -1;            	 // -1: run all frames in input dataset
    //minpix = 0; maxpix = 256*256-1;  // Work on this set of pixels
    minpix = 32741; maxpix = 32842;  // Work on this set of pixels

    // Load calibration library
    gSystem->Load("libTOTCalib.so");

    TString folder_local = "/Users/thomasbilloud/workspace/DATA/MAFOutput/2017-02-GaAs500-E06W0203/";
    
    // Process one source ///////////////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pCu = new TOTCalib(folder_local+"MAFOutput_MPXNtuple_Cu_-50V.root","Cu-fluo", minpix, maxpix, 200, nTotalFrames);
    pCu->SetVerboseLevel(TOTCalib::__VER_INFO);
    pCu->SetKernelBandWidth(20.); // Kernel Bandwidth
    pCu->Loop(); // Run this source
    pCu->SavePixelResolution();    


}
