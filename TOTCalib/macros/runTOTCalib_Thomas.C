/**
 *  Author: John Idarraga <idarraga@cern.ch>
 *  TOT calibration for a Timepix/Timepix3 device
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

void runTOTCalib_Thomas (  ) {

    nTotalFrames = -1; // 5000            	 // -1: run all frames in input dataset
    //minpix = 0; maxpix = 256*256-1;  // Work on this set of pixels
    minpix = 0; maxpix = 5000;  // Work on this set of pixels

    TString folder_local = "/Users/thomasbilloud/workspace/DATA/MAFOutput/2017-02-GaAs500-E06W0203/";
    int verb = TOTCalib::__VER_QUIET;
    
    // Process one source in the linear region ////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pAm = new TOTCalib(folder_local+"MAFOutput_MPXNtuple_Am_-300V.root","Am241", minpix, maxpix, 400, nTotalFrames);
    pAm->SetVerboseLevel(verb);
    pAm->SetHistoRebinning(2);
    pAm->SetKernelBandWidth(15); // must be close to the sigma of the peak 
    pAm->OptimizeOnePeak(0.3); // remove all found peaks with amplitude lower than "highest peaks amplitude * thl"
    pAm->Loop(); // Run this source

    // Process the 2nd source in the linear region /////////////////////////////////////////////////////////////////////////////
    TOTCalib * pCd = new TOTCalib(folder_local+"MAFOutput_MPXNtuple_Cd_-300V.root","Cd_fluo_TPXGaAs", minpix, maxpix, 250, nTotalFrames);
    pCd->SetVerboseLevel(verb);
    pCd->SetKernelBandWidth(10); // Kernel Bandwidth
    pCd->SetHistoRebinning(2); 
    pCd->OptimizeOnePeak(0.3);    
    pCd->Loop(); // Run this source    
    
    // Process the low energy source ///////////////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pCu = new TOTCalib(folder_local+"MAFOutput_MPXNtuple_Cu_-300V.root","Cu_fluo", minpix, maxpix, 120, nTotalFrames);
    pCu->SetVerboseLevel(verb);
    pCu->SetHistoRebinning(2);     
    pCu->SetKernelBandWidth(10); // Kernel Bandwidth
    pCu->SetLowEnFit_Params(80.,10.,200.,1.); // constant, sigma, c, t
    pCu->Loop(); // Run this source
    
    // Mix and produce calibration //////////////////////////////////////////////////////////////////////////////////
    // Blend the results and get the surrogate functions.
    // This member function prints out the set of parameters
    // a,b,c,t to separate output files.
    pAm->Choose2LinearPeaks(59.5409,23.1); // in keV
    pAm->SetLowEnergyPeak(8.05);
    //pAm->SetLowEnergyPeak(5.899);    
    pAm->Blender2(pCd,pCu,"test");

    // Partial report ///////////////////////////////////////////////////////////////////////////////////////////////
    // Draw full info for a few pixels.  You can run this interactively in CINT too
    // root [2] pFe55->DrawFullPixelCalib(50051);
//    int reportCntr = 0;
//    for(int i = minpix; i <= maxpix ; i++){
//        pAm->DrawFullPixelCalib2( i );
//        if(reportCntr++ > 2) break; // make only 5 reports
//    }

}
