/**
 *  Author: John Idarraga <idarraga@cern.ch>
 *  TOT calibration for a Timepix/Timepix3 device
 */

// Tip: execute macro using -b (batch mode) option:
// $ root -l -q -b runTOTCalib.C)
// to avoid problems with X11 when running with tmux from server. (else
// progam crashes at the end if connection was lost and root file is not
// written correctly).

#include <TH2I.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include "TOTCalib.h"


#include <vector>

using namespace std;

class TOTCalib;

// Prototypes and a few harmless globals
int nTotalFrames;
int minpix; int maxpix;

R__LOAD_LIBRARY(libTOTCalib.so); // ROOT6 (cling)

void runTOTCalib (  ) {

    nTotalFrames = -1;            	 // -1: run all frames in input dataset
    //minpix = 0; maxpix = 256*256-1;  // Work on this set of pixels
    minpix = 3000; maxpix = 3005;  // Work on this set of pixels

    // Load calibration library
    //gSystem->Load("libTOTCalib.so");

    TString server = "/lcg/storage13/atlas/medipix/DATA/2017/2017_02_PRAGUE_GaAs500_E06W0203/";
    TString local = " /Users/thomasbilloud/workspace/DATA/MAFOutput/2017-02-GaAs500-E06W0203/";

    // Process one source ///////////////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pCd = new TOTCalib(local+"MAFOutput_MPXNtuple_Cd_-300V.root","Cd-fluo",minpix, maxpix, 180, nTotalFrames);
    //pCd->SetVerboseLevel(TOTCalib::__VER_DEBUG_LOOP);
    pCd->SetVerboseLevel(TOTCalib::__VER_INFO);
    //pCd->SetVerboseLevel(TOTCalib::__VER_DEBUG);
    //pCd->SetVerboseLevel(TOTCalib::__VER_QUIET);
    pCd->SetKernelBandWidth(10.); // Kernel Bandwidth - Should be as close as possible to gaussian sigma
    pCd->Loop(); // Run this source
//    pCd->SavePixelResolution();
//    //pCd->GetInputStats();
//    //pCd->SetGlobalThresholdEnergyAndErr(3.93,1.0);

    // Process one source ///////////////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pAm = new TOTCalib(local+"MAFOutput_MPXNtuple_Am_-300V.root","Am241", minpix, maxpix, 350, nTotalFrames);
    pAm->SetKernelBandWidth(20.); // Kernel Bandwidth
    pAm->Loop(); // Run this source
    //pAm->GetInputStats();

    // Process one source ///////////////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pFe = new TOTCalib(local+"MAFOutput_MPXNtuple_Fe_-300V.root","Fe55", minpix, maxpix, 100, nTotalFrames);
    pFe->SetKernelBandWidth(10.); // Kernel Bandwidth
    pFe->Loop(); // Run this source


    // Mix and produce calibration //////////////////////////////////////////////////////////////////////////////////
    // Blend the results and get the surrogate functions.
    // This member function prints out the set of parameters
    // a,b,c,t to separate output files.
    pCd->Blender(pAm,pFe,"TestGaAs500", TOTCalib::__jakubek);

//    // Partial report ///////////////////////////////////////////////////////////////////////////////////////////////
//    // Draw full info for a few pixels.  You can run this interactively in CINT too
      //pCd->DrawFullPixelCalib_2(3000);
//    //pCd->DrawFullPixelCalib(32746);
    int reportCntr = 0;
    for(int i = minpix; i <= maxpix ; i++){
        pCd->DrawFullPixelCalib( i );
        if(reportCntr++ > 5) break; // make only 5 reports
    }

    // Use the blender object to write a root file with fit information
    pCd->Finalize();

}

