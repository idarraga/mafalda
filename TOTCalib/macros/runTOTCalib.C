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

void runTOTCalib (  ) {

    nTotalFrames = -1;            	 // -1: run all frames in input dataset
    //minpix = 0; maxpix = 256*256-1;  // Work on this set of pixels
    minpix = 32741; maxpix = 32760;  // Work on this set of pixels

    // Load calibration library
    gSystem->Load("libTOTCalib.so");

    // Process one source ///////////////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pFe55 = new TOTCalib("/home/idarraga/analysis/mafalda/mafalda_framework/MAFOutput_MPXNtuple_Timepix3UofH_Calib_Fe55_set9.root",
                                    "Fe55", minpix, maxpix, 20, nTotalFrames);
    //pFe55->SetVerboseLevel(TOTCalib::__VER_DEBUG);
    pFe55->SetKernelBandWidth(2); // Kernel Bandwidth
    pFe55->Loop(); // Run this source
    pFe55->GetInputStats();
    pFe55->SetGlobalThresholdEnergyAndErr(4.14, 0.414);

    // Process one source ///////////////////////////////////////////////////////////////////////////////////////////
    TOTCalib * pAm241PlusSn = new TOTCalib("/home/idarraga/analysis/mafalda/mafalda_framework/MAFOutput_MPXNtuple_Timepix3UofH_Calib_Am241PlusSn_set9.root",
                                     "Am241CutLow_Plus_SnFluo", minpix, maxpix, 50, nTotalFrames);
    //pAm241->SetVerboseLevel(TOTCalib::__VER_DEBUG);
    pAm241PlusSn->SetKernelBandWidth(3); // Kernel Bandwidth
    pAm241PlusSn->Loop(); // Run this source
    pAm241PlusSn->GetInputStats();

    // Mix and produce calibration //////////////////////////////////////////////////////////////////////////////////
    // Blend the results and get the surrogate functions.
    // This member function prints out the set of parameters
    // a,b,c,t to separate output files.
    pFe55->Blender(pAm241PlusSn, "Timepix3");

    // Partial report ///////////////////////////////////////////////////////////////////////////////////////////////
    // Draw full info for a few pixels.  You can run this interactively in CINT too
    // root [2] pFe55->DrawFullPixelCalib(50051);
    int reportCntr = 0;
    for(int i = minpix; i < maxpix ; i++){
        pFe55->DrawFullPixelCalib( i );
        if(reportCntr++ > 5) break; // make only 5 reports
    }

    // Use the blender object to write a root file with fit information
    //pFe55->Finalize();

}
