 /*
 * 	Copyright 2011 John Idarraga, Mathieu Benoit
 *
 * 	This file is part of MAFalda.

    MAFalda is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    MAFalda is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MAFalda.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
using namespace std;

R__LOAD_LIBRARY(libMediPixAnalysisCore) // ROOT6 (cling)

void runFEI3() {

  /* Load the MediPix analysis lib */
  gSystem->Load("libMediPixAnalysisCore.so");
  
  /* Get an instance of the analysis manager */
  //AnalysisManager mpxAnalysis("testdata/MPXNtuple_241Am_106Ru_137Cs_13keV_30deg_784_100V.root");
  //AnalysisManager mpxAnalysis("/home/idarraga/geant4/allpix/MPXNtuple_AllPix_det_104.root");
  //AnalysisManager mpxAnalysis("/home/idarraga/geant4/allpix/positron4GeV_test3/MPXNtuple_AllPix_det_102.root");
  //AnalysisManager mpxAnalysis("/home/idarraga/analysis/pixeldataconverter/MPXNtuple_USBPixRaw.root");

  //AnalysisManager mpxAnalysis("/home/idarraga/analysis/mpxdataconverter/MPXNtuple_USBPixRaw_LALCosmics.root");
  //AnalysisManager mpxAnalysis("testdata/MPXNtuple_12C_TimePix.root");

  //AnalysisManager mpxAnalysis("testdata/MPXNtuple_USBPixRaw_FEI3_Sr90.root");
  AnalysisManager mpxAnalysis("/home/idarraga/analysis/mpxdataconverter/MPXNtuple_USBPixRaw.root");

  //string outfile = "output.root";

  /* Meta Data */
  MetaData * md = new MetaData;
  mpxAnalysis.ConnectAlgo("MetaData", md);

  /* Blobs Finder */
  BlobsFinder * bf = new BlobsFinder;
  mpxAnalysis.ConnectAlgo("BlobsFinder", bf);
  bf->changeOutputLevel(MSG::INFO);
  bf->SetBorderExclusion(0);
  bf->SetDiscontinuityTolerance(2);
  bf->ReadConfiguration();

  /* Patter Recognition basic species */
  PRBasicSpecies * pr = new PRBasicSpecies;
  mpxAnalysis.ConnectAlgo("PRBasicSpecies", pr);
  pr->changeOutputLevel(MSG::INFO);
  pr->SetHeavyBlobOnOff(true);
  pr->ReadConfiguration();

  FEI3Mips * fe = new FEI3Mips;
  mpxAnalysis.ConnectAlgo("FEI3Mips", fe);
  fe->changeOutputLevel(MSG::INFO);
  fe->ReadConfiguration();

  /* This is an special algorithm that works as a
     frames viewer */
  MPXViewer * v1 = new MPXViewer;
  v1->changeOutputLevel(MSG::DEBUG);
  mpxAnalysis.ConnectAlgo("MPXViewer", v1);
  v1->SetCuts(10, 10);
  v1->SetFrameTitle("FE-I3 @ DESY | back side tilt");

  /* Get a list of the algorithms connected */
  mpxAnalysis.DumpAlgoList();

  /* If you want to setup an output filename yourself 
   * otherwise MAFalda builds something convenient
   */
  //mpxAnalysis.SetOutputNtupleFilename(outfile.c_str());
  
  /* run ! */
  mpxAnalysis.Run(0, 20000);         // all frames
  //mpxAnalysis.Run(0, 100); // range of frames
  //mpxAnalysis.Run(1);      // ony one frame
  
}
