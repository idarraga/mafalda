/*
 * 	Copyright 2014 John Idarraga
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

void runTOTCalibrationPreparation(TString fn){

  /**
   * Author: John Idarraga <idarraga@cern.ch> 
   * LAL, Orsay.
   *
   * C++ ROOT macro for analysis of MediPix data using multiformat
   * ROOT-based files and MAFalda - the MediPix Analysis Framework
   * 
   */

  /* Load the MediPix analysis lib */
  gSystem->Load("libMediPixAnalysisCore.so");
  
  /* Get an instance of the analysis manager */
  AnalysisManager mpxAnalysis( fn );

  mpxAnalysis.SetStoreGateOutputLevel(MSG::ERROR);

  //string outfile = "output.root";

  /* Blobs Finder */
  SimpleClustering * ac = new SimpleClustering;
  mpxAnalysis.ConnectAlgo("SimpleClustering", ac);
  //ac->changeOutputLevel(MSG::DEBUG);

  //BlobsFinder * bf = new BlobsFinder;
  //mpxAnalysis.ConnectAlgo("BlobsFinder", bf);
  //bf->changeOutputLevel(MSG::ERROR);
  //bf->SetBorderExclusion(0);
  //bf->SetDiscontinuityTolerance(0);
  //bf->ReadConfiguration();

  TOTCalibrationPreparation * cp = new TOTCalibrationPreparation;
  mpxAnalysis.ConnectAlgo("TOTCalibrationPreparation", cp);
  //cp->changeOutputLevel(MSG::DEBUG);

  /* This is an special algorithm that works as a
     frames viewer */
  MPXViewer * v1 = new MPXViewer;
  //v1->changeOutputLevel(MSG::DEBUG);
  //mpxAnalysis.ConnectAlgo("MPXViewer", v1);
  v1->SetCuts(0, 0);
  v1->SetFrameTitle("calib preparation");

  /* Get a list of the algorithms connected */
  mpxAnalysis.DumpAlgoList();

  /* If you want to setup an output filename yourself 
   * otherwise MAFalda builds something convenient
   */
  //mpxAnalysis.SetOutputNtupleFilename(outfile.c_str());
  
  /* run ! */
  mpxAnalysis.Run();         // all frames
  //mpxAnalysis.Run(0, 100); // range of frames
  //mpxAnalysis.Run(1);      // ony one frame
  
}
