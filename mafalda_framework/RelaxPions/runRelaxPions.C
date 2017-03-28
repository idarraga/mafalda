/*
 * 	Copyright 2012 John Idarraga, Mathieu Benoit
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

/**
 * Created automatically with MAFalda (Fri Nov 23 19:24:07 CET 2012)
 *
 * An example of a top layer file used to steer the job.
 * To run in a shell issue the command:
 *  $ root -l runRelaxPions.C
 */

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
using namespace std;

R__LOAD_LIBRARY(libMediPixAnalysisCore) // ROOT6 (cling)

void runRelaxPions(){

	// Load the MediPix analysis lib
	gSystem->Load("libMediPixAnalysisCore.so");

	// Get an instance of the analysis manager and load data
	//AnalysisManager mpxAnalysis("testdata/MPXNtuple_12C_TimePix.root");
	//AnalysisManager mpxAnalysis("/home/idarraga/storage/NIKHEF/MPXNtuple_2012-11-23_set2_-55V_run001.root");
	AnalysisManager mpxAnalysis("/media/TOSHIBA EXT/Testbeam_Nov_2012/Nov2012_edgelessDetectors_Nikhef/MPXNtuple_Nov2012_edgelessDetectors_Nikhef_2012-11-23_set2_-55V_run008.root");

	// Blobs Finder.  This is a Clustering Algorithm
	BlobsFinder * bf = new BlobsFinder;
	mpxAnalysis.ConnectAlgo("BlobsFinder", bf);
	bf->changeOutputLevel(MSG::INFO);
	bf->SetBorderExclusion(1);
	bf->SetDiscontinuityTolerance(8);
	//bf->ReadConfiguration();

	PRBasicSpecies * pr = new PRBasicSpecies;
	mpxAnalysis.ConnectAlgo("PRBasicSpecies", pr);

	// Your algorithm --> RelaxPions
	RelaxPions * ac = new RelaxPions;
	mpxAnalysis.ConnectAlgo("RelaxPions", ac);
	// Setting up the calibration from this algorithm by loading
	//  the 4 calibration files. If you fail loading any of the 4 calibration
	//  files the calibration won't be available but you can still run the job.
	//ac->SetCalibrationConfigFile_a("/__set_your_own_path__/calib_a.txt");
	//ac->SetCalibrationConfigFile_b("/__set_your_own_path__/calib_b.txt");
	//ac->SetCalibrationConfigFile_c("/__set_your_own_path__/calib_c.txt");
	//ac->SetCalibrationConfigFile_t("/__set_your_own_path__/calib_t.txt");
	//ac->ReadConfiguration();
	//ac->OutputLevel(MSG::DEBUG);

	// This is an special algorithm that works as a frames viewer
	MPXViewer * v1 = new MPXViewer;
	v1->changeOutputLevel(MSG::DEBUG);
	mpxAnalysis.ConnectAlgo("MPXViewer", v1);
	v1->SetCuts(50,50); // A minimum cut to skip uninteresting frames
	v1->SetFrameTitle("Offline Analysis");

	// Get a list of the algorithms connected
	mpxAnalysis.DumpAlgoList();

	// If you want to setup an output filename yourself
	//  otherwise MAFalda builds something convenient
	//mpxAnalysis.SetOutputNtupleFilename("outputfile.root");

	// Run !
	mpxAnalysis.Run();         // all frames
	//mpxAnalysis.Run(0, 100); // range of frames
	//mpxAnalysis.Run(1);      // ony one frame

}

