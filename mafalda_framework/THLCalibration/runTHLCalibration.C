 /*
 * 	Copyright 2013 John Idarraga
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
 * Created automatically with MAFalda (Thu Apr  4 10:01:28 CEST 2013)
 *
 * An example of a top layer file used to steer the job.
 * To run in a shell issue the command:
 *  $ root -l runTHLCalibration.C
 */

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
using namespace std;

R__LOAD_LIBRARY(libMediPixAnalysisCore) // ROOT6 (cling)

void runTHLCalibration (TString fn) {

	// Load the MediPix analysis lib
	gSystem->Load("libMediPixAnalysisCore.so");

	// Get an instance of the analysis manager and load data
	AnalysisManager mpxAnalysis( fn.Data() );

	// Blobs Finder.  This is a Clustering Algorithm
	SimpleClustering * ac = new SimpleClustering;
	//mpxAnalysis.ConnectAlgo("SimpleClustering", ac);
	//ac->changeOutputLevel(MSG::DEBUG);

	BlobsFinder * ac2 = new BlobsFinder;
	mpxAnalysis.ConnectAlgo("BlobsFinder", ac2);
	ac2->SetBorderExclusion(0);
	//ac2->SetMaxOcc(0.5);
	ac2->changeOutputLevel(MSG::FATAL);


	// Your algorithm --> THLCalibration
	THLCalibration * thlc = new THLCalibration;
	mpxAnalysis.ConnectAlgo("THLCalibration", thlc);
	//thlc->UseExtendedMaskFile("/media/datum/NASA_Assemblies_March2013/K09_W0214/K09_W0214_THLCalibration_ExtendedMask.txt");
	//thlc->UseExtendedMaskFile("/media/datum/NASA_Assemblies_March2013/K09_W0214/K09_W0214_THLCalibration_Freq_ExtendedMask.txt");

	// Setting up the calibration from this algorithm by loading
	//  the 4 calibration files and the clock used for calibration. 
	//  If you fail loading any of the 5 following pieces then the 
	//  calibration won't be available but you can still run the job.
	//ac->SetCalibrationConfigFile_a("/__set_your_own_path__/calib_a.txt");
	//ac->SetCalibrationConfigFile_b("/__set_your_own_path__/calib_b.txt");
	//ac->SetCalibrationConfigFile_c("/__set_your_own_path__/calib_c.txt");
	//ac->SetCalibrationConfigFile_t("/__set_your_own_path__/calib_t.txt");
	//ac->SetCalibClk(9.6); // MHz
	//ac->ReadConfiguration();
	//ac->OutputLevel(MSG::DEBUG);

	////////////////////////////////////////////////////////////////////
	// Two algos to get the exteded mask
	// Blobs Finder
	//BlobsFinder * bf = new BlobsFinder;
	//mpxAnalysis.ConnectAlgo("BlobsFinder", bf);
	//bf->changeOutputLevel(MSG::INFO);
	//bf->SetBorderExclusion(0);
	//bf->SetDiscontinuityTolerance(0);
	//CalibrationValidation * cv = new CalibrationValidation;
	//mpxAnalysis.ConnectAlgo("CalibrationValidation", cv);
	////////////////////////////////////////////////////////////////////

	// This is an special algorithm that works as a frames viewer
	MPXViewer * v1 = new MPXViewer;
	v1->changeOutputLevel(MSG::DEBUG);
	mpxAnalysis.ConnectAlgo("MPXViewer", v1);
	//v1->SetCuts(10,10); // A minimum cut to skip uninteresting frames
	v1->SetFrameTitle("THL calibration");

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

