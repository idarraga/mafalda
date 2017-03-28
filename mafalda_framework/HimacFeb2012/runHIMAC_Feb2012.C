/*
 * 	Copyright 2011 John Idarraga
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

void runHIMAC_Feb2012(TString fn, TString detector, TString config){

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
	AnalysisManager mpxAnalysis(fn);
	//string outfile = "output.root";

	// Meta Data
	MetaData * md = new MetaData;
	mpxAnalysis.ConnectAlgo("MetaData", md);

	// Blobs Finder
	BlobsFinder * bf = new BlobsFinder;
	mpxAnalysis.ConnectAlgo("BlobsFinder", bf);
	bf->changeOutputLevel(MSG::INFO);
	bf->SetBorderExclusion(1);
	bf->SetDiscontinuityTolerance(0);
	//bf->ReadConfiguration();

	// Patter Recognition basic species
	//PRBasicSpecies * pr = new PRBasicSpecies;
	//mpxAnalysis.ConnectAlgo("PRBasicSpecies", pr);
	//pr->changeOutputLevel(MSG::DEBUG);
	//pr->SetHeavyBlobOnOff(true);
	//pr->ReadConfiguration();

	// HIMAC Feb2012
	HimacFeb2012 * ac = new HimacFeb2012;
	mpxAnalysis.ConnectAlgo("HimacFeb2012", ac);
	ac->ReadConfiguration(config);
	//ac->OutputLevel(MSG::DEBUG);
	if (detector == "G03") {
		// G03 - flight unit
		cout << "[INFO] (top macro) loading G03 calibration files" << endl;
		ac->SetCalibrationConfigFile_a("/home/idarraga/storage/HIMAC-Feb2012/calib/G03_J03Cal/calib_G03-W0094_a.txt");
		ac->SetCalibrationConfigFile_b("/home/idarraga/storage/HIMAC-Feb2012/calib/G03_J03Cal/calib_G03-W0094_b.txt");
		ac->SetCalibrationConfigFile_c("/home/idarraga/storage/HIMAC-Feb2012/calib/G03_J03Cal/calib_G03-W0094_c.txt");
		ac->SetCalibrationConfigFile_t("/home/idarraga/storage/HIMAC-Feb2012/calib/G03_J03Cal/calib_G03-W0094_t.txt");

	} else if (detector == "J03") {
		// J03 - flight unit
		cout << "[INFO] (top macro) loading J03 calibration files" << endl;
		ac->SetCalibrationConfigFile_a("/home/idarraga/storage/HIMAC-Feb2012/calib/G03_J03Cal/calib_J03-W0094_a.txt");
		ac->SetCalibrationConfigFile_b("/home/idarraga/storage/HIMAC-Feb2012/calib/G03_J03Cal/calib_J03-W0094_b.txt");
		ac->SetCalibrationConfigFile_c("/home/idarraga/storage/HIMAC-Feb2012/calib/G03_J03Cal/calib_J03-W0094_c.txt");
		ac->SetCalibrationConfigFile_t("/home/idarraga/storage/HIMAC-Feb2012/calib/G03_J03Cal/calib_J03-W0094_t.txt");
	} else if (detector == "G04") {
		// G4 - Lawrence' detector
		cout << "[INFO] (top macro) loading G04 calibration files" << endl;
		ac->SetCalibrationConfigFile_a("/home/idarraga/storage/HIMAC-Feb2012/calib/Box_053/Calibration/FBK_128/E_not_shifted/Calib_2point_a.txt");
		ac->SetCalibrationConfigFile_b("/home/idarraga/storage/HIMAC-Feb2012/calib/Box_053/Calibration/FBK_128/E_not_shifted/Calib_2point_b.txt");
		ac->SetCalibrationConfigFile_c("/home/idarraga/storage/HIMAC-Feb2012/calib/Box_053/Calibration/FBK_128/E_not_shifted/Calib_2point_c.txt");
		ac->SetCalibrationConfigFile_t("/home/idarraga/storage/HIMAC-Feb2012/calib/Box_053/Calibration/FBK_128/E_not_shifted/Calib_2point_t.txt");
		ac->SetCalibClkFactor(20.0/9.6); // calib @ 20MHz, data at 9.6MHz
	} else {
		cout << "[WARNING] working without calibration " << endl;
	}

	ac->SetMaxFrameWidth(256);
	ac->SetMaxFrameHeight(256);

	/* This is an special algorithm that works as a
     frames viewer */
	MPXViewer * v1 = new MPXViewer;
	v1->changeOutputLevel(MSG::DEBUG);
	mpxAnalysis.ConnectAlgo("MPXViewer", v1);
	//v1->SetCuts(10,10);
	v1->SetFrameTitle("Offline Analysis");


	/* Get a list of the algorithms connected */
	mpxAnalysis.DumpAlgoList();

	/* If you want to setup an output filename yourself
	 * otherwise MAFalda builds something convenient */
	//mpxAnalysis.SetOutputNtupleFilename(outfile.c_str());

	/* run ! */
	mpxAnalysis.Run();         // all frames
	//mpxAnalysis.Run(0, 100); // range of frames
	//mpxAnalysis.Run(1);      // ony one frame

}
