/*
 * 	Copyright 2012 John Idarraga
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
 * MAFalda Top layer file in the form of a C++ ROOT macro
 *
 * - Run this script in the following way
 *   $ root -l runExample.C
 * - For a full debug run look at share/runCintAttachGdbEmacs.C
 */

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
using namespace std;

R__LOAD_LIBRARY(libMediPixAnalysisCore) // ROOT6 (cling)

void runExample(TString fn = "") {

	//gSystem->Load("libMediPixAnalysisCore"); // ROOT5 (CINT)

	// Get an instance of the analysis manager, including input data
	// hint: it can also be a text file containing a list of input files
	if ( fn == "" ) fn = "./testdata/MPXNtuple_12C_TimePix.root";
	AnalysisManager mpxAnalysis( fn.Data() );

	// A very simple example algorithm
	AnalysisExample * e1 = new AnalysisExample;
	mpxAnalysis.ConnectAlgo("SimpleExample", e1);
	//e1->UseExtendedMaskMatrixFile("/home/idarraga/storage/maskbits.txt");

	/*
	 * Clustering:
	 *  BlobsFinder is a clustering algorithm.
	 *  It supports discontinuities and other
	 *  exceptions.  It is suitable for most
	 *  applications.
	 */
	BlobsFinder * bf = new BlobsFinder;
	mpxAnalysis.ConnectAlgo("BlobsFinder", bf);
	bf->changeOutputLevel(MSG::INFO);
	bf->SetBorderExclusion(0);
	bf->SetDiscontinuityTolerance(1);
	//bf->ReadConfiguration();

	PRBasicSpecies * pbs = new PRBasicSpecies;
	mpxAnalysis.ConnectAlgo("PRBasicSpecies", pbs);
	

	// This is an special algorithm which works as a frames viewer
	MPXViewer * v1 = new MPXViewer;
	v1->changeOutputLevel(MSG::DEBUG);
	mpxAnalysis.ConnectAlgo("MPXViewer", v1);
	//v1->SetCuts(10, 10);
	v1->SetFrameTitle("add a title");

	// Get a list of the algorithms connected
	mpxAnalysis.DumpAlgoList();

	// Use SetOutputNtuple() If you want to change the output filename.
	// Otherwise MAFalda will do it based on the input filename above
	//mpxAnalysis.SetOutputNtuple("myoutput.root");

	// run !
	mpxAnalysis.Run();         // all frames
	//mpxAnalysis.Run(0, 100); // range of frames
	//mpxAnalysis.Run(1);      // single frame

}
