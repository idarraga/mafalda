echo "[INFO] creating a top layer file to drive the run (ROOT macro)"
echo "       ./run${algoBaseName}.C"

# Special case of two algorithms
if [ $algoTemplate == $example3 ]; then

cat <<EOF> ./run${algoBaseName}.C

/**
 * Created automatically with MAFalda (${thedate})
 * MAFalda Author:  John Idarraga <idarraga@cern.ch>
 *
 * An example of a top layer file used to steer the job.
 * To run in a shell issue the command:
 *  $ root -l run${algoBaseName}.C
 */

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
using namespace std;

R__LOAD_LIBRARY(libMediPixAnalysisCore) // ROOT6 (cling)

void run${algoBaseName}(TString fn = "testdata/MPXNtuple_12C_TimePix.root"){

	// Load the MediPix analysis lib
	gSystem->Load("libMediPixAnalysisCore.so");

	// Get an instance of the analysis manager and load data
	AnalysisManager mpxAnalysis( fn );

	// Blobs Finder.  This is a Clustering Algorithm
	BlobsFinder * bf = new BlobsFinder;
	mpxAnalysis.ConnectAlgo("BlobsFinder", bf);
	bf->changeOutputLevel(MSG::INFO);
	bf->SetBorderExclusion(1);
	bf->SetDiscontinuityTolerance(0);
	//bf->ReadConfiguration();

	// Your First algorithm --> ${algoA}
	${algoA} * ac_A = new ${algoA};
	mpxAnalysis.ConnectAlgo("${algoA}", ac_A);
	// Setting up the calibration from this algorithm by loading
	//  the 4 calibration files and the clock used for calibration. 
	//  If you fail loading any of the 5 following pieces then the 
	//  calibration won't be available but you can still run the job.
	//ac_A->SetCalibrationConfigFile_a("/__set_your_own_path__/calib_a.txt");
	//ac_A->SetCalibrationConfigFile_b("/__set_your_own_path__/calib_b.txt");
	//ac_A->SetCalibrationConfigFile_c("/__set_your_own_path__/calib_c.txt");
	//ac_A->SetCalibrationConfigFile_t("/__set_your_own_path__/calib_t.txt");
	//ac_A->SetCalibClk(9.6); // MHz
	//ac_A->ReadConfiguration();
	//ac_A->OutputLevel(MSG::DEBUG);

	// Your Second algorithm --> ${algoB}
	${algoB} * ac_B = new ${algoB};
	mpxAnalysis.ConnectAlgo("${algoB}", ac_B);

	// This is an special algorithm that works as a frames viewer
	MPXViewer * v1 = new MPXViewer;
	v1->changeOutputLevel(MSG::DEBUG);
	mpxAnalysis.ConnectAlgo("MPXViewer", v1);
	//v1->SetCuts(10,10); // A minimum cut to skip uninteresting frames
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

EOF

else

cat <<EOF> ./run${algoBaseName}.C

/**
 * Created automatically with MAFalda (${thedate})
 * MAFalda Author:  John Idarraga <idarraga@cern.ch>
 *
 * An example of a top layer file used to steer the job.
 * To run in a shell issue the command:
 *  $ root -l run${algoBaseName}.C
 */

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
using namespace std;

R__LOAD_LIBRARY(libMediPixAnalysisCore) // ROOT6 (cling)

void run${algoBaseName}(TString fn = "testdata/MPXNtuple_12C_TimePix.root"){

	// Load the MediPix analysis lib
	gSystem->Load("libMediPixAnalysisCore.so");

	// Get an instance of the analysis manager and load data
	AnalysisManager mpxAnalysis( fn );

	// Blobs Finder.  This is a Clustering Algorithm
	BlobsFinder * bf = new BlobsFinder;
	mpxAnalysis.ConnectAlgo("BlobsFinder", bf);
	bf->changeOutputLevel(MSG::INFO);
	bf->SetBorderExclusion(1);
	bf->SetDiscontinuityTolerance(0);
	//bf->ReadConfiguration();

	// Your algorithm --> ${algoBaseName}
	${algoBaseName} * ac = new ${algoBaseName};
	mpxAnalysis.ConnectAlgo("${algoBaseName}", ac);
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

	// This is an special algorithm that works as a frames viewer
	MPXViewer * v1 = new MPXViewer;
	v1->changeOutputLevel(MSG::DEBUG);
	mpxAnalysis.ConnectAlgo("MPXViewer", v1);
	//v1->SetCuts(10,10); // A minimum cut to skip uninteresting frames
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

EOF

fi
