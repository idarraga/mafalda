#!/bin/bash

# Author: John Idarraga <idarraga@cern.ch>
# Medipix Group, Universite de Montreal

# Starting up a new algorithm in MAFalda

echo
echo "[-->?] Enter the name of your new algorithm"
echo "       avoiding spaces ex: \"AlphaCalculator\""

read algoName

echo "[-->?] Choose between the available templates:"
echo "       1) Simple matrix read and processing - no Clustering.  Includes calibration."
echo "       2) Analysis after Clustering.  Includes calibration."
echo "       3) Example of communication between two algorithms."
echo "       4) Analysis after Clustering & Pattern recognition."
echo "       5) Empty analysis.  Code skeleton only."

example1=1
example2=2
example3=3
example4=4
example5=5

read algoTemplate

if [ $algoTemplate -lt $example1 ]; then
    echo "[ERROR] pick up a number withing the templates range."
    exit
fi
if [ $algoTemplate -gt $example5 ]; then
    echo "[ERROR] pick up a number within the options listed."
    exit
fi

if [ -d ./$algoName ]; then
    echo "[WARNING] directory ./$algoName already exists !!!"
    echo "          hit any key to continue and OVERRIDE, or Ctrl+C to quit"
    read q1
fi

echo "[INFO] creating directory for new algorithm --> ./${algoName}"
if [ ! -d ./$algoName ]; then
    mkdir ./$algoName
fi

# save the user defined name
algoBaseName=$algoName

# Proceed bulding the algo(s)
cd ./$algoName

echo "[INFO] creating header and implementation files"
echo "       ./${algoName}.h and ./${algoName}.cpp"
echo "        using template $algoTemplate ..."

thedate=`date`

##############################################
# template 1
if [ $algoTemplate == $example1 ]; then

cat <<EOF> ./${algoName}.h

/**
 * Created automatically with MAFalda (${thedate})
 * MAFalda Author:  John Idarraga <idarraga@cern.ch>
 *
 * An example of simple processing without clustering.
 */

#ifndef __${algoName}_h
#define __${algoName}_h

#include "MPXAlgo/MediPixAlgo.h"
#include "CalibrationLoader/CalibrationLoader.h"

class ${algoName} : public MediPixAlgo , public CalibrationLoader {

public:

  ${algoName}();
  virtual ~${algoName}() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

private:

  // Variables for output
  // First in the form of a C-style structure
  typedef struct {
    Int_t totalHits;
    Int_t totalCharge;
  } score;
  score m_myhits;
  // In the form of a vector
  vector<double> m_energyPerPixel;

  ClassDef(${algoName}, 1)
};

#endif
EOF


cat <<EOF> ./${algoName}.cpp

/**
 * Created automatically with MAFalda (${thedate})
 * MAFalda Author:  John Idarraga <idarraga@cern.ch>
 *
 * An example of simple processing without clustering.
 */

#ifndef __${algoName}_cpp
#define __${algoName}_cpp

#include "${algoName}.h"

using namespace MSG;

ClassImp(${algoName})

${algoName}::${algoName}() : MediPixAlgo(), CalibrationLoader(this) {

  m_myhits.totalHits = 0;
  m_myhits.totalCharge = 0;

}

void ${algoName}::Init() {

  // You will get an ntuple file containing a TTree with the name of this
  //  Algorithm.  The branches registered through getMyTree() get registered
  //  in that tree so you can fill them each time Execute() gets called.
  getMyTree()->Branch("myhits", &m_myhits , "totalHits/I:totalCharge");
  getMyTree()->Branch("energyPerPixel", &m_energyPerPixel);

  Log << MSG::INFO << "Init function !" << endreq;

}

void ${algoName}::Execute(){

  // REMINDER: ${algoName}::Execute() runs once per frame
  //  you may need ro reinitialize variables.
  m_myhits.totalCharge = 0;
  m_myhits.totalHits = 0;
  
  int xDim = GetMatrixXdim();
  int yDim = GetMatrixYdim();

  int rowItr = 0, colItr = 0, tot = 0;
  double calib_edep = 0.;

  for(colItr = 0; colItr < xDim ; colItr++) {

      for(rowItr = 0; rowItr < yDim ; rowItr++) {

	  if(GetMatrixElement(colItr, rowItr) > 0) {

	      // Get the TOT
	      tot = GetMatrixElement(colItr, rowItr);

	      // Get the TOT or counts
	      m_myhits.totalCharge += tot;
	      m_myhits.totalHits++;

	      // Get the calibrated energy if the calibration files are loaded in the 
              //  top layer (ROOT Macro used to drive the run).
	      calib_edep = CalculateAndGetCalibEnergy(colItr, rowItr, tot);

  	      // Push back in the vector for output.  There is no notion of "cluster"
	      //  in this example.  We are simply storing the energy from all pixels.
	      m_energyPerPixel.push_back( calib_edep );

	    }

	}
    }

  Log << MSG::INFO << "Total charge/counts     : " << m_myhits.totalCharge << endreq ;
  Log << MSG::INFO << "Number of active pixels : " << m_myhits.totalHits << endreq ;

  getMyTree()->Fill();

  // Get output variables ready for next frame after Filling the Tree
  m_energyPerPixel.clear();

}

void ${algoName}::Finalize(){

  Log << MSG::INFO << "Finalize function !" << endreq;

}

#endif
EOF

##############################################
# template 2
elif [ $algoTemplate == $example2 ]; then

cat <<EOF> ./${algoName}.h

/**
 * Created automatically with MAFalda (${thedate})
 * MAFalda Author:  John Idarraga <idarraga@cern.ch>
 *
 * An example of how to use the clustering results
 * for further processing.
 */

#ifndef __${algoName}_h
#define __${algoName}_h

#include "MPXAlgo/MediPixAlgo.h"
#include "CalibrationLoader/CalibrationLoader.h"

// This algorithm is using an object put in the StoreGate
//  by BlobsFinder. I need to know the object.
#include "BlobsFinder/BlobsFinder.h"

class ${algoName} : public MediPixAlgo , public CalibrationLoader {

public:

  ${algoName}();
  virtual ~${algoName}() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

private:

  AllBlobsContainer * m_aB;
  Int_t m_minNPixels;

  // for output
  vector<int> m_clusterTOT;
  vector<double> m_clusterEnergy;

  ClassDef(${algoName}, 1)
};

#endif
EOF

cat <<EOF> ./${algoName}.cpp

/**
 * Created automatically with MAFalda (${thedate})
 * MAFalda Author:  John Idarraga <idarraga@cern.ch>
 *
 * An example of how to use the clustering results
 * for further processing.
 */

#ifndef __${algoName}_cpp
#define __${algoName}_cpp

#include "${algoName}.h"
#include "MAFTools.h"

using namespace MSG;

ClassImp(${algoName})

${algoName}::${algoName}() : MediPixAlgo(), CalibrationLoader(this) {

  // This value will be overridden by the configuration because it'll registered
  //  as a ConfigurationValue in the Init member of this class.
  m_minNPixels = 5;

}

void ${algoName}::Init(){

	Log << MSG::INFO << "Init function !" << endreq;


	// You will get an ntuple file containing a TTree with the name of this
	//  Algorithm.  The branches registered through getMyTree() get registered
	//  in that tree so you can fill them each time Execute() gets called.
	getMyTree()->Branch("clusterTOT", &m_clusterTOT);
	getMyTree()->Branch("clusterEnergy", &m_clusterEnergy);

	// A configuration value that can be tuned from the Viewer
	RegisterConfigurationValue(&m_minNPixels, "minNPixels");

}

void ${algoName}::Execute(){

	// Ask the store gate if the previous algorithm (BlobsFinder --> reponsible for clustering)
	//  sent any objects to the StoreGate.
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");
	if(lastObject == 0)
		return;

	// If so, get the pointer to the last object.  BlobsFinder, as it comes out of the box, sends
	//  a single object containing all clusters.
	m_aB = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject-1);

	// AllBlobsContainer is a box full of Blobs(Clusters). Now we can iterate over all clusters inside.
	vector<blob> blobsVector = m_aB->GetBlobsVector();
	Log << MSG::INFO << "Number of blobs from clustering = " << (Int_t) blobsVector.size() << endreq;
	vector<blob>::iterator blobsItr = blobsVector.begin(); //allBlobs.begin();

	cluster cl;
	for( ; blobsItr != blobsVector.end() ; blobsItr++) {

		cl = *blobsItr;

		// Limit all this to clusters with a minimum size.
		// Note that m_minNPixels can be configured through the Viewer
		//  so you can reprocess and check results online.
		if(cl.bP.nPixels < m_minNPixels)
			continue;
		
		// If the cluster passes the minimum requirements loop over the
		// the constituents of the cluster --> ClusterDescription
		list< pair < pair<int, int>, int > > cl_des = cl.GetClusterDescription();
		list< pair < pair<int, int>, int > >::iterator i = cl_des.begin();

		// Store the cluster TOT for output
		m_clusterTOT.push_back( cl.bP.clusterTOT );

		double calib_edep = 0.0, clusterEdep = 0.;
		int tot = 0;
		pair<int, int> pix;

		for ( ; i != cl_des.end() ; i++) {

			// pixel coordinates and tot
			pix = (*i).first;
			tot = (*i).second;

			// Use calibration to obtain E = Surrogate(TOT) for this pixel
			calib_edep = CalculateAndGetCalibEnergy(pix, tot);

			// Calculate the energy of the cluster
			clusterEdep += calib_edep;

		}

		// Store the cluster Energy calculated in the previous loop
		m_clusterEnergy.push_back( clusterEdep );

	}

	// Fill the output tree of this algorithm
	getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_clusterEnergy.clear();
	m_clusterTOT.clear();

}

void ${algoName}::Finalize() {

	Log << MSG::INFO << "Finalize function !" << endreq;

}

#endif
EOF

##############################################
# template 3
elif [ $algoTemplate == $example3 ]; then

    source ../share/StoreGateExample.sh

##############################################
# template 4
elif [ $algoTemplate == $example4 ]; then

    echo "not available yet.  Contact <idarraga@cern.ch> if you need this."
    exit

##############################################
# template 5
elif [ $algoTemplate == $example5 ]; then

    echo "not available yet.  Contact <idarraga@cern.ch> if you need this."
    exit

fi

##############################################
# Creating a top layer
cd ..
source ./share/TopLayerExample.sh

##############################################
# Preparing to build

# Special case of two algorithms
if [ $algoTemplate == $example3 ]; then

echo "[INFO] appending \"${algoA}\" and \"${algoB}\" to buildList.Algo and buildList.LinkDef"
echo "$algoBaseName" >> buildList.Algo
echo "$algoA" >> buildList.LinkDef
echo "$algoB" >> buildList.LinkDef

else

echo "[INFO] appending \"${algoBaseName}\" to buildList.Algo and buildList.LinkDef"
echo "$algoBaseName" >> buildList.Algo
echo "$algoBaseName" >> buildList.LinkDef

fi

##############################################
# Report to developper

echo ""
echo "[DONE] README !!! The new algorithm \"${algoName}\" is ready. "
echo "       You need to do the following: "
echo "       1) To compile, type --> make"
echo "       2) To run a job run the ROOT-Macro ./run${algoBaseName}.C like this"
echo "          root -l ./run${algoBaseName}.C"
echo "       4) If you want to run without the Viewer, comment out the line"
echo "          \"mpxAnalysis.ConnectAlgo("MPXViewer", v1);\" and run again"
echo "          You will get the results of the running algorithms in an output"
echo "          file called MAFOutput_?.root"

