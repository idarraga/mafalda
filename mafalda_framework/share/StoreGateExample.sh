################################################################
# StoreGate communication example
userAlgoName=${algoName}

################################################################
# Algorithm A 
algoName=${userAlgoName}_A
algoA=$algoName

cat <<EOF> ./${algoName}.h

/**
 * Created automatically with MAFalda (${thedate})
 * MAFalda Author:  John Idarraga <idarraga@cern.ch>
 *
 * An example on how to implement the communication
 * between two algorithms using the StoreGate.
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
 * An example on how to implement the communication
 * between two algorithms using the StoreGate.
 */

#ifndef __${algoName}_cpp
#define __${algoName}_cpp

#include "${algoName}.h"
#include "MAFTools.h"

using namespace MSG;

ClassImp(${algoName})

${algoName}::${algoName}() : MediPixAlgo() , CalibrationLoader( this ) {

  // This value will be overridden by the configuration because it'll registered
  //  as a ConfigurationValue in the Init member of this class.
  m_minNPixels = 5;

}

void ${algoName}::Init() {

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

################################################################
# Algorithm B
algoName=${userAlgoName}_B
algoB=$algoName

cat <<EOF> ./${algoName}.h

/**
 * Created automatically with MAFalda (${thedate})
 * MAFalda Author:  John Idarraga <idarraga@cern.ch>
 *
 * This is the second algorithm in a set of two communicating 
 * algorithms. ${algoA} puts an object in the StoreGate and
 * ${algoB} retreives it.
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
