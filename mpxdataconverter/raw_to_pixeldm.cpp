//***************************************************************************
//
//                     raw_to_cluster.cpp - description
//                     ----------------------------
//   simple ROOT macro to read raw source scan data of a single FE-I3 and 
//   generate a TTree with the clustered hit information.
//   Additional software can be used to load the calibration data and to 
//   analyse the data. 
//   More Inforamtion in calibTree.cpp and clusterAnalysis.cpp
//
//   author: Chrsitian Gallrapp, Christian.Gallrapp@cern.ch
//           2 September 2010
//
// ***************************************************************************
//	How to use it:
//		.x raw_to_cluster.cpp("yourfile.raw")
//
//	Output:
//		out.root with a tree of clustered hits
//
// ***************************************************************************

#include <vector>
#include <iostream>

#include <math.h>
#include <stdlib.h>

#include <TFile.h>
#include <TTree.h>

#include "MediPixWriteToEntuple.h"
#include "allpix_dm.h"

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float>>+;
#endif

using namespace std;

// TTree and TBranches
vector<float> row, col, tot, lv1;
int trigger, clustersize, BCIDnumOfHits, BCIDnumOfCluster;

TFile *f;
TTree *tree;


void plotBCIDwindow(vector<float> in0, vector<float> in1, vector<float> in2, 
		vector<float> in3, vector<float> in4, vector<float> in5,
		vector<float> in6){
	// plot data from one BCID Window
	// this is a temporary function to check the data from one BCID window

	cout << "vLVL1 lsB   col   row   sca1  sca2  tot" << endl;
	for(int it=0; it < (int)in0.size(); it++){
		cout << in0[it] << "     " << in1[it] << "     " << in2[it] << "     " 
				<< in3[it] << "     " << in4[it] << "     " << in5[it] << "     " << in6[it] << endl;
	}
	cout << endl;
}

void plotHitData(vector<float> hitDataCol, vector<float> hitDataRow, 
		vector<float> hitDataTOT, vector<float> hitDataLVL1,
		vector<float> hitDataClustered){
	// plot hit data from one BCID Window
	// this is a temporary function to check only the hits from one BCID window
	// 

	cout << "col   row   tot   virtualLVL1  clustered" << endl;
	for(int it=0; it < (int)hitDataCol.size(); it++){
		cout << hitDataCol[it] << "     " << hitDataRow[it] << "     " << hitDataTOT[it] << "     " 
				<< hitDataLVL1[it] << "     " << hitDataClustered[it] << "     " << endl;
	}
	cout << endl;
}
bool checkCluster(float col1, float col2, float row1, float row2, float lv11, float lv12){
	// checkCluster takes two hit vectors and checks if they are geometricaly connected 
	// and if there is a temporally connection via the Level1 trigger 
	// data structure if one hit: col = [0], row = [1], tot = [2], lv1 = [3]
	// 

	int deltaCol = 2; // column distance between two pixel to build a cluster
	int deltaRow = 2; // row distance between two pixel to build a cluster
	int deltaLVL1 = 4; // time distance between two pixel to build a cluster

	if(fabs(col1-col2)<deltaCol && fabs(row1-row2)<deltaRow && fabs(lv11-lv12)<deltaLVL1){
		return true;
	}
	else{
		return false;
	}
}

void FillThisEvent(vector<float> in0, vector<float> /*in1*/, vector<float> in2,
		vector<float> in3, vector<float> /*in4*/, vector<float> /*in5*/,
		vector<float> in6, int /*BCIDnumber*/, FramesHandler * frame){

	for(int it=0; it < (int)in0.size(); it++){
			if(in3[it]<160){
				// set pixel info, (x, y, tot)
				frame->LoadFramePixel((int)in2[it], (int)in3[it], (int)in6[it], 0., 0.);
				// set lvl1 (x, y, lvl1)
				frame->SetLVL1((int)in2[it], (int)in3[it], (int)in0[it]);
			}
		}

}

void clustering(vector<float> in0, vector<float> /*in1*/, vector<float> in2,
		vector<float> in3, vector<float> /*in4*/, vector<float> /*in5*/,
		vector<float> in6, int BCIDnumber, TTree *tree){
	// clustering takes the BCIDnumber which is the trigger number of one BCIDwindow and the TTree 
	// where the resulting data should be safed. The BCIDwindo data is stored in the global 
	// variable BCIDwindow, which gets reset for every new BCID window. In the end it puts all the
	// cluster from the BCIDwindow in the TTree.
	// TTree name: cluster
	// Branches of tree: clustersize, trigger, <col>, <row>, <tot>, <lv1>, BCIDnumOfHits, BCIDnumOfCluster
	// 

	vector<float> hitDataCol, hitDataRow, hitDataTOT, hitDataLVL1, hitDataClustered;
	for(int it=0; it < (int)in0.size(); it++){
		if(in3[it]<160){
			hitDataCol.push_back(in2[it]); // column
			hitDataRow.push_back(in3[it]); // row
			hitDataTOT.push_back(in6[it]); // tot
			hitDataLVL1.push_back(in0[it]); // virtual LVL1
			hitDataClustered.push_back(0); // virtual LVL1
		}
	}

	float clusterNumber = 0;
	int unsorted = 1; // sorting variable
	//sort hits to clusters
	for(int i=0;i < (int)hitDataCol.size();i++){
		if(hitDataClustered[i]==0){ // unsorted hit found
			clusterNumber++;
			hitDataClustered[i] = clusterNumber; // add hit to new cluster 
			for(int j=0;j< (int)hitDataCol.size();j++){ // find connected hits which are not yet connected
				if(hitDataClustered[j]==0 && checkCluster(hitDataCol[i], hitDataCol[j], 
						hitDataRow[i], hitDataRow[j], hitDataLVL1[i], hitDataLVL1[j])){
					hitDataClustered[j] = clusterNumber;
					// swap possitions and sort connected hits to cluster
					swap(hitDataCol[unsorted],hitDataCol[j]);
					swap(hitDataRow[unsorted],hitDataRow[j]);
					swap(hitDataTOT[unsorted],hitDataTOT[j]);
					swap(hitDataLVL1[unsorted],hitDataLVL1[j]);
					swap(hitDataClustered[unsorted],hitDataClustered[j]);
					unsorted++;
				}
			}
		}
		else if(hitDataClustered[i]!=0){ //sorted hit found
			for(int j=0;j < (int)hitDataCol.size();j++){ // find connected hits which are not yet connected
				if(hitDataClustered[j]==0 && checkCluster(hitDataCol[i], hitDataCol[j], 
						hitDataRow[i], hitDataRow[j], hitDataLVL1[i], hitDataLVL1[j])){
					hitDataClustered[j] = clusterNumber;
					// swap possitions and sort connected hits to cluster
					swap(hitDataCol[unsorted],hitDataCol[j]);
					swap(hitDataRow[unsorted],hitDataRow[j]);
					swap(hitDataTOT[unsorted],hitDataTOT[j]);
					swap(hitDataLVL1[unsorted],hitDataLVL1[j]);
					swap(hitDataClustered[unsorted],hitDataClustered[j]);
					unsorted++;
				}
			}
		}
	}	// end for

	col.clear();
	row.clear();
	tot.clear();
	lv1.clear();

	int clusterNum=0;
	trigger = BCIDnumber;

	BCIDnumOfHits = hitDataCol.size(); // number of hits in one BCID window
	BCIDnumOfCluster = (int)clusterNumber; // number of cluster in one BCID window

	for(int i=0;i < (int)hitDataCol.size();i++){	// write clusters to Tree

		if(clusterNum==0){
			clusterNum = (int)hitDataClustered[i];
			col.push_back(hitDataCol[i]);
			row.push_back(hitDataRow[i]);
			tot.push_back(hitDataTOT[i]);
			lv1.push_back(hitDataLVL1[i]);
		}
		else if(clusterNum!=hitDataClustered[i]){
			clustersize = col.size();
			tree->Fill();	// add one cluster to the tree
			col.clear();
			row.clear();
			tot.clear();
			lv1.clear();
			clusterNum = (int)hitDataClustered[i];
			col.push_back(hitDataCol[i]);
			row.push_back(hitDataRow[i]);
			tot.push_back(hitDataTOT[i]);
			lv1.push_back(hitDataLVL1[i]);
		}
		else{
			col.push_back(hitDataCol[i]);
			row.push_back(hitDataRow[i]);
			tot.push_back(hitDataTOT[i]);
			lv1.push_back(hitDataLVL1[i]);
		}
	}
	clustersize = col.size();
	tree->Fill(); // add the last cluster to the tree
}



void raw_to_cluster(const char *fname){
	// raw_to_cluster uses a raw data file to read out the hit information and generates a TTree with clustered hits.
	// It takes an entire BCID window to do the clustering for each trigger.

	// root-file, with TTree and Branches to store the clustered hits
	f = TFile::Open("out.root", "RECREATE");
	tree = new TTree("cluster","cluster data::row:col:tot:lv1");
	tree->Branch("row",&row);
	tree->Branch("col",&col);
	tree->Branch("tot",&tot);
	tree->Branch("lv1",&lv1);
	tree->Branch("trigger",&trigger);
	tree->Branch("clustersize",&clustersize);
	tree->Branch("BCIDnumOfHits",&BCIDnumOfHits);
	tree->Branch("BCIDnumOfCluster",&BCIDnumOfCluster);

	// MAFaldaDM
	TString dataset = "USBPixRaw";
	TString tempScratchDir = "";
	WriteToNtuple MPXnTuple(dataset, tempScratchDir);
	FramesHandler * frame = new FramesHandler(dataset);
	frame->SetnX(18);
	frame->SetnY(164);
	frame->SetCurrentFrameId(0);

	// input raw-file
	FILE *in = fopen(fname,"r");
	if(in==0){
		printf("Can't open data file %s\n", fname);
		return;
	}
	char line[2000];
	char triggerErrorCode[3];
	int input1, input2, input3, input4, input5, input6;
	vector<float> in0, in1, in2, in3, in4, in5, in6;

	int BCIDnumber = 0;
	int flagNewTrigger = 0;
	int errorFlag = 0;
	int triggerErrorFlag = 0;
	int hitNum = 0;
	int errorNumber = 0;
	float lvl1 = 0;

	int cntr = 0;
	int linecntr = 0;

	// read till end of file
	while(fgets(line,2000,in)!=0){
		linecntr++;
		// data line
		if(sscanf(line,"%d %d %d %d %d %d",&input1, &input2, &input3, &input4, &input5, &input6)==6){
			flagNewTrigger = 0;
			in0.push_back(lvl1);
			in1.push_back(input1);
			in2.push_back(input2); // col
			in3.push_back(input3); // row
			in4.push_back(input4);
			in5.push_back(input5);
			in6.push_back(input6); // tot
			if(input3>159){
				lvl1++;
			}
			else{
				hitNum++;
			}
		}
		// text line
		else if(strncmp("Raw data: 0x",line,12)==0){
			flagNewTrigger++;
			if(flagNewTrigger==1){
				sscanf(line,"%*s %*s %*8s%2s", &triggerErrorCode);
			}
			else if(flagNewTrigger==2){
				if(BCIDnumber>0 && errorFlag==0){ // new BCID window
					if(triggerErrorFlag==1){
						triggerErrorFlag = 0;
						errorNumber++;
					}
					else{
						clustering(in0, in1, in2, in3, in4, in5, in6, BCIDnumber, tree);
						FillThisEvent(in0, in1, in2, in3, in4, in5, in6, BCIDnumber, frame);
						MPXnTuple.fillVars(frame);
						frame->SetnX(18);
						frame->SetnY(160);
						frame->IncreaseCurrentFrameId();
						if(++cntr%1000 == 0){cout << "frame : " << cntr << endl;}

					}
				}
				else if(BCIDnumber>0 && errorFlag==1){ // cluster with at least one error
					errorFlag = 0;
					errorNumber++;
				}
				if(strncmp(triggerErrorCode,"01",2)!=0 && strncmp(triggerErrorCode,"02",2)!=0 && strncmp(triggerErrorCode,"03",2)!=0){
					triggerErrorFlag = 1;
				}
				in0.clear();
				in1.clear();
				in2.clear();
				in3.clear();
				in4.clear();
				in5.clear();
				in6.clear();
				lvl1=0;
				BCIDnumber++;
			}
		}
		// error line
		else if(strncmp("ERROR",line,5)==0){
			errorFlag = 1;
		}
		// last data line
		else if(strncmp("End of dfifo block read",line,23)==0){
			if(errorFlag==0 && triggerErrorFlag==0){
				clustering(in0, in1, in2, in3, in4, in5, in6, BCIDnumber, tree);
			}
			else{
				errorNumber++;
			}
			break;
		}
		// comments
		else if(strncmp("#",line,1)==0){}
		else {
			cout << "Congratulation you found a case, which is not covered by this code!" << endl;
			cout << "Please contact christian.gallrapp@cern.ch in order to fix this bug!" << endl;
		}
	} // end while
	cout<< "Number of Triggers = " << BCIDnumber <<endl;
	cout<< "Number of Hits     = " << hitNum <<endl;
	cout<< "Number of Errors   = " << errorNumber <<endl;

	MPXnTuple.closeNtuple();

	f->Write();
	f->Close();
}

int main(int argc, char ** argv){

	if(argc < 2){
		cout << "use: " << endl;
		cout << "   " << argv[0] << " raw_file(string)" << endl;
		exit(1);
	}

	raw_to_cluster(argv[1]);

	return 0;
}
