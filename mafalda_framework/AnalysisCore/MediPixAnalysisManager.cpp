/*
 * 	Copyright 2008 John Idarraga
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

#ifndef MediPixAnalysisManager_cpp
#define MediPixAnalysisManager_cpp

#include "MediPixAnalysisManager.h"


AnalysisManager::AnalysisManager(const Char_t * filename){

	// Set the Log service
	Log.setAlgoName("AnalysisManager");

	// Get the TTree
	if(!InputFilesOK(filename)){
		Log << MSG::FATAL << "found problems openning input file(s) ... leaving" << endreq;
		exit(1);
	}

	// Creation of the AnalysisCore Instance
	//   done by the Manager, not by the user.
	analysisCore = new MediPixAnalysisCore(m_mpxTree);
	analysisCore->ListenToTheManager(this);
	Log << MSG::INFO << "Creating MediPixAnalysisCore ("  << analysisCore << ")" << endreq;

	SetOutputNtupleFilename(m_outputFilename.c_str());


	// Creation of a StoreGate intance (by default for now)
	Log << MSG::INFO  << "Creating StoreGate" << endreq;
	storeGate = new MediPixStoreGate();

	// style
	Log << MSG::INFO << "Applying Medipix style settings" << endreq ;
	gROOT->ProcessLine(".x styles/MedipixStyle.C");

}

AnalysisManager::~AnalysisManager(){

	// watch out !, I can't delete the analysisCore ptr from here
	//delete analysisCore;
}

Bool_t AnalysisManager::InputFilesOK(const Char_t * filename){

	m_mpxFile = new TFile(filename,"READ");
	//TTree * T;

	// in case several files are available
	TChain * chain = new TChain("MPXTree");
	string outputfilename = "";
	string basename = "";
	Int_t nFiles = 0;
	size_t npos = 0;
	string onefilename;

	if(m_mpxFile->IsZombie())
	{
		Log << MSG::INFO << "Trying loading a list of ROOT files from \""
				<< filename << "\""<< endreq;
		m_mpxFile->Close();

		ifstream ifs( filename , ifstream::in );
		char temp[512];
		string tempS;

		while (ifs.good())
		{
			ifs.getline(temp,512);
			tempS = temp;

			//size_t npos = 0;
			npos = tempS.find(".root");

			if((int)npos != -1)
			{

				if(tempS[0] == '#'){ // handle comments, not considering the line
					//cout << "Skipping        : [" << tempS << "] " << endl;
					continue;
				}

				//string onefilename;
				npos = tempS.find_last_of("/");
				if((int)npos != -1)
				{
					onefilename = tempS.substr(npos+1, tempS.length()-1);
				}
				else
				{
					onefilename = tempS;
				}

				// Before adding to the chain get the number
				// of entries in this file and information to
				// files map
				CheckOneFile(tempS.c_str(), onefilename);

				// add to the chain
				chain->Add(tempS.c_str());

				if(nFiles == 0) Log << MSG::INFO << "Openning        : [" << nFiles << "] "
						<< onefilename << endreq;
				if(nFiles > 0) Log << MSG::INFO << "Adding to chain : [" << nFiles << "] "
						<< onefilename << " | Cumulative entries = " << chain->GetEntries() << endreq;
				nFiles++;

				// set output filename
				outputfilename = "MAFOutput_";
				outputfilename += onefilename;
				// Store this value ... is going to be used after
				// initialization of MediPixAnalysisCore
				m_outputFilename = outputfilename;

			}
			else
			{
				// A non valid filename in the list
				if((int)tempS.length() > 0)
				{
					Log << MSG::INFO << "\"" << tempS << "\""
							<< " does not seem to be a valid data file ..."
							<< endreq;
					return false;
				}
				// a carriage return has length 0.
			}
		}

		ifs.close();

	}
	else
	{
		string tempS = filename;
		//string onefilename;
		size_t npos = 0;
		npos = tempS.find_last_of("/");

		if((int)npos != -1){
			onefilename = tempS.substr(npos+1, tempS.length()-1);
		}
		else{
			onefilename = tempS;
		}

		// Before adding to the chain get the number
		// of entries in this file and information to
		// files map
		CheckOneFile(tempS.c_str(), onefilename);

		// get the Tree
		chain->Add(filename);
		Log << MSG::INFO << "Input file is a single ROOT file, trying to load \"MPXTree\" from "
				<< onefilename << " | entries = " << chain->GetEntries() << endreq;

		// set output filename
		outputfilename = "MAFOutput_";
		outputfilename += onefilename;
		// store this value ... is going to be used after initialization
		// of MediPixAnalysisCore
		m_outputFilename = outputfilename;

	}

	//else {
	//	return false;
	//}

	// Store the basename, i.e. name without prefix and without ending .root
	// Chop out the .root ending
	npos = onefilename.find_last_of(".root");
	m_outputBaseName = onefilename.substr(0, npos-4);
	// See of the prefix MPXNtuple_ is present.  If so, chop it.
	npos = m_outputBaseName.find_first_of("MPXNtuple_");
	if( npos != string::npos ){
		m_outputBaseName = m_outputBaseName.substr(npos + 10, m_outputBaseName.length()-1);
	}

	//Log << MSG::INFO << "Number of files in the files map : "<< m_filesMap.size() << endreq;

	// deliver a TTree
	m_mpxTree = (TTree *) chain;

	Log << MSG::INFO << "Output filename : " << outputfilename.c_str() << endreq;

	if(!m_mpxTree)
		return false;

	return true;
}

void AnalysisManager::CheckOneFile(const Char_t * fn, TString sfn){

	// fn: has the full path
	// sfn: is the stripped version
	TFile * f = new TFile(fn, "READ");
	TTree * T = (TTree *)f->Get("MPXTree");

	// fill map
	m_filesMap[sfn] = (Int_t)T->GetEntries();

}


void AnalysisManager::DumpAlgoList(){

	std::map<TString,MediPixAlgo *>::iterator it;
	Log << MSG::INFO << "[INFO] Algorithms list:" << endreq;
	for ( it=algosMap.begin() ; it != algosMap.end(); it++ ){
		Log << MSG::INFO << "       " << (*it).first << " => " << (*it).second << endreq;
	}
}

Bool_t AnalysisManager::ConnectAlgo(TString algoName, MediPixAlgo * oneAlgo){

	/* Keep a local copy of the algos list */
	algosMap[algoName] = oneAlgo;
	/* Report them to the analysisCore */
	analysisCore->ConnectAlgo(algoName.Data(), oneAlgo);
	/* Tell the algo who the manager is, but the algo can not directly access
     this information.  Like this a proper interface to the data can
     be defined */
	oneAlgo->ListenToTheManager(this);

	return true;
}

/**
 * The algorithm is informing through MedipixAlgo that it should
 * be disconnected starting from the next event.  The disconnection
 * must not happen at this point.  The analysis core will erase it
 * from the schedule list at the end of the event.
 */
Bool_t AnalysisManager::DisconnectAlgo(TString algoName){

	// Find the algo in the map
	map<TString, MediPixAlgo *>::iterator itr;
	itr = algosMap.find(algoName);

	if(itr == algosMap.end()){
		Log << MSG::INFO << "Trying to disconnect algo " << algoName <<
				" which couldn't be found in the list.  Doing nothing." << endreq;
		return false;
	}

	analysisCore->ScheduleAlgoForDisconnection(algoName.Data());

	return true;
}


void AnalysisManager::Run(){
	// all
	analysisCore->Loop(0, analysisCore->fChain->GetEntries()-1);
}

void AnalysisManager::Run(Int_t initFrame){
	analysisCore->Loop(initFrame);
}

void AnalysisManager::Run(Int_t initFrame, Int_t endFrame){
	analysisCore->Loop(initFrame, endFrame);
}

void AnalysisManager::SetOutputNtupleFilename(const Char_t * fileName){
	analysisCore->SetOutputNtupleFilename(fileName);
}

void AnalysisManager::SetOutputNtupleFilename(){
	analysisCore->SetOutputNtupleFilename("MPXResults.root"); // if no output name is given
}

#endif
