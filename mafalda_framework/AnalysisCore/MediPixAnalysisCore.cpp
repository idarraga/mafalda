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

#ifndef MediPixAnalysisCore_cxx
#define MediPixAnalysisCore_cxx

#ifndef MediPixAnalysisCore_h
#include "MediPixAnalysisCore.h"
#endif 

#include "TRint.h"

#include "WriteToEntuple/MediPixWriteToEntuple.h"
#include "MPXAlgo/MediPixAlgo.h"

MediPixAnalysisCore::MediPixAnalysisCore(TTree *tree)
{

	// Set the Log service
	Log.setAlgoName("MediPixAnalysisCore");
	Log.OutputLevel = MSG::INFO;

	//g_theApp = gROOT->GetApplication();
	//g_theApp->SetReturnFromRun(true);
	//if(!g_theApp)

	g_theApp = new TApplication("MAFalda", 0, NULL);
	//TRint l;

	//g_theApp = gApplication;
	//std::cout << g_theApp << std::endl;

	Init(tree);

	/* timers */
	algoTimers = new MediPixAlgoTimer;
	/* signals */
	algoSignals.nextFrame = DONT_JUMP;
	algoSignals.direction = SEEK_FORWARD;

	outputNtuple = 0;
	m_outputFilename = "MPXResults.root"; // output filename by default

	// fast matrix
	m_frameXC_Fast = 0x0;
	m_frameXToA_Fast = 0x0;
	m_frameXFastToA_Fast = 0x0;

	// TOT cut
	m_totMaxCut = 1000000; // timepix cut, settable through toplayer
}

MediPixAnalysisCore::~MediPixAnalysisCore()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

void MediPixAnalysisCore::Loop(Long64_t initFrame_i, Long64_t lastFrame_i)
{

	m_initFrame = initFrame_i;
	m_lastFrame = lastFrame_i;

	MediPixAnalysisCore::Loop();
}

void MediPixAnalysisCore::Loop(Long64_t oneFrame_i)
{

	m_initFrame = oneFrame_i;
	m_lastFrame = oneFrame_i;

	MediPixAnalysisCore::Loop();
}

void MediPixAnalysisCore::Loop()
{

	//   In a ROOT session, you can do:
	//      Root > .L MediPixAnalysisCore.C
	//      Root > MediPixAnalysisCore t
	//      Root > t.GetEntry(12); // Fill t data members with entry number 12
	//      Root > t.Show();       // Show values of entry 12
	//      Root > t.Show(16);     // Read and show values of entry 16
	//      Root > t.Loop();       // Loop on all entries
	//

	//    This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch

	/* calling algos to pre initialize */
	std::map<Int_t, string>::iterator algScheduleItr;
	string tempName = "";
	for ( algScheduleItr=m_algosSchedule.begin() ; algScheduleItr != m_algosSchedule.end(); algScheduleItr++ )
	{
		tempName = (*algScheduleItr).second;
		// timers to ntuple information initialization
		m_timePerFrameMap[tempName] = 0.;
		// actual pre-init
		m_algosMap[tempName]->InitMPXAlgo(tempName ,this);
	}

	/* Include algorithms in the Output Ntuple */
	SetOutputNtuple();

	/* Initialization */
	if (fChain == 0) return;
	Long64_t nbytes = 0, nb = 0;

	for ( algScheduleItr=m_algosSchedule.begin() ; algScheduleItr != m_algosSchedule.end(); algScheduleItr++ )
	{
		tempName = (*algScheduleItr).second;
		m_algosMap[tempName]->PreInitMPXAlgo();
		m_algosMap[tempName]->Init();
		m_algosMap[tempName]->PostInitMPXAlgo();
	}


	Log << MSG::INFO << ">------------- | run starts | frame: " << m_initFrame << " -------------<" << endreq;

	for (Long64_t jentry = m_initFrame ; jentry <= m_lastFrame ; )
	{

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		// At this point the Baskets have been opened.
		// Reading the pixel matrix (std::map) might be slow at this point.
		// Converting to a C-type array in order to fasten the reading
		RemakePixelMap();

		// Calling algos for execution
		for (algScheduleItr = m_algosSchedule.begin() ; algScheduleItr != m_algosSchedule.end() ; algScheduleItr++ ) {

			tempName = (*algScheduleItr).second;

			/* Set counters first and the call Execute() */
			m_algosMap[tempName]->SetGlobalCounters(jentry);

			algoTimers->ContinueATimer(tempName);
			m_algosMap[tempName]->Execute();
			algoTimers->StopATimer(tempName);
			// timer info for ntuple
			m_timePerFrameMap[tempName] = algoTimers->GetLastSlotRealTime(tempName); // [s]

			/* Implementing the direction.
			 */
			if(m_algosMap[tempName]->SignalFlag())
			{
				algoSignals = m_algosMap[tempName]->SignalToSteering();
				if(algoSignals.nextFrame != DONT_JUMP) /* if we want to jump to other frame */
				{
					jentry = algoSignals.nextFrame;// - 1;  /* Because it will increment for next loop */
					//algoSignals.nextFrame = DONT_JUMP; /* Jump just once */
				}
			}

		}

		// Finised to run all Algos (Execute).  Tell the manager to do admin stuff
		//  like taking care of the StoreGate
		//Int_t nObjs = 0;
		theManager->GetStoreGate()->CleanUpAllStoreGateExcept(MPXDefs::CONF);

		// if (Cut(ientry) < 0) continue;
		if (algoSignals.direction == SEEK_FORWARD && algoSignals.nextFrame == DONT_JUMP) {
			jentry++;
		} else if (algoSignals.direction == SEEK_BACKWARDS && algoSignals.nextFrame == DONT_JUMP) {
			jentry--;
		} else if (algoSignals.direction == DONT_MOVE)
		{ ; }
		//else{
		// jentry++;}

		// Check here if any algorithm has been scheduled for disconnection
		if ( !m_scheduledForDisconnection.empty() ) {

			vector<string>::iterator itrD = m_scheduledForDisconnection.begin();
			for( ; itrD != m_scheduledForDisconnection.end() ; itrD++ ) {
				Log << MSG::INFO << "Disconnecting algorithm " << *itrD << endreq;
				DisconnectAlgo(*itrD);
			}
		}
		// and then clear the scheduleForDisconnection vector
		m_scheduledForDisconnection.clear();

		// Clean up stuff cleanable on an event per event basis
		theManager->CleanUpCalibMap();

	}

	Log << MSG::INFO << ">------------- | run ends | frame: " << m_lastFrame << " -------------<" << endreq;

	for ( algScheduleItr=m_algosSchedule.begin() ; algScheduleItr != m_algosSchedule.end(); algScheduleItr++ )
	{
		tempName = (*algScheduleItr).second;
		m_algosMap[tempName]->Finalize();
	}

	// Finised to run all Algos (Finalize).  Tell the manager to do admin stuff
	//  like taking care of the StoreGate
	//Int_t nObjs = 0;
	theManager->GetStoreGate()->CleanUpAllStoreGateExcept(MPXDefs::CONF);//CleanUpAllStoreGate(); // clean everything !
	//std::cout << "[AnalysisCore] " << nObjs << " obj(s) erased from StoreGate" << std::endl;

	if(outputNtuple)
		closeNtuple();

	if(algoTimers)
		algoTimers->DumpTimers();

	Log << MSG::INFO << "done." << endreq;
	//g_theApp->Run();
}

Bool_t MediPixAnalysisCore::ConnectAlgo(string algoName, MediPixAlgo * oneAlgo){

	m_algosMap[algoName] = oneAlgo;
	m_algosSchedule[algosCntr] = algoName;
	//oneAlgo->SetAlgoName(algoName); // done in pre initialization
	algosCntr++;

	/* create timers for each algo */
	if(algoTimers){
		algoTimers->SetupATimer(algoName);
	}

	return true;
}

Bool_t MediPixAnalysisCore::DisconnectAlgo(string algoName){

	// erase from algosMap
	map<string, MediPixAlgo *>::iterator itrM;
	itrM = m_algosMap.find(algoName); // already checked if it is there (Manager)
	m_algosMap.erase(itrM);

	// erase from algosSchefule --> need to find index
	map<Int_t, string>::iterator itrS = m_algosSchedule.begin();
	Int_t scheduleIndx = -1;
	for ( ; itrS != m_algosSchedule.end() ; itrS++) {
		if( (*itrS).second == algoName ) {
			scheduleIndx = (*itrS).first;
		}
	}
	itrS = m_algosSchedule.find(scheduleIndx);
	m_algosSchedule.erase(itrS);

	algosCntr--;

	if(algoTimers){
		algoTimers->StopATimer(algoName);
	}

	return true;
}

/**
 * The managers informs that the algorithm "algoName" should be
 * disconnected starting from the next event.  The disconnection
 * must not happen at this point.  The analysis core will erase it
 * from the schedule list at the end of the event.
 */
Bool_t MediPixAnalysisCore::ScheduleAlgoForDisconnection(string algoName){

	// This is the only action for now.  Schedule algoName for disconnection.
	m_scheduledForDisconnection.push_back(algoName);

	return true;
}

void MediPixAnalysisCore::DumpAlgoList(){

	std::map<string,MediPixAlgo *>::iterator it;
	std::cout << "[INFO] Algorithms list:" << std::endl;
	for ( it=m_algosMap.begin() ; it != m_algosMap.end(); it++ )
	{
		std::cout << "       " << (*it).first << " => " << (*it).second << std::endl;
	}
}


void MediPixAnalysisCore::InitMessage(){

	std::cout << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;
	std::cout << "-                                          MAFalda -" << std::endl;
	std::cout << "-           Analysis Framework for pixel detectors -" << std::endl;
	std::cout << "-                 John Idarraga <idarraga@cern.ch> -" << std::endl;
	std::cout << "----------------------------------------------------" << std::endl;
	std::cout << std::endl;

}

void MediPixAnalysisCore::SetOutputNtupleFilename(const Char_t * fileName){

	m_outputFilename = fileName;

}

void MediPixAnalysisCore::SetOutputNtuple(){

	// if no file has bee setup
	outputNtuple = new WriteToNtuple(m_outputFilename.c_str());

	Log << MSG::INFO << "Including a Tree for each algo listed." << endreq;

	/* calling algos to be added as brances in the OutputNtuple */
	std::map<string, MediPixAlgo *>::iterator algIt;
	for ( algIt=m_algosMap.begin() ; algIt != m_algosMap.end(); algIt++ )
	{
		//std::cout << (*algIt).first << " ,";
		outputNtuple->includeTree((*algIt).first);
		// make some default branches
		//Log << MSG::INFO << ((*algIt).second)->GetAlgoName() << endreq;
		((*algIt).second)->getMyTree()->Branch("timePerFrame", &m_timePerFrameMap[(*algIt).first], "timePerFrame/D");
	}
	// include a default tree with MetaData
	//outputNtuple->includeTree("MetaData");
	//outputNtuple->getTree("MetaData")->SetBranch();

	//closeNtuple();

}

void MediPixAnalysisCore::closeNtuple(){

	std::map<string,MediPixAlgo *>::iterator algIt;
	for ( algIt=m_algosMap.begin() ; algIt != m_algosMap.end(); algIt++ )
	{
		outputNtuple->writeTree((*algIt).first);
	}
	outputNtuple->closeNtuple();
	Log << MSG::INFO << "ntuple closed" << endreq;

}

void MediPixAnalysisCore::RemakePixelMap(){

	if ( m_frameXC_Fast == 0x0 ) { // first time

		Log << MSG::ALWAYS << "Creating fast matrix map" << endreq;
		m_frameXC_Fast = new int[fWidth*fHeight];

	} else if ( fWidth != fWidth_previous || fHeight != fHeight_previous ) { // in case matrix size changes during run

		delete m_frameXC_Fast;
		m_frameXC_Fast = new int[fWidth*fHeight];

	}

	// ToA
	if ( m_frameXToA_Fast == 0x0 ) { // first time

		Log << MSG::ALWAYS << "Creating fast matrix ToA map" << endreq;
		m_frameXToA_Fast = new int[fWidth*fHeight];

	} else if ( fWidth != fWidth_previous || fHeight != fHeight_previous ) { // in case matrix size changes during run

		delete m_frameXToA_Fast;
		m_frameXToA_Fast = new int[fWidth*fHeight];

	}

	// Fast ToA
	if ( m_frameXFastToA_Fast == 0x0 ) { // first time

		Log << MSG::ALWAYS << "Creating fast matrix FastToA map" << endreq;
		m_frameXFastToA_Fast = new int[fWidth*fHeight];

	} else if ( fWidth != fWidth_previous || fHeight != fHeight_previous ) { // in case matrix size changes during run

		delete m_frameXFastToA_Fast;
		m_frameXFastToA_Fast = new int[fWidth*fHeight];

	}

	// clean up always
	for( int i = 0 ; i < fWidth*fHeight ; i++ )  {
		m_frameXC_Fast[i] = 0;
		m_frameXToA_Fast[i] = 0;
		m_frameXFastToA_Fast[i] = 0;
	}

	// copy information from the map
	std::map<int,int>::iterator itr = m_frameXC.begin();
	int arraylength = fWidth*fHeight;

	for ( ; itr != m_frameXC.end() ; itr++) {

		if( (*itr).first > arraylength || (*itr).first < 0 ) { // out of pad
			std::cout << "[{AnalysisCore} WARNING] Attempt to fill up entries outside of the pad --> inconsistent data !" << std::endl;
			std::cout << "                         Trying to recover by avoiding this entry.  Position : " << (*itr).first << std::endl;
			std::cout << "                         in a " << fWidth << "x" << fHeight << " matrix." << std::endl;
			continue;
		}

		// Fill counts
		m_frameXC_Fast[ (*itr).first ] = (*itr).second;

	}

	fWidth_previous = fWidth;
	fHeight_previous = fHeight;

	// ToA
	itr = m_frameXToA.begin();
	arraylength = fWidth*fHeight;
	for ( ; itr != m_frameXToA.end() ; itr++) {
		if( (*itr).first > arraylength || (*itr).first < 0 ) { // out of pad
			std::cout << "[{AnalysisCore} WARNING] Attempt to fill up entries outside of the pad --> inconsistent data !" << std::endl;
			std::cout << "                         Trying to recover by avoiding this entry.  Position : " << (*itr).first << std::endl;
			std::cout << "                         in a " << fWidth << "x" << fHeight << " matrix." << std::endl;
			continue;
		}
		// Fill counts
		m_frameXToA_Fast[ (*itr).first ] = (*itr).second;
	}

	// FastToA
	itr = m_frameXFastToA.begin();
	arraylength = fWidth*fHeight;
	for ( ; itr != m_frameXFastToA.end() ; itr++) {
		if( (*itr).first > arraylength || (*itr).first < 0 ) { // out of pad
			std::cout << "[{AnalysisCore} WARNING] Attempt to fill up entries outside of the pad --> inconsistent data !" << std::endl;
			std::cout << "                         Trying to recover by avoiding this entry.  Position : " << (*itr).first << std::endl;
			std::cout << "                         in a " << fWidth << "x" << fHeight << " matrix." << std::endl;
			continue;
		}
		// Fill counts
		m_frameXFastToA_Fast[ (*itr).first ] = (*itr).second;
	}


}

#endif
