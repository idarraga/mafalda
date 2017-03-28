/*
 * 	Copyright 2007 John Idarraga
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

/*
 * Started from atomatic generation with TTree::MakeClass The
 * automated message follows.
 *
 * This class has been automatically generated on
 * Tue Jul 24 12:15:39 2007 by ROOT version 5.13/02
 * from TTree MPXTree/Medi/TimePix data
 * found on file: MPXNtuple_TOT3.root
 */

#ifndef MediPixAnalysisCore_h
#define MediPixAnalysisCore_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TApplication.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TVector3.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "WriteToEntuple/MediPixWriteToEntuple.h"
#include "AlgoTimer/MediPixAlgoTimer.h"
#include "MPXViewer/Viewer_defs.h"
#include "MPXAlgo/OutputMng.h"

using namespace std;

//#include "AnalysisCore/MediPixAnalysisManager.h"
class AnalysisManager;
class MediPixAlgo;

/** 
 * Signals I send back to the Steering
 */

class AlgoSignalsHandler {
public:
	inline AlgoSignalsHandler();
	Long64_t nextFrame;
	Int_t direction;
};

AlgoSignalsHandler::AlgoSignalsHandler(){
	nextFrame = 0;
	direction = SEEK_FORWARD;
}

const Int_t kMaxm_primaryVertex = 19;

class MediPixAnalysisCore {

public :

	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain

	//Long64_t        jFrameId; // current frame

	// Declaration of leave types
	std::map<int,int> m_frameXC;
	std::map<int, int> m_frameXToA;
	std::map<int, int> m_frameXFastToA;

	// FIXME -->  I am unable to deserialize here <int, double>
	std::map<int,int> m_frameXC_TruthE;
	std::map<int,int> m_frameXC_E;

	// Fast versions
	int * m_frameXC_Fast;
	int * m_frameXToA_Fast;
	int * m_frameXFastToA_Fast;

	std::map<int, int>   m_lvl1;
	//Int_t           m_frameMatrix[256][256];
	//Double_t        m_frameMatrixMCEdep[MAX_FRAME_ROW][MAX_FRAME_COL]; // MC //
	Int_t           m_nEntriesPad;
	Int_t           m_nHitsInPad;
	Int_t           m_nChargeInPad;
	Bool_t          m_isMCData;

	Int_t           fFormat;
	// size of the matrix
	Int_t           fWidth;
	Int_t           fHeight;
	// check previous size
	Int_t           fWidth_previous;
	Int_t           fHeight_previous;

	Int_t           fAcq_mode;
	Double_t        fAcq_time;
	TString         fApplied_filters;
	Double_t        fAuto_erase_interval;
	Int_t           fAutoerase_interval_counter;
	Bool_t          fBS_active;
	TString         fChipboardID;
	Double_t        fCoinc_live_time;
	UChar_t         fCoincidence_delay;
	UChar_t         fCoincidence_mode;
	std::vector<Int_t>   fCounters;
	std::vector<Int_t>   fDACs;
	Double_t        fHV;
	Int_t           fHw_timer;
	TString         fInterface;
	Double_t        fMpx_clock;
	Int_t           fMpx_type;
	Int_t           fPolarity;
	Double_t        fStart_time;
	TString         fStart_timeS;
	Double_t         fTimepix_clock;
	Double_t        fTrigger_time;

	std::vector<Double_t>        m_primaryVertex_x;   //[m_primaryVertex_]
	std::vector<Double_t>        m_primaryVertex_y;   //[m_primaryVertex_]
	std::vector<Double_t>        m_primaryVertex_z;   //[m_primaryVertex_]

	Long_t          fFrameId;
	TString         fMPXDataSetNumber;


	// List of branches
	//TBranch        *b_FramesData_m_frameMatrix;   //!
	//TBranch        *b_FramesData_m_frameMatrixMCEdep;   //!
	TBranch		   *b_FramesData_m_frameXC;
	TBranch		   *b_FramesData_m_frameXToA;
	TBranch		   *b_FramesData_m_frameXFastToA;
	TBranch		   *b_FramesData_m_frameXC_TruthE;
	TBranch		   *b_FramesData_m_frameXC_E;
	TBranch		   *b_FramesData_m_lvl1;
	TBranch        *b_FramesData_m_nEntriesPad;   //!
	TBranch        *b_FramesData_m_nHitsInPad;   //!
	TBranch        *b_FramesData_m_nChargeInPad;   //!
	TBranch        *b_FramesData_m_isMCData;
	TBranch        *b_FramesData_fFormat;   //!
	TBranch        *b_FramesData_fWidth;   //!
	TBranch        *b_FramesData_fHeight;   //!
	TBranch        *b_FramesData_fAcq_mode;   //!
	TBranch        *b_FramesData_fAcq_time;   //!
	TBranch        *b_FramesData_fApplied_filters;   //!
	TBranch        *b_FramesData_fAuto_erase_interval;   //!
	TBranch        *b_FramesData_fAutoerase_interval_counter;   //!
	TBranch        *b_FramesData_fBS_active;   //!
	TBranch        *b_FramesData_fChipboardID;   //!
	TBranch        *b_FramesData_fCoinc_live_time;   //!
	TBranch        *b_FramesData_fCoincidence_delay;   //!
	TBranch        *b_FramesData_fCoincidence_mode;   //!
	TBranch        *b_FramesData_fCounters;   //!
	TBranch        *b_FramesData_fDACs;   //!
	TBranch        *b_FramesData_fHV;   //!
	TBranch        *b_FramesData_fHw_timer;   //!
	TBranch        *b_FramesData_fInterface;   //!
	TBranch        *b_FramesData_fMpx_clock;   //!
	TBranch        *b_FramesData_fMpx_type;   //!
	TBranch        *b_FramesData_fPolarity;   //!
	TBranch        *b_FramesData_fStart_time;   //!
	TBranch        *b_FramesData_fStart_timeS;   //!
	TBranch        *b_FramesData_fTimepix_clock;   //!
	TBranch        *b_FramesData_fTrigger_time;   //!

	TBranch        *b_m_primaryVertex_x;   //!
	TBranch        *b_m_primaryVertex_y;   //!
	TBranch        *b_m_primaryVertex_z;   //!

	TBranch        *b_FramesData_fFrameId;   //!
	TBranch        *b_FramesData_fMPXDataSetNumber;   //!

	Long64_t m_initFrame;
	Long64_t m_lastFrame;

	TApplication * g_theApp;
	WriteToNtuple * outputNtuple;

	// log service
	OutputMng Log;
	MSG::Endreq endreq;

	/* FIXME !
      Check all these inline definitions.
	 */

	MediPixAnalysisCore(TTree *tree=0);
	virtual ~MediPixAnalysisCore();

	//inline virtual Int_t    Cut(Long64_t entry);
	inline virtual Int_t    GetEntry(Long64_t entry);
	inline virtual Long64_t LoadTree(Long64_t entry);
	inline virtual void     Init(TTree *tree);
	virtual void            Loop();
	virtual void            Loop(Long64_t);
	virtual void            Loop(Long64_t, Long64_t);
	inline virtual Bool_t   Notify();
	inline virtual void     Show(Long64_t entry = -1);

	void InitMessage();
	Bool_t ConnectAlgo(string, MediPixAlgo *);
	Bool_t DisconnectAlgo(string);
	Bool_t ScheduleAlgoForDisconnection(string);
	void SetOutputNtupleFilename(const Char_t *);
	void SetOutputNtuple();
	void closeNtuple();
	void DumpAlgoList();

	inline void ListenToTheManager(AnalysisManager * mngr){theManager = mngr;};

	inline void SetTOTMaxCut(int c){m_totMaxCut = c;};

private:

	// Copy the frameXC map to a C-array.
	// This fasten things up
	void RemakePixelMap();


	MediPixAlgoTimer * algoTimers;
	std::map<string, MediPixAlgo *> m_algosMap;
	std::map<Int_t, string> m_algosSchedule;
	vector<string> m_scheduledForDisconnection;
	Int_t algosCntr;

	AlgoSignalsHandler algoSignals;
	// Have a pointer to the manager.
	AnalysisManager * theManager;

	// output filename (ntuple)
	std::string m_outputFilename;

	// default branches
	std::map<string, Double_t> m_timePerFrameMap;

	// TOT cuts
	int m_totMaxCut;

};

Int_t MediPixAnalysisCore::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t MediPixAnalysisCore::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->IsA() != TChain::Class()) return centry;
	TChain *chain = (TChain*)fChain;
	if (chain->GetTreeNumber() != fCurrent) {
		fCurrent = chain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void MediPixAnalysisCore::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normaly not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	//fChain->SetBranchAddress("m_frameMatrix[256][256]", m_frameMatrix, &b_FramesData_m_frameMatrix);
	//fChain->SetBranchAddress("m_frameMatrixMCEdep[256][256]", m_frameMatrixMCEdep, &b_FramesData_m_frameMatrixMCEdep);

	fChain->SetBranchAddress("m_frameXC", &m_frameXC, &b_FramesData_m_frameXC);
	fChain->SetBranchAddress("m_frameXToA", &m_frameXToA, &b_FramesData_m_frameXToA);
	fChain->SetBranchAddress("m_frameXFastToA", &m_frameXFastToA, &b_FramesData_m_frameXFastToA);

	fChain->SetBranchAddress("m_frameXC_TruthE", &m_frameXC_TruthE, &b_FramesData_m_frameXC_TruthE);
	fChain->SetBranchAddress("m_frameXC_E", &m_frameXC_E, &b_FramesData_m_frameXC_E);
	fChain->SetBranchAddress("m_lvl1", &m_lvl1, &b_FramesData_m_lvl1);
	fChain->SetBranchAddress("m_nEntriesPad", &m_nEntriesPad, &b_FramesData_m_nEntriesPad);
	fChain->SetBranchAddress("m_nHitsInPad", &m_nHitsInPad, &b_FramesData_m_nHitsInPad);
	fChain->SetBranchAddress("m_nChargeInPad", &m_nChargeInPad, &b_FramesData_m_nChargeInPad);
	fChain->SetBranchAddress("m_isMCData", &m_isMCData, &b_FramesData_m_isMCData);
	fChain->SetBranchAddress("fFormat", &fFormat, &b_FramesData_fFormat);
	fChain->SetBranchAddress("fWidth", &fWidth, &b_FramesData_fWidth);
	fChain->SetBranchAddress("fHeight", &fHeight, &b_FramesData_fHeight);
	fChain->SetBranchAddress("fAcq_mode", &fAcq_mode, &b_FramesData_fAcq_mode);
	fChain->SetBranchAddress("fAcq_time", &fAcq_time, &b_FramesData_fAcq_time);
	// TString // fChain->SetBranchAddress("fApplied_filters", &fApplied_filters, &b_FramesData_fApplied_filters);
	fChain->SetBranchAddress("fAuto_erase_interval", &fAuto_erase_interval, &b_FramesData_fAuto_erase_interval);
	fChain->SetBranchAddress("fAutoerase_interval_counter", &fAutoerase_interval_counter, &b_FramesData_fAutoerase_interval_counter);
	fChain->SetBranchAddress("fBS_active", &fBS_active, &b_FramesData_fBS_active);
	// TString // fChain->SetBranchAddress("fChipboardID", &fChipboardID, &b_FramesData_fChipboardID);
	fChain->SetBranchAddress("fCoinc_live_time", &fCoinc_live_time, &b_FramesData_fCoinc_live_time);
	fChain->SetBranchAddress("fCoincidence_delay", &fCoincidence_delay, &b_FramesData_fCoincidence_delay);
	fChain->SetBranchAddress("fCoincidence_mode", &fCoincidence_mode, &b_FramesData_fCoincidence_mode);
	fChain->SetBranchAddress("fCounters", &fCounters, &b_FramesData_fCounters);
	fChain->SetBranchAddress("fDACs", &fDACs, &b_FramesData_fDACs);
	fChain->SetBranchAddress("fHV", &fHV, &b_FramesData_fHV);
	fChain->SetBranchAddress("fHw_timer", &fHw_timer, &b_FramesData_fHw_timer);
	// TString // fChain->SetBranchAddress("fInterface", &fInterface, &b_FramesData_fInterface);
	fChain->SetBranchAddress("fMpx_clock", &fMpx_clock, &b_FramesData_fMpx_clock);
	fChain->SetBranchAddress("fMpx_type", &fMpx_type, &b_FramesData_fMpx_type);
	fChain->SetBranchAddress("fPolarity", &fPolarity, &b_FramesData_fPolarity);
	fChain->SetBranchAddress("fStart_time", &fStart_time, &b_FramesData_fStart_time);
	// TString // fChain->SetBranchAddress("fStart_timeS", &fStart_timeS, &b_FramesData_fStart_timeS);
	fChain->SetBranchAddress("fTimepix_clock", &fTimepix_clock, &b_FramesData_fTimepix_clock);
	fChain->SetBranchAddress("fTrigger_time", &fTrigger_time, &b_FramesData_fTrigger_time);
	//fChain->SetBranchAddress("m_primaryVertex", &m_primaryVertex, &b_FramesData_m_primaryVertex);
	fChain->SetBranchAddress("m_primaryVertex_x", &m_primaryVertex_x, &b_m_primaryVertex_x);
	fChain->SetBranchAddress("m_primaryVertex_y", &m_primaryVertex_y, &b_m_primaryVertex_y);
	fChain->SetBranchAddress("m_primaryVertex_z", &m_primaryVertex_z, &b_m_primaryVertex_z);

	fChain->SetBranchAddress("fFrameId", &fFrameId, &b_FramesData_fFrameId);
	// TString //fChain->SetBranchAddress("fMPXDataSetNumber", &fMPXDataSetNumber, &b_FramesData_fMPXDataSetNumber);

	Notify();

	m_initFrame = 0;
	if (fChain == 0) return;
	m_lastFrame = fChain->GetEntriesFast();

	InitMessage();

}

Bool_t MediPixAnalysisCore::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normaly not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void MediPixAnalysisCore::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
//
//Int_t MediPixAnalysisCore::Cut(Long64_t entry)
//{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
//   return 1;
//}


#endif
