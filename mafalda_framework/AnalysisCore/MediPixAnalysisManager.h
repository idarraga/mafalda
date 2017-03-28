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

#ifndef MediPixAnalysisManager_h
#define MediPixAnalysisManager_h

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TDirectory.h>
#include <TVector3.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include "MediPixAnalysisCore.h"
#include "MPXAlgo/MediPixAlgo.h"
#include "MPXAlgo/OutputMng.h"
#include "MPXStoreGate/MPXStoreGate.h"
#include "MPXStoreGate/CandidateContainer.h"

using namespace std;

/**
 * AnalysisManager is the interface between the 
 *   core of the framework and the algorithms.
 * It basically protects the data from the user ;)
 *
 */

class AnalysisManager {

public:

	AnalysisManager(){};
	AnalysisManager(const Char_t *);
	~AnalysisManager();
	OutputMng Log;
	MSG::Endreq endreq;

	Bool_t ConnectAlgo(TString, MediPixAlgo *);
	Bool_t DisconnectAlgo(TString);
	void DumpAlgoList();
	void SetOutputNtupleFilename(const Char_t *);
	void CheckOneFile(const Char_t *, TString);
	map<TString, Int_t> GetInputFilesMap(){return m_filesMap;};
	void SetOutputNtupleFilename();
	void Run();
	void Run(Int_t);
	void Run(Int_t, Int_t);

	Bool_t InputFilesOK(const Char_t *);

	// StoreGate related members
	MediPixStoreGate * GetStoreGate(){return storeGate;};
	void SetStoreGateOutputLevel(MSG::Level l){storeGate->SetOutputLevel(l);};

	/*
	 //Retreive data safely
	inline Int_t GetMatrixElement(Int_t col, Int_t row)
	{
		int xplace = row*analysisCore->fWidth + col;
		if(analysisCore->m_frameXC.find(xplace)
				== analysisCore->m_frameXC.end()){
			return 0;
		}else{
			return analysisCore->m_frameXC[xplace];
		}

		return 0;
	};
	 */
	inline Int_t GetMatrixElement(Int_t col, Int_t row)
	{
		int xplace = row*analysisCore->fWidth + col;
		return analysisCore->m_frameXC_Fast[xplace];
	};

	inline Int_t GetToA(Int_t col, Int_t row)
	{
		int xplace = row*analysisCore->fWidth + col;
		return analysisCore->m_frameXToA_Fast[xplace];
	};

	inline Int_t GetFastToA(Int_t col, Int_t row)
	{
		int xplace = row*analysisCore->fWidth + col;
		return analysisCore->m_frameXFastToA_Fast[xplace];
	};


	inline Int_t GetLVL1(Int_t col, Int_t row){
		int xplace = row*analysisCore->fWidth + col;
		if(analysisCore->m_lvl1.find(xplace)
				== analysisCore->m_lvl1.end()){
			return -1;
		}else{
			return analysisCore->m_lvl1[xplace];
		}
	};

	inline Double_t GetMatrixElementMCEdep(Int_t col, Int_t row) {
		return analysisCore->m_frameXC_TruthE[row*analysisCore->fWidth + col];
	};

	void SetCalibEnergy(Int_t xy, Double_t e) {
		m_frameCalibEnergyMap[xy] = (Int_t)e;
	};
	Int_t GetCalibEnergy(Int_t xy) {
		return m_frameCalibEnergyMap[xy];
	};

	void CleanUpCalibMap(){
		m_frameCalibEnergyMap.clear();
	}

	std::map<int,int> GetTOTMap()
	{return analysisCore->m_frameXC;};

	inline Int_t GetnEntriesPad()
	{return analysisCore->m_nEntriesPad;};

	inline Int_t GetnHitsInPad()
	{return analysisCore->m_nHitsInPad;};

	inline Int_t GetnChargeInPad()
	{return analysisCore->m_nChargeInPad;};

	inline Bool_t GetisMCData()
	{return analysisCore->m_isMCData;};

	inline Int_t GetfFormat()
	{return analysisCore->fFormat;};

	inline Int_t GetfWidth()
	{return analysisCore->fWidth;};

	inline Int_t GetfHeight()
	{return analysisCore->fHeight;};

	inline Int_t GetfAcq_mode()
	{return analysisCore->fAcq_mode;};

	inline Double_t GetfAcq_time()
	{return analysisCore->fAcq_time;};

	inline TString GetfApplied_filters()
	{return analysisCore->fApplied_filters;};

	inline Double_t GetfAuto_erase_interval()
	{return analysisCore->fAuto_erase_interval;};

	inline Int_t GetfAutoerase_interval_counter()
	{return analysisCore->fAutoerase_interval_counter;};

	inline Bool_t GetfBS_active()
	{return analysisCore->fBS_active;};

	inline TString GetfChipboardID()
	{return analysisCore->fChipboardID;};

	inline Double_t GetfCoinc_live_time()
	{return analysisCore->fCoinc_live_time;};

	//UChar_t --> Byte_t
	inline UChar_t GetfCoincidence_delay()
	{return analysisCore->fCoincidence_delay;};

	//UChar_t --> Byte_t
	inline UChar_t GetfCoincidence_mode()
	{return analysisCore->fCoincidence_mode;};

	inline std::vector<Int_t> GetfCounters()
																		  {return analysisCore->fCounters;};

	inline std::vector<Int_t> GetfDACs()
																		  {return analysisCore->fDACs;};

	inline Int_t GetTHL(){
		std::vector<Int_t> vDAC = analysisCore->fDACs;
		if((Int_t)vDAC.size() < 7) return -1;
		return vDAC[6];
	};

	inline Double_t GetfHV()
	{return analysisCore->fHV;};

	inline Int_t GetfHw_timer()
	{return analysisCore->fHw_timer;};

	inline TString GetfInterface()
	{return analysisCore->fInterface;};

	inline Double_t GetfMpx_clock()
	{return analysisCore->fMpx_clock;};

	inline Int_t GetfMpx_type()
	{return analysisCore->fMpx_type;};

	inline Int_t GetfPolarity()
	{return analysisCore->fPolarity;};

	inline Double_t GetfStart_time()
	{return analysisCore->fStart_time;};

	inline TString GetfStart_timeS()
	{return analysisCore->fStart_timeS;};

	inline Double_t GetfTimepix_clock()
	{return analysisCore->fTimepix_clock;};

	inline Double_t GetfTrigger_time()
	{return analysisCore->fTrigger_time;};

	//inline std::vector<TVector3> GetPrimaryMCVertex()
	//{return analysisCore->m_primaryVertex;};

	inline Int_t GetPrimaryMCVertex_N()
	{return analysisCore->m_primaryVertex_x.size();};
	inline Double_t GetPrimaryMCVertex_X(int i)
	{return analysisCore->m_primaryVertex_x[i];};
	inline Double_t GetPrimaryMCVertex_Y(int i)
	{return analysisCore->m_primaryVertex_y[i];};
	inline Double_t GetPrimaryMCVertex_Z(int i)
	{return analysisCore->m_primaryVertex_z[i];};

	inline Long_t GetfFrameId()
	{return analysisCore->fFrameId;};

	inline TString GetfMPXDataSetNumber()
	{return analysisCore->fMPXDataSetNumber;};

	inline Int_t GetnEntriesChain()
	{return analysisCore->fChain->GetEntries();};

	inline TApplication * GetApplication()
    {return analysisCore->g_theApp;};

	string GetOutputBaseName(){ return m_outputBaseName; };

	//void SetCalibEnergy(Int_t col, Int_t row, Double_t e);
	//void SetCalibEnergy(Int_t xy, Double_t e);
	//Int_t GetCalibEnergy(Int_t col, Int_t row);
	//Int_t GetCalibEnergy(Int_t xy);

private:

	// Data after Calibration.  This is held here.  Does not belong to the Core as it
	// doesn't come from data
	std::map<int,int> m_frameCalibEnergyMap;

	MediPixAnalysisCore * analysisCore;
	map<TString, MediPixAlgo *> algosMap;
	// maps files with number of entries in each file
	map<TString, Int_t> m_filesMap;
	MediPixStoreGate * storeGate;
	TTree * m_mpxTree;
	TFile * m_mpxFile;
	string m_outputFilename;
	string m_outputBaseName;
	//string m_outputFilenameStripped; // striped filename, without '.root'

};


#endif
