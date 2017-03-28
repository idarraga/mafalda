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

#ifndef MediPixAlgo_h
#define MediPixAlgo_h


#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>

#include <TROOT.h>
#include <TTree.h>
#include <TString.h>
#include <TVector3.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TF2.h>

#include "AnalysisCore/MediPixAnalysisCore.h"
#include "AnalysisCore/MediPixAnalysisManager.h"
#include "MPXStoreGate/CandidateContainer.h"
#include "MPXStoreGate/MPXStoreGate.h"
#include "AnalysisCore/AnalysisCore_defs.h"

#include "ConfigurationValue.h"
#include "OutputMng.h"

using namespace std;

class Highlighter;
class blob;

/** 
 * Algos must inheriting from this class.  I will interface standard
 *  services here.
 */

class MediPixAlgo {

 public:

  MediPixAlgo();
  virtual ~MediPixAlgo() { };

  OutputMng Log;
  void changeOutputLevel(MSG::Level newL){Log.OutputLevel = newL;};
  void OutputLevel(MSG::Level newL){changeOutputLevel(newL);}; // same as changeOutputLevel
  void InitMPXAlgo(TString, MediPixAnalysisCore *);
  void PreInitMPXAlgo();
  void PostInitMPXAlgo();
  TString GetAlgoName(){return algoName;};
  void SetAlgoName(TString an){algoName = an;};
  void ListenToTheManager(AnalysisManager *);
  void SetConfigurationFilename(const Char_t *);

  virtual void Execute() = 0;
  virtual void Init() = 0;
  virtual void Finalize() = 0;
  TString authorSignature;

  // Extended mask
  void UseExtendedMaskInput(bool in){ m_inputExtendedMask = in; };
  void UseExtendedMaskFile(TString extendedMaskFilename){ m_extendedMaskFile = extendedMaskFilename; m_inputExtendedMask = true; };
  void LoadExtendedMask();
  void DumpExtendedMask();
  bool PixelInExtendedMask(pair<int, int> pix, int sizex);
  bool PixelInExtendedMask(int x, int y, int sizex);
  bool PixelInExtendedMaskBelongsToCluster(blob, int, int);
  bool UsingExtendedMask() { return m_inputExtendedMask; };
  vector<int> GetExtendedMask() { return m_extendedMask; };

  /** ***********************************************
   * DM Interface ! 
   */

  void DisconnectMeFromRun();

  Int_t GetMatrixXdim();
  Int_t GetMatrixWidth(){return GetMatrixXdim();}; // same
  Int_t GetWidth(){return GetMatrixXdim();}; // same

  std::map<int,int> GetTOTMap();

  Int_t GetMatrixYdim();
  Int_t GetMatrixHeight(){return GetMatrixYdim();}; // same
  Int_t GetHeight(){return GetMatrixYdim();}; // same

  Int_t GetMatrixElement(Int_t col, Int_t row);
  Int_t GetMatrixElement(std::pair<Int_t, Int_t>);
  Int_t GetToA(Int_t col, Int_t row);
  Int_t GetToA(std::pair<Int_t, Int_t>);
  Int_t GetFastToA(Int_t col, Int_t row);
  Int_t GetFastToA(std::pair<Int_t, Int_t>);

  // Set pixel map after calibration
  Int_t GetCalibEnergy(Int_t col, Int_t row);
  void SetCalibEnergy(std::pair<Int_t, Int_t> pix, Double_t e);
  void SetCalibEnergy(Int_t col, Int_t row, Double_t e);
  //void CalibrateCluster(cluster b);

  Int_t GetLVL1(Int_t col, Int_t row);
  Int_t GetLVL1(std::pair<Int_t, Int_t>);

  Double_t GetMatrixElementMCEdep(Int_t col, Int_t row);
  Double_t GetMatrixElementMCEdep(std::pair<Int_t, Int_t>);

  UInt_t GetMatrixCropAsBitWord(Int_t, Int_t, Int_t, Int_t);
  std::vector<Int_t> GetRowAsVector(Int_t);
  std::vector<Int_t> GetColAsVector(Int_t);

  Int_t GetEntriesPad();
  Int_t GetHitsInPad();
  Int_t GetChargeInPad(); // alias of GetTotalTOT
  Int_t GetTotalTOT();
  Bool_t GetIsMCData();
  Bool_t IsMCData(); // same than above

  void DumpDACs(MSG::Level);

  Int_t GetFormat();
  Int_t GetFrameWidth();
  Int_t GetFrameHeight();
  Int_t GetAcqMode();
  Double_t GetAcqTime();
  TString GetAppliedFilters();
  Double_t GetAutoEraseInterval();
  Int_t GetAutoeraseIntervalCounter();
  Bool_t GetBSActive();
  TString GetChipboardID();
  Double_t GetCoincLiveTime();
  UChar_t GetCoincidenceDelay();
  UChar_t GetCoincidenceMode();
  std::vector<Int_t> GetCounters();
  std::vector<Int_t> GetDAQs();
  Int_t GetTHL();
  Double_t GetHV();
  Double_t GetBiasVoltage();
  Int_t GetHwTimer();
  TString GetInterface();
  Double_t GetMpxClock();
  Int_t GetMpxType();
  Int_t GetPolarity();
  Double_t GetStartTime(); 
  TString GetStartTimeS();
  Double_t GetTimepixClock();
  Double_t GetTriggerTime();
  Long_t GetFrameId();
  TString GetDataSetNumber();
  Int_t GetnEntriesChain();
  Int_t GetNFrames();
  //std::vector<TVector3> GetPrimaryMCVertex();
  Int_t GetPrimaryMCVertex_N();
  Double_t GetPrimaryMCVertex_X(int);
  Double_t GetPrimaryMCVertex_Y(int);
  Double_t GetPrimaryMCVertex_Z(int);

  // error handling --> TODO
  void __FATALMC_ERROR();

  /* ****************************************** */

  /* Other Gets */
  TApplication * GetApplication();
  // Get the TTree corresponding to this algo
  TTree * getMyTree();
  TTree * GetMyTree(){return getMyTree();};
  // Get the ROOTFile;
  TFile * getMyROOTFile();
  TFile * GetMyROOTFile(){return getMyROOTFile();};
  // output
  MSG::Endreq endreq;

  inline void IncreaseGlobalCounters(){m_globalCounter++;};
  inline void DecreaseGlobalCounters(){m_globalCounter--;};
  inline void SetGlobalCounters(Int_t i){m_globalCounter = i;};
  inline Long_t GetGlobalCounters(){return m_globalCounter;};
  inline map<TString, Int_t > GetInputFilesMap(){return m_inputFilesMap;};

  /** Signals to the steer **/
  /*  Just return the same thing
   *  If you can return something different
   *  overload this method
   */
  virtual AlgoSignalsHandler SignalToSteering();
  virtual Bool_t SignalFlag(){return false;};
  AlgoSignalsHandler GetSignals(){return mySignals;};
  void ResetSignals(AlgoSignalsHandler);
  void JumpToFrame(Long64_t, Int_t);

  //////////////////////////
  // Store Gate Connection
  Bool_t PullToStoreGateAccess(CandidateContainer *, MPXDefs::Flags);
  Bool_t PushToStoreGate(CandidateContainer *, MPXDefs::Flags);
  Int_t GetNumberOfSpecialObjects(MPXDefs::SpecialObjs);
  Int_t GetNumberOfSpecialObjectsLessEqualThan(MPXDefs::SpecialObjs);
  Int_t GetNumberOfObjectsWithAuthor(TString);
  CandidateContainer * GetObjectFromAuthor(TString, Int_t);
  vector<CandidateContainer *> GetObjectsSpecial(MPXDefs::SpecialObjs);
  vector<CandidateContainer *> GetObjectsSpecialWithAuthor(MPXDefs::SpecialObjs, TString);
  string GetOutputBaseName();

  ////////////////////////////
  // Configuration
  void TouchConfiguraton(std::string, std::string, std::string, std::string);
  void BadConfigValue(const Char_t *);
  //template <class T> void RegisterConfigurationValue(T *, Char_t *);
  void RegisterConfigurationValue(Double_t *, const Char_t *);
  void RegisterConfigurationValue(Float_t *, const Char_t *);
  void RegisterConfigurationValue(Int_t *, const Char_t *);
  void RegisterConfigurationValue(Bool_t *, const Char_t *);

  Bool_t HasConfiguration(){return m_hasConfiguration;};
  void ReadConfiguration(const Char_t *);
  void ReadConfiguration();

  // Extra drawing
  TCanvas * DrawInSeparateWindow(vector<TGraph2D *>, MSG::Level vl = MSG::INFO, TString = "tri1", TString = "");
  TCanvas * DrawInSeparateWindow(vector<TGraph *>, TString, MSG::Level vl = MSG::INFO);
  TCanvas * DrawInSeparateWindow(vector<vector<TGraph * > >, TString, MSG::Level vl = MSG::INFO);

  TCanvas * DrawInSeparateWindow(vector<TGraph *>, vector<double>, TString, MSG::Level vl = MSG::INFO);
  TCanvas * DrawInSeparateWindow(vector<TGraph *>, vector<TF1 *>, TString, MSG::Level vl = MSG::INFO);
  TCanvas * DrawInSeparateWindow(vector<TH1 *>, vector<TF1 *>, TString, MSG::Level vl = MSG::INFO);
  TCanvas * DrawInSeparateWindow(vector<TH1 *>,  MSG::Level vl = MSG::INFO);
  TCanvas * DrawInSeparateWindow(vector<TH2 *>,  MSG::Level vl = MSG::INFO);
  TCanvas * DrawInSeparateWindow(vector<TF2 *>,  MSG::Level vl = MSG::INFO);

  void FillValuesForDisplay(Highlighter * hl, blob bl);
  void SetMarkerFontSizeDrawingSeparateWindow(double s) { m_MarkerFontSizeDrawingSeparateWindow = s; };

  /////////////////////////////
  // Report objects to the TBrowser
  //void ConnectTBrowser(TBrowser * b1_i)
  //{m_b1 = b1_i;};
  //void AddObjectToBrowser(MAFTArrow * a1_i)
  //{if(m_b1) m_b1->Add((TObject *)a1_i, "MAFTArrow",0);};

 private:

  MediPixAnalysisCore * m_MPXdata;
  TString algoName;
  AnalysisManager * m_myManager;

  /* steer capabilities, available for all algos */
  AlgoSignalsHandler mySignals;

  /* Global frame counter  */
  Long_t m_globalCounter;
  
  /* files map */
  map<TString, Int_t> m_inputFilesMap;

  //////////////////////////
  // Viewer Connection
  ConfigurationValue<Double_t> * confValue_double;
  ConfigurationValue<Float_t> * confValue_float;
  ConfigurationValue<Int_t> * confValue_int;
  ConfigurationValue<Bool_t> * confValue_bool;
  Bool_t m_hasConfiguration;
  Bool_t m_readConfiguration;
  string m_configurationFilename;

  // Other configuration
  double m_MarkerFontSizeDrawingSeparateWindow;

  // Extended mask
  bool m_inputExtendedMask;
  ifstream * infile_extendedmask;
  vector<int> m_extendedMask;
  TString m_extendedMaskFile;

  ClassDef(MediPixAlgo, 0)
};

#endif

