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

#ifndef ViewerControl_h
#define ViewerControl_h

#include <iostream>
#include <fstream>
#include <string>

#include <TROOT.h>
//#include <TApplication.h>
#include <TVirtualX.h>
#include <TVirtualPadEditor.h>
#include <TGResourcePool.h>
#include <TGListBox.h>
#include <TGListTree.h>
#include <TGFSContainer.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGMsgBox.h>
#include <TGMenu.h>
#include <TGCanvas.h>
#include <TGComboBox.h>
#include <TGTab.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>
#include <TGFileDialog.h>
#include <TGTextEdit.h>
#include <TGShutter.h>
#include <TGProgressBar.h>
#include <TGColorSelect.h>
#include <RQ_OBJECT.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TEnv.h>
#include <TFile.h>
#include <TKey.h>
#include <TGDockableFrame.h>
#include <TGFontDialog.h>
#include <TG3DLine.h>
#include <TGSimpleTable.h>
#include <TGListBox.h>

#include <TControlBar.h>
#include <TGTextEntry.h>
#include <TGFrame.h>
#include <TROOT.h>

#include <iostream>
#include "TGWindow.h"
#include <RQ_OBJECT.h>
#include <TString.h>

#include "MPXAlgo/MediPixAlgo.h"
#include "MPXAlgo/OutputMng.h"

using namespace std;

// enum for DAQ falues
typedef enum {
  _frameId = 0,
  _frameIdLocal,
  _fileName,
  _acqTime,
  _nHits,
  _nCounts,
  _Thl,
  _StS
} BoxesIndex;

typedef struct {
  Long_t frameId; 
  Long_t frameIdLocal;
  TString currentFile;
  Int_t currentFileId;
  Long_t nEntries;
  Long_t nEntriesLocal;
  Float_t acqTime;
  Int_t nHitsInPad;
  Long_t nChargeInPad;
  Int_t vTHL;
  TString startTimeS;
  Double_t startTimeD;
} ControlUpdateInfo;

typedef struct {
  Int_t direction;
  Int_t oldDirection;
} ViewerSteer ;

typedef enum {
	__TOT = 0,
	__Calib,
	__LVL1,
	__MC,
	__TRUTH_MC
} plotVar;

class ControlFeedbackInfo {

public:

  Int_t maxHisto; 
  Int_t minHisto; 
  Long64_t jumpTo; 
  Bool_t autoAdjust; 
  Bool_t m_plotTOT;
  Bool_t m_plotCalib;
  Bool_t m_plotLvl1;
  Bool_t m_plotMC;
  Bool_t m_plotTruthMC;
  Bool_t m_disconnect;

  Bool_t operator==(const ControlFeedbackInfo & ) const;
  Bool_t operator!=(const ControlFeedbackInfo & ) const;

  Bool_t operator!() const;

  void SetPlotExclusive(plotVar i){
	  if(i == __TOT)      { m_plotTOT = true;  m_plotCalib = false; m_plotLvl1 = false; m_plotMC = false; m_plotTruthMC = false; }
	  if(i == __Calib)      { m_plotTOT = false;  m_plotCalib = true; m_plotLvl1 = false; m_plotMC = false; m_plotTruthMC = false; }
	  if(i == __LVL1)     { m_plotTOT = false; m_plotCalib = false; m_plotLvl1 = true;  m_plotMC = false; m_plotTruthMC = false; }
	  if(i == __MC)       { m_plotTOT = false; m_plotCalib = false; m_plotLvl1 = false; m_plotMC = true;  m_plotTruthMC = false; }
	  if(i == __TRUTH_MC) { m_plotTOT = false; m_plotCalib = false; m_plotLvl1 = false; m_plotMC = false; m_plotTruthMC = true; }
  };

  inline void operator=(ControlFeedbackInfo * oldFeedback){
    maxHisto = oldFeedback->maxHisto;
    minHisto = oldFeedback->minHisto;
    jumpTo = oldFeedback->jumpTo;
  }

};

// class TextMargin defined at the end of this file
class TextMargin;

class MPXViewerControl : public TGMainFrame {

RQ_OBJECT("MPXViewerControl")

  protected:

  // Log
  OutputMng Log;
  MSG::Endreq endreq;

  // ListBox of DAQ values
  TGListBox * fListBox; 

  TextMargin ** fconfBoxes;

  ////////////////////////////////////////////


  TGCompositeFrame  * fMatrixControl;
  TGCompositeFrame  * fSteeringControl;

  // Composite Frame inside the fMain frame to put algorithm info
  TGCompositeFrame  * fAlgoControl;

  // auto-adjust
  TGCheckButton * fdisableAutoAdjust;

  // mapping between keys "Algo:ConfigurationName" and TextMargin boxes pointers 
  std::map<TString, TextMargin *> m_confMapping;
  // maps of Algorithms and number of parameters to show
  std::map<TString, Int_t> m_algorithmMap;

   // Containers for algorithm Configuration
   TGNumberEntry ** fAlgoNumberEntry_double;
   TGNumberEntry ** fAlgoNumberEntry_float;
   TGNumberEntry ** fAlgoNumberEntry_int;

   TGLabel ** fAlgoLabels_double;
   TGLabel ** fAlgoLabels_float;
   TGLabel ** fAlgoLabels_int;

   TGTextButton      * fButtonReprocessEvent;

   TGNumberEntry     * fHistoMin;
   TGNumberEntry     * fHistoMax;
   TGTextButton      * fButtonSetMinMax;

   TGCompositeFrame  * fStatusFrameMedaData;
   TGCompositeFrame  * fStatusHistos;

   TGTextButton      * fButtonQuit;

   TGTextButton      * fButtonSetHisto1;
   TGTextButton      * fButtonSeekBack;
   TGTextButton      * fButtonSeekForward;

   TGTextButton      * fSaveConfiguration;

   TGNumberEntry     * fJumpToFrame;
   TGTextButton      * fButtonJumpToFrame;

   TGCanvas          * fCanvas;
   const TGWindow    * m_p;

   TGLabel * fLbl_main;
   TGLabel * fLbl_frameId;
   TGLabel * fLbl_acqTime;
   TGLabel * fLbl_nHits;
   TGLabel * fLbl_nCounts;
   TGLabel * fLbl_THL;

   TGHorizontalFrame *fH1;
   ControlFeedbackInfo * feedBackInfoControl;

   // Plot switches
   TGRadioButton * m_plotSwitchTOT;
   TGRadioButton * m_plotSwitchCalib;
   TGRadioButton * m_plotSwitchLvl1;
   TGRadioButton * m_plotSwitchMC; // energy with digitizer effects
   TGRadioButton * m_plotSwitchTruthMC; // only hits energy

   // special objects from algorithms
   std::vector<CandidateContainer *> m_objs_conf_int;
   std::vector<CandidateContainer *> m_objs_conf_float;
   std::vector<CandidateContainer *> m_objs_conf_double;
   std::vector<CandidateContainer *> m_objs_conf_bool;
   Int_t m_objs_conf_itr;

 private:
   TApplication * g_theApp;
   ViewerSteer * vSteerControl;
   int m_matrixSizeX;
   int m_matrixSizeY;

 public:
  MPXViewerControl(TApplication *, ViewerSteer *, std::vector<CandidateContainer *> *, int, int);
  virtual ~MPXViewerControl();
  
  void Build();
  //TGCompositeFrame * GetAlgoControl(){return fAlgoControl;};
  void GenerateAlgorithmMap();
  
  /* slots for GUI events */
  void PrintEventStats();
  void seekForward();
  void seekBack();
  void ShowResults();
  void jump100F();
  void jump100B();
  void hardQuit();
  void softQuit();
  void ReprocessCurrentFrame();
  void DisconnectAndRun();
  void SetHighlightersEnabled(Bool_t);
  void DoLeftMargin(char *);

  TString BuildKeyNameTextMargin(TString a, TString b){
    TString name = a;
    name += ":";
    name += b;
    return name;
  };

  ControlFeedbackInfo * ControlUpdate(ControlUpdateInfo *, ControlFeedbackInfo *);
  
  void DoSetMin();
  void DoSetMax();
  void ToogleAutoAdjust();
  void TooglePlotValue_TOT();
  void TooglePlotValue_Calib();
  void TooglePlotValue_lvl1();
  void TooglePlotValue_MC();
  void TooglePlotValue_TMC();
  
  void SetJumpTo();
  void DoJumpTo();
  
  void SetConfigAlgoValue();
  void SaveConfiguration();
  
  TString CreateStartTimeString(Double_t);

  ClassDef(MPXViewerControl,0)
};


//////////// auxilary class ///////////////////////////////////////////////////

#include <TGNumberEntry.h>
#include <TGToolTip.h>
#include <TBox.h>

class TextMargin : public TGHorizontalFrame {
  
protected:
  
  TGNumberEntry *fEntry;
  
public:

  TextMargin(const TGWindow *p, 
	     const char *name,
	     TGNumberFormat::EStyle style = TGNumberFormat::kNESReal, 
	     TGNumberFormat::EAttribute attr = TGNumberFormat::kNEAAnyNumber, 
	     TGNumberFormat::ELimit limits = TGNumberFormat::kNELNoLimits, 
	     Double_t min = 0, 
	     Double_t max = 1) : TGHorizontalFrame(p) 
  {
    
    fEntry = new TGNumberEntry(this, 0, 6, -1, style, attr, limits, min, max);
    AddFrame(fEntry, new TGLayoutHints(kLHintsLeft));

    TGLabel * label = new TGLabel(this, name);
    AddFrame(label, new TGLayoutHints(kLHintsLeft, 10));

    // try to have a ToolTip for the TGLabel object. See TGTextButton class.
    //TGToolTip * toolTip = new TGToolTip(this, new TBox(1,1,2,2), name, 350);
    //toolTip->Hide();
    //AddFrame(toolTip, new TGLayoutHints(kLHintsLeft, 10));

  }

  TGTextEntry * GetEntry() const { return fEntry->GetNumberEntry(); }
  TGNumberEntry * GetTGNumberEntry() {return fEntry;};
  
  ClassDef(TextMargin, 0)
};

#endif
