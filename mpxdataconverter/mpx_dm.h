/**
 *  Author: John Idarraga <idarraga@cern.ch>
 *  Universite de Montreal
 *
 *  Classes containing the frames
 */

#ifndef mpx_dm_h
#define mpx_dm 1

#include <iostream>
#include <map>
#include <vector>
#include <stdlib.h>
#include "TROOT.h"
#include "TObject.h"
#include "TH2.h"
#include "mpx_dm_consts.h"


/** Contains only the pixels matrix 
 *  and some members that may go to the
 *  ROOT file.
 */

class FrameContainer {
  
 private:
  
  Int_t m_frameMatrix[MAX_FRAME_ROW][MAX_FRAME_COL];
  Int_t m_nEntriesPad ;
  Int_t m_nHitsInPad;
  Int_t m_nChargeInPad;
  Bool_t m_isMCData; // in case is data coming from mpxGeant4

 public:
  FrameContainer();
  ~FrameContainer(){};
  void FillOneElement(Int_t, Int_t, Int_t);
  void FillElements(Int_t **);
  
  void ResetCountersPad();
  void CleanUpMatrix();
  void SetFrameAsMCData(){m_isMCData = true;};
  Int_t GetEntriesPad(){return m_nEntriesPad;};
  Int_t GetHitsInPad(){return m_nHitsInPad;};
  Int_t GetChargeInPad(){return m_nChargeInPad;};
  Int_t * GetPixelsMatrix();

  ClassDef(FrameContainer,2)
};


/** 
 * Contains the meta data and an object for the frame matrix.  It
 *  inherits the methods to fill up the frameMatrix from class
 *  FrameContainer
 */

class FrameStruct : public FrameContainer {
  
 private:
  
  /* head info */
  Int_t    fFormat;
  Int_t    fWidth;
  Int_t    fHeight;

  /* all metadata */
  Int_t    fAcq_mode;
  Double_t fAcq_time;
  TString  fApplied_filters;
  Double_t fAuto_erase_interval;
  Int_t    fAutoerase_interval_counter;
  Bool_t   fBS_active;
  TString  fChipboardID;
  Double_t fCoinc_live_time;
  Byte_t   fCoincidence_delay;
  Byte_t   fCoincidence_mode;
  std::vector<Int_t> fCounters;
  std::vector<Int_t> fDACs;
  Double_t fHV;
  Int_t    fHw_timer;
  TString  fInterface;
  Double_t fMpx_clock;
  Int_t    fMpx_type;
  Int_t    fPolarity;  
  Double_t fStart_time;
  TString  fStart_timeS;
  Byte_t   fTimepix_clock;
  Double_t fTrigger_time;

  /* joint */
  Long_t   fFrameId;
  TString  fMPXDataSetNumber;

 public:
  FrameStruct(TString);
  ~FrameStruct(){};
  
  // filling Tree methods
  void FillMetaData(TString, Int_t);
  void IncreaseId(){ fFrameId++; }
  void SetId(Int_t id){ fFrameId=id; }
  TString GetDataSet(){ return fMPXDataSetNumber; };
  Int_t GetFrameId(){return fFrameId;};
  void RewindMetaDataValues();
  
  ClassDef(FrameStruct,4)
};


/** 
 * Containst a FrameStruct and implements the methods needed to read
 *  the frame files, load up information to the frame.
 */

class FramesHandler {
   
 private:
   
  Int_t m_nFrames;
  FrameStruct * m_aFrame;

  Bool_t getAFrameMatrix_flag;
  Bool_t getAFrameHist_flag;
  
  Int_t m_metaBit;
  Int_t m_metaCode;
  Bool_t m_ParseAndHold;
  
  Int_t nFrames256x256;
  Int_t nFramesXYC;

  // in case I am dealing with more than one detector
  Int_t m_detID;
  
 public:
  
  FramesHandler(TString);
  ~FramesHandler(){};
  Int_t IdentifyTypeOfInput(TString);

  /***************************************************** 
   *  Filling frame methods 
   */
  /* Reads from txt & dsc files.  Completely fills m_aFrame. */ 
  Bool_t readOneFrame(TString, TString);
  /* load a single frame pixel (X,Y,C), fills m_aFrame */ 
  Bool_t LoadFramePixel(Int_t, Int_t, Int_t);
  /* load a single frame MetaData entry (METAString, Int_t metaCode)*/
  Bool_t LoadFrameMetaData(TString, Int_t);
  void IncreaseCurrentFrameId();
  void SetCurrentFrameId(Int_t id){m_aFrame->SetId(id);};
  void SetAsMCData(){m_aFrame->SetFrameAsMCData();};
  /* *************************************************** */

  Int_t GetDetectorId(){return m_detID;};
  void  SetDetectorId(Int_t id){m_detID = id;};
  void RewindAll();
  TH2I * getAFrameHist(TString, TString, TString);
  Int_t ** getAFrameMatrix(TString, TString);
  TH2I * getHistFrame(Int_t, Int_t *);
  FrameStruct * getFrameStructObject(){return m_aFrame;};
  Bool_t frameSupervisor();
  void parseMetaLine(Char_t *, Int_t);
  
};

#endif
