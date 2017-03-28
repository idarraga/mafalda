//////////////////////////////////////////
// writeToEntuple 
// John Idarraga 08/2006
//////////////////////////////////////////

#ifndef MediPixWriteToEntuple_h
#define MediPixWriteToEntuple_h 1

#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TMath.h>
#include "allpix_dm_consts.h"


/** Implements what is needed to write the 
 *  frames info and the MetaData 
 *  to an Ntuple ROOT file.
 */

class FramesHandler;
class FrameStruct;

class WriteToNtuple {

public:

  WriteToNtuple(TString, TString);
  ~WriteToNtuple();
  void fillVars(FramesHandler *, bool rmd = true);
  void closeNtuple();
  TString GetNtupleFileName(){return m_ntupleFileName;};

private:
  
  FrameStruct * m_frame;

  TH2 * h1;
  TFile * nt;
  TTree * t2;
  TString m_MPXDataSetNumber;

  TString m_ntupleFileName;

  //ClassDef(WriteToNtuple,1)
};

#endif
