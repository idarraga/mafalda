/////////////////////////////////////////
// John Idárraga
// For the MediPix collaboration 06/2007

#include <vector>
#include <iostream>

#include "MediPixWriteToEntuple.h"
//#include "mpx_dm.h"
#include "allpix_dm.h"

WriteToNtuple::WriteToNtuple(TString dataSet, TString tempScratchDir){

  m_MPXDataSetNumber = dataSet;
  m_ntupleFileName = "";

  // change dir if requested
  if(tempScratchDir.Length() > 0){
    m_ntupleFileName += tempScratchDir;
    m_ntupleFileName += "/";
  }
  
  m_ntupleFileName += "MPXNtuple_"+m_MPXDataSetNumber+".root";
  nt = new TFile(m_ntupleFileName, "RECREATE");
  t2 = new TTree("MPXTree","Medi/TimePix data");

  m_frame = new FrameStruct(m_MPXDataSetNumber);
  t2->Branch("FramesData", "FrameStruct", &m_frame, 128000, 2);
  
}

WriteToNtuple::~WriteToNtuple(){

	if(t2) t2->Delete(); // First the TTree
	if(nt) nt->Delete();

	if(m_frame) delete m_frame;

}

void WriteToNtuple::fillVars(FramesHandler * frameHandlerObj, bool rewind_metadata){

  // Variables(class) in the Tree
  m_frame = frameHandlerObj->getFrameStructObject();
  nt->cd();
  // fill the Tree
  t2->Fill();
  // clean up
  frameHandlerObj->RewindAll(rewind_metadata);

}

void WriteToNtuple::closeNtuple()
{

  t2->Write();
  nt->Close();

}
