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

#ifndef MediPixWriteToEntuple_h
#define MediPixWriteToEntuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TTree.h>
#include <TMath.h>
#include <vector>
#include <map>
#include "WriteToEntuple_defs.h"

/** Implements what is needed to write the 
 *  frames info and the MetaData 
 *  to an Ntuple ROOT file.
 */

class WriteToNtuple {

public:

  WriteToNtuple();
  WriteToNtuple(TString);
  virtual ~WriteToNtuple(){};
  //Int_t fillVars(std::vector<Int_t>, std::vector<Int_t>, std::vector<Int_t>);
  void includeTree(TString);
  void writeTree(TString);
  void closeNtuple();
  //void cleanUpArrays();

  TTree * getTree(TString);
  TFile * getFile(){return outputROOTFile;};

private:
  
  TFile * outputROOTFile;
  std::map<TString ,TTree *> outputTrees;

  ClassDef(WriteToNtuple,1)
};


#endif

