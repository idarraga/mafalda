/////////////////////////////////////////
// John Idárraga
// For the MediPix collaboration 06/2007

#include <TROOT.h>
#include <TString.h>


#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

using namespace std;

/** Handles the list of files in a given directory
 *  that corresponds to frame files to be read.  It calculates
 *  the number the number of files and returns a vector<string>
 *  with the filenames.
 */ 

#define FRAME_FILES 10
#define DSC_FILES   20
#define IDX_FILES   30

#define __max_length_dosepix  10000
#define __max_length_timepix3 	1000
#define __max_length_dexter 	10000

class ListHandler {

public:
  ListHandler(Char_t *);
  ~ListHandler(){};
  int getListToLoop (TString tempScratchDir, vector<string> &, vector<string> &, vector<string> & );
  int getListToLoopDosepix (TString /*tempScratchDir*/, vector<string> & files);
  int getListToLoopTimepix3 (TString /*tempScratchDir*/, vector<string> & files);
  int getListToLoopDexterTXT(TString /*tempScratchDir*/, vector<string> & files);
  void list_order(vector<string> & );

  vector<string> getListToLoop(TString, Int_t typeF=FRAME_FILES);
  vector<string> OrderFilesPairMatch(vector<string>, vector<string>);

private:
  Long_t m_nFiles;
  string m_dirPath;
  vector<string> vectorOfFiles;
};
