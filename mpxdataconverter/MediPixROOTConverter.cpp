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

#include <iostream>
#include <string>
#include <vector>

#include "TH2.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TColor.h"
#include "TPaletteAxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPaveLabel.h"
#include "TControlBar.h"
#include "TObject.h"
#include "TString.h"

#include "MediPixViewer.h"
#include "listHandler.h"
#include "control.h"
#include "MediPixWriteToEntuple.h"
#include "allpix_dm.h"

using namespace std;

void checkParameters(int, char**, TString *);
void histoProperties(TH2I *, TCanvas *);
void ProcessDosePixFiles(vector<string>, Long_t, FramesHandler *, WriteToNtuple *);
bool ProcessOneDosePixFile(string fn, FramesHandler *, WriteToNtuple * );
void ProcessTimepix3Files(vector<string>, FramesHandler * frames, WriteToNtuple * ntup);
void ProcessDexterFiles(vector<string> datfiles, FramesHandler * frames, WriteToNtuple * ntup);
bool ProcessOneTimepix3File(string fn, FramesHandler *, WriteToNtuple * );
bool ProcessOneDexterFile(string fn, FramesHandler *, WriteToNtuple * );

TApplication * g_theApp = nullptr;
Int_t g_direction = 0;

int main(int argc, char ** argv){
	cout << "[START] Main entry point\n";
	TApplication * g_theApp = new TApplication("Output", 0, NULL);

	/////////////////////////////////////
	// check flags
	TString tempScratchDir = "";
	checkParameters(argc, argv, &tempScratchDir);

	cout << "[DEBUG] Flags checked\n";

	// create ntuple and FramesHandler
	TString dataset = argv[2];
	cout << "[DEBUG] Dataset = " << dataset << "\n";

	WriteToNtuple * MPXnTuple = new WriteToNtuple(dataset, tempScratchDir);
	cout << "[DEBUG] Made WriteToNTuple object\n";

	FramesHandler frames(dataset);
	cout << "[DEBUG] Made FrameHandler for the dataset\n";

	Long_t skipFrames = 0;
	if(argc == 5) skipFrames = atoi(argv[4]);

	/////////////////////////////////////
	// Handle the list of files
	ListHandler handlerL(argv[1]);

	vector<string> listOfFiles;
	vector<string> listOfDSCFiles;
	vector<string> listOfIDXFiles;

	vector<string>::iterator itrF;
	int nPixelmanFiles = handlerL.getListToLoop (tempScratchDir, listOfFiles, listOfDSCFiles, listOfIDXFiles);

	cout << "[DEBUG] Handled list of files. Starting search\n";

	// Search for Timepix3 files
	int nTimepix3 = 0;
	if ( nPixelmanFiles == 0 ) {
		nTimepix3 = handlerL.getListToLoopTimepix3(tempScratchDir, listOfFiles);
	}

	int nDosepixFiles = 0;
	if( nPixelmanFiles == 0 && nTimepix3 == 0) { // Search for Dosepix files
        nDosepixFiles = handlerL.getListToLoopDosepix(tempScratchDir, listOfFiles);
	}

	int nDexterFiles = 0;
	if ( nPixelmanFiles == 0 && nTimepix3 == 0 && nDosepixFiles == 0) {
	        nDexterFiles = handlerL.getListToLoopDexterTXT(tempScratchDir, listOfFiles);
	}

	// Handle dosepix files if any
	if ( nDosepixFiles != 0 ) ProcessDosePixFiles(listOfFiles, skipFrames, &frames, MPXnTuple);

	// Handle Timepix3 files if any
	if ( nTimepix3 != 0 ) ProcessTimepix3Files(listOfFiles, &frames, MPXnTuple);

       // Handle Dexter files
       if ( nDexterFiles != 0 ) ProcessDexterFiles(listOfFiles, &frames, MPXnTuple);

	// If no files at all
	if( nPixelmanFiles==0 && nDosepixFiles==0 && nTimepix3==0 && nDexterFiles==0) {
	// Nothing at all
        cout << "[NONE] Could not find any files. Nothing to be done." << endl;
	delete MPXnTuple;
	return 0;
	}

	// If only Timepix3
	if ( nTimepix3 !=0 ) {
		MPXnTuple->closeNtuple();
		cout << "[ OK ] done with the conversion.  Check your output file --> "
				<< MPXnTuple->GetNtupleFileName() << endl;
		return 0;
	}

	cout << "[INFO] there are " << listOfDSCFiles.size() << " DSC files in directory" << endl;
	cout << "[INFO] there are " << listOfFiles.size() << " frame files in directory" << endl;
	cout << "[INFO] there are " << listOfIDXFiles.size() << " IDX files in directory" << endl;

	if ( listOfDSCFiles.empty() && nDexterFiles == 0) {
		cout << "[ERROR] I can't find the DSC files in the directory. Giving up." << endl;
		delete MPXnTuple;
		return 1;
	}

	if ( listOfFiles.empty() ) {
		cout << "[ERROR] I can't find some files in the specified directory or the dir does not exist at all." << endl;
		delete MPXnTuple;
		return 1;
	}

	// make sure frames list and DSC list match alright
	//if(listOfDSCFiles.size() != 0) {
	//listOfDSCFiles = handlerL.OrderFilesPairMatch(listOfFiles, listOfDSCFiles);
	//}

	//
	if( (listOfFiles.size() != listOfDSCFiles.size()) && nDexterFiles == 0 ) {
		std::cout << "[ERROR] you have a different amount of dsc and frame files." << std::endl;
		std::cout << "        You probably have an extra file that doesn't correspond" << std::endl;
		std::cout << "        to the data taking ... check it out.  Giving up." << std::endl;
		exit(1);
	}

	Long_t filesItr;
	string oneFileName = "";
	string oneDSCFileName = "";
	string oneIDXFileName = "";

	int ftype;

	for(filesItr = 0 ; filesItr < (Int_t)listOfFiles.size() ; filesItr++) {

		if(filesItr < skipFrames) continue; // skip if needed

		oneFileName = listOfFiles[filesItr];
		oneDSCFileName = listOfDSCFiles[filesItr];
		if(listOfIDXFiles.size() > 0) oneIDXFileName = listOfIDXFiles[filesItr];

        //frames.getAFrameMatrix((TString) oneFileName, (TString) oneDSCFileName);
		int nread = frames.readOneFrame((TString) oneFileName, (TString) oneDSCFileName, &ftype);

		// Write what was read
		if(nread == 1) {
			MPXnTuple->fillVars(&frames); // rewind metadata for the next frame
		}

		// Deal with multi-frame here
		if(nread > 1) {
			frames.ProcessMultiframe(oneFileName, oneDSCFileName, oneIDXFileName, MPXnTuple, ftype);
			// rewind the frame at the very end in this case.  Meta-data wasn't rewinded before
			frames.RewindAll(false);
		}

	}

	MPXnTuple->closeNtuple();

	cout << "[ OK ] done with the conversion.  Check your output file --> "
			<< MPXnTuple->GetNtupleFileName() << endl;

	/*
	// erasing possible remaining working files
	std::string command0 = "/bin/bash -c 'if [ -e ";
	if(tempScratchDir.Length() > 0) { command0 += tempScratchDir.Data(); command0 += "/"; }
	command0 += "listOfFiles.txt ] ; then rm -f ";
	if(tempScratchDir.Length() > 0) { command0 += tempScratchDir.Data(); command0 += "/"; }
	command0 += "listOfFiles.txt ; fi'";

	std::string command02 = "/bin/bash -c 'if [ -e ";
	if(tempScratchDir.Length() > 0) { command02 += tempScratchDir.Data(); command02 += "/"; }
	command02 += "listOfFiles.dsc.txt ] ; then rm -f ";
	if(tempScratchDir.Length() > 0) { command02 += tempScratchDir.Data(); command02 += "/"; }
	command02 += "listOfFiles.dsc.txt ; fi'";

	std::string command03 = "/bin/bash -c 'if [ -e ";
	if(tempScratchDir.Length() > 0) { command03 += tempScratchDir.Data(); command03 += "/"; }
	command03 += "listOfFiles.idx.txt ] ; then rm -f ";
	if(tempScratchDir.Length() > 0) { command03 += tempScratchDir.Data(); command03 += "/"; }
	command03 += "listOfFiles.idx.txt ; fi'";


	system(command0.c_str());
	system(command02.c_str());
	system(command03.c_str());
	 */
	return 0;
}

void checkParameters(int argc, char ** argv, TString * tempScratchDir){

	if(argc < 3)
	{
		std::cout << "  use:" << std::endl;
		std::cout << "      " << argv[0] << " pathToData  dataSetName  tempScratchDir(optional.  ./ default)  skip(optional)" << std::endl;
		std::cout << "       pathToData   : Directory where your have stored your data " << std::endl;
		std::cout << "                      (uncompressed). Can be an absolute or " << std::endl;
		std::cout << "                      relative path." << std::endl;
		std::cout << "       dataSetName  : A short string to easily identify your data later on." << std::endl;
		std::cout << "                      Example --> MediPix_SPS_TOTmode_25-06-2007" << std::endl;
		std::cout << "                      please do not use special character or spaces." << std::endl;
		std::cout << "       tempFilesDir : In case you can only read from the directory" << std::endl;
		std::cout << "                      containing the data but can't write, you can" << std::endl;
		std::cout << "                      tell " << argv[0] << "to use a different temp dir."<< std::endl;
		std::cout << "       skip         : Number of events to skip optional (int)" << std::endl;

		exit(1);
	}

	TString dataSet = argv[2];

	if (dataSet == "" || dataSet.Contains(' '))
	{
		std::cout << "[ERROR] badformat in dataSet !" << std::endl;
		std::cout << "        dataSetNumber: example --> MediPix_SPS_TOTmode_25-06-2007" << std::endl;
		std::cout << "        please do not use special character or spaces." << std::endl;
		exit(1);
	}

	// user decided to use a different scratch dir
	// optional
	if ( argc > 2 ) {
		*tempScratchDir += argv[3];
		if(tempScratchDir->Length() == 0) *tempScratchDir = "./";
		cout << "[INFO] working directory is " << *tempScratchDir << endl;
	}

}


void histoProperties(TH2I * h1, TCanvas * c){

	h1->SetStats(false);
	gStyle->SetPalette(1,0);

	c->SetBorderMode(0);
	c->SetFillColor(kWhite);

	c->Update();

}

void ProcessDexterFiles(vector<string> datfiles, FramesHandler * frames, WriteToNtuple * ntup) {

    cout << "[ OK ] Process Dexter files, first filename:" << datfiles[0] << "\n";
    Long_t filesItr;
    string oneFileName;

    for(filesItr = 0 ; filesItr < (Long_t) datfiles.size() ; filesItr++) {

        oneFileName = datfiles[filesItr];
        // Parse one file
        ProcessOneDexterFile(oneFileName, frames, ntup);
    }

}

void ProcessTimepix3Files(vector<string> datfiles, FramesHandler * frames, WriteToNtuple * ntup) {

	Long_t filesItr;
	string oneFileName;

	for(filesItr = 0 ; filesItr < (Long_t) datfiles.size() ; filesItr++) {

		oneFileName = datfiles[filesItr];
		// Parse one file
		ProcessOneTimepix3File(oneFileName, frames, ntup);
	}

}


void ProcessDosePixFiles(vector<string> csvfiles, Long_t skipFrames, FramesHandler * frames, WriteToNtuple * ntup) {

	Long_t filesItr;
	string oneFileName;

	for(filesItr = 0 ; filesItr < (Long_t) csvfiles.size() ; filesItr++) {

		if(filesItr < skipFrames) continue; // skip if needed

		oneFileName = csvfiles[filesItr];

		// open one file and start extraction
		ProcessOneDosePixFile(oneFileName, frames, ntup);

	}

	ntup->closeNtuple();
	cout << "[ OK ] done with the conversion.  Check your output file --> " << ntup->GetNtupleFileName() << endl;

}

bool ProcessOneDexterFile(string fn, FramesHandler *, WriteToNtuple * ) {

    std::filebuf fb;
    fb.open (fn.c_str(), std::ios::in);
    if( ! fb.is_open() ) {
        cout << "[ERROR] Couldn't open file " << fn.c_str() << endl;
        return false;
    } else {
        cout << "Opening file " << fn.c_str() << endl;
    }

    // Dexter clear text Parser MATRIX
    char temp[__max_length_timepix3];

    if ( fb.is_open() ) {

        std::istream is(&fb);
        while ( is.good() ) {

            // check end of file
            if ( is.eof() ) {
                fb.close();
                break;
            }

            // read a line
            is.getline(temp, __max_length_dexter);
            string toparse = string(temp);
            //cout << toparse << endl;

        }

    }

    fb.close();

}

bool ProcessOneTimepix3File(string fn, FramesHandler * frame, WriteToNtuple * ntup) {

	std::filebuf fb;
	fb.open (fn.c_str(), std::ios::in);
	if( ! fb.is_open() ) {
		cout << "[ERROR] Couldn't open file " << fn.c_str() << endl;
		return false;
	} else {
		cout << "Opening file " << fn.c_str() << endl;
	}

	// Timepix3 clear text Parser (.dat file)
	int stampcntr = 0, inlineCntr = 0;
	char temp[__max_length_timepix3];
	size_t pos;
	vector<unsigned int> oneStamp; // x, y, ToT, ToA16, ToA14, fToA, calculatedTime(double);
	string tempStr;
	bool endFrame = false;

	if ( fb.is_open() ) {

		std::istream is(&fb);
		while ( is.good() ) {

			// read a line
			is.getline(temp, __max_length_timepix3);
			string toparse = string(temp);

			// check end of file
			if ( is.eof() ) {
				fb.close();
				break;
			}

			// Initial comment. Do nothing
			if ( toparse[0] == '#' and stampcntr == 0 ) continue;

			// This is the end of a "frame"
			if ( toparse[0] == '#' and stampcntr != 0 ) {
				endFrame = true;
			}

			//cout << "[PARS] --> " << toparse << endl;

			// Rewind line
			inlineCntr = 0;
			while ( toparse.size() > 1 && !endFrame ) {

				// Find a tab
				pos = toparse.find_first_of( 0x9 ); // Tab

				// Otherwise find the end of the line
				if( pos == string::npos ) {
					pos = toparse.find_first_of( 0x0a );
				}
				if( pos == string::npos ) {
					pos = toparse.find_first_of( 0x0d );
				}
				// Otherwise this is the last value break
				if( pos == string::npos ) {
					pos = toparse.size() - 1;
				}

				// Get the substring
				tempStr = toparse.substr(0, pos+1);

				// shorten parsable string
				toparse.erase(0, pos+1);

				// The value is available here
				// Take x, y, ToT, ToA16, ToA14, fToA (6 values)
				if( inlineCntr > __timepix3_stampLength ) break;

				oneStamp.push_back( atoi ( tempStr.c_str() ) );
				inlineCntr++;

			}

			if ( endFrame ) {

				// At this point one frame has been read
				frame->readOneStampTimepix3( oneStamp );
				// send to the ntup
				ntup->fillVars(frame);

				// rewind
				oneStamp.clear();
				endFrame = false;
			}

			stampcntr++;

		}

	}

	cout << "[INFO] processed " << stampcntr << " stamps" << endl;

	return true;
}

bool ProcessOneDosePixFile(string fn, FramesHandler * frame, WriteToNtuple * ntup) {

	std::filebuf fb;
	fb.open (fn.c_str(), std::ios::in);
	if( ! fb.is_open() ) {
		cout << "[ERROR] Couldn't open file " << fn.c_str() << endl;
		return false;
	} else {
		cout << "Opening file " << fn.c_str() << endl;
	}
	char temp[__max_length_dosepix];
	size_t pos;
	string timestamp;
	string tempval;

	vector<int> oneframe;
	int framecntr = 0;

	if ( fb.is_open() ) {

		std::istream is(&fb);
		while ( is.good() ) {

			// read a line
			is.getline(temp, __max_length_dosepix);
			string toparse = string(temp);

			// check end of file
			if ( is.eof() ) {
				fb.close();
				break;
			}

			// Search for the timestamp.  Find first ','
			pos = toparse.find_first_of( ',' );
			timestamp = toparse.substr(0, pos);
			// Get rid of the timestamp
			toparse.erase(0, pos+1);

			bool foundcolon = true;

			while ( foundcolon ) {

				pos = toparse.find_first_of( ',' );
				// get the val
				tempval = toparse.substr(0, pos);
				// and remove it
				toparse.erase(0, pos+1);

				// push the value
				oneframe.push_back( atoi ( tempval.c_str() ) );

				//cout << tempval << " ";

				if ( pos == string::npos ) {  // end of the frame (same as end of the line)
					foundcolon = false;
				}

			}

			// At this point one frame has been read
			frame->readOneFrameDosepix(oneframe, timestamp);
			// send to the ntup
			ntup->fillVars(frame);

			// rewind
			oneframe.clear();

			framecntr++;
		}

	}

	cout << "[INFO] processed " << framecntr << " frames" << endl;

	return true;
}

