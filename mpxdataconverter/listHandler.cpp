/*
 *      Copyright 2006 John Idarraga
 *
 *      This file is part of MAFalda.

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

#include "listHandler.h"
#include "algorithm"

ListHandler::ListHandler (Char_t * dirPath) {

	m_dirPath = dirPath;
	std::cout << "[INFO] will search inside " << m_dirPath << std::endl;

}

std::vector<std::string> ListHandler::OrderFilesPairMatch(std::vector<std::string> listFrames, std::vector<std::string> listDSC){


	std::vector<std::string> orderedString;

	std::vector<std::string>::iterator itFrames = listFrames.begin();
	std::vector<std::string>::iterator itDSD = listDSC.begin();

	TString tempName = "";
	TString tempName2 = "";

	for( ; itFrames != listFrames.end() ; itFrames++ )
	{
		// first possible name
		tempName = ((TString)(*itFrames));
		tempName+=".dsc";

		// second possible name
		tempName2 = ((TString)(*itFrames));
		int size = tempName2.Sizeof();
		bool marktobreak = false;
		for(int i = size-1 ; i >= 0 ; i--) {
			if(tempName2[i] == '.'){
				marktobreak = true;
			}
			// remove first
			tempName2.Remove(i);
			// then break if necessary
			if(marktobreak) break;
		}
		tempName2+=".dsc";

		/*
		if(((TString)(*itFrames)).Contains(".txt")) {
			tempName2.Remove(TString::kTrailing, 't');
			tempName2.Remove(TString::kTrailing, 'x');
			tempName2.Remove(TString::kTrailing, 't');
			tempName2.Remove(TString::kTrailing, '.');
		}
		tempName2+=".dsc";

		// third possible name
		tempName2 = ((TString)(*itFrames));
		if(((TString)(*itFrames)).Contains(".dat")) {
			tempName2.Remove(TString::kTrailing, 't');
			tempName2.Remove(TString::kTrailing, 'x');
			tempName2.Remove(TString::kTrailing, 't');
			tempName2.Remove(TString::kTrailing, '.');
		}
		tempName2+=".dsc";
		 */

		/* look for the corresponding .dsc file
	 	 that hopefully is not so far in the list.
	 	 No need to optimize this now, is fast enough.
		 */

		for(itDSD = listDSC.begin() ; itDSD != listDSC.end() ; itDSD++ ) {
			if((TString)(*itDSD) == tempName || (TString)(*itDSD) == tempName2) {
				/* get the good filename in the list and erase the entry
		 	 	   in  listDSC.  So it gets faster next time. */
				orderedString.push_back((*itDSD));
				listDSC.erase(itDSD);
				//std::cout << " !!! new size = " << listDSC.size() << std::endl;
				//std::cout << " !!! found --> " << (*itDSD) << std::endl;
				break;
			}
		}

	}

	return orderedString;
}

void ListHandler::list_order(vector<string> & f){

	sort (f.begin(), f.end());

}

int ListHandler::getListToLoopTimepix3 (TString /*tempScratchDir*/, vector<string> & files) {

	DIR * dp;
	struct dirent *dirp;
	if((dp  = opendir(m_dirPath.c_str())) == NULL) {
		cout << "Error(" << errno << ") opening " << m_dirPath << endl;
		return errno;
	}

	string tempfn;
	int nDat = 0;
	while ((dirp = readdir(dp)) != NULL) {
		tempfn = string(dirp->d_name);
		if( tempfn.find(".dat") != string::npos ) { // find number of csv files
			nDat++;
		}
	}
	// rewind
	rewinddir(dp);

	cout << "[INFO] found " << nDat << " dat files" << endl;

	// Prepare to retrieve all files
	string dat_string = ".dat";

	// Find an get dsc files now
	while ((dirp = readdir(dp)) != NULL) {

		tempfn = string(dirp->d_name);

		if( tempfn.find( dat_string.c_str() ) != string::npos ) { // found .dsc.dsc
			files.push_back( m_dirPath + "/" + tempfn );

		}
	}

	// These files are not in order.  Attempt to oder
	list_order(files);

    return nDat;
}

int ListHandler::getListToLoopDexterTXT(TString, vector<string> &files)
{

    DIR * dp;
    struct dirent *dirp;
    if((dp  = opendir(m_dirPath.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << m_dirPath << endl;
        return errno;
    }

    // Dexter files come with a Json file containing the metadata
    string tempfn;
    int nJson = 0;
    while ((dirp = readdir(dp)) != NULL) {
        tempfn = string(dirp->d_name);
        if( tempfn.find(".json") != string::npos ) { // find number of csv files
            nJson++;
        }
    }
    // rewind
    rewinddir(dp);

    cout << "[INFO] found " << nJson << " json Dexter files" << endl;

    // In Dexter there won't be one json file per frame.  It is not
    //  necesary.  If the conditions didn't change then only one
    //  json file will be present.

    // We need to pickup the txt files
    string csv_string = ".txt";

    // Find an get dsc files now
    while ((dirp = readdir(dp)) != NULL) {

        tempfn = string(dirp->d_name);

        if( tempfn.find( csv_string.c_str() ) != string::npos ) { // found .dsc.dsc
            files.push_back( m_dirPath + "/" + tempfn );
        }
    }

    // These files are not in order.  Attempt to oder
    list_order(files);

    return (int)files.size();
}

int ListHandler::getListToLoopDosepix (TString /*tempScratchDir*/, vector<string> & files) {

	DIR * dp;
	struct dirent *dirp;
	if((dp  = opendir(m_dirPath.c_str())) == NULL) {
		cout << "Error(" << errno << ") opening " << m_dirPath << endl;
		return errno;
	}

	string tempfn;
	int nCSV = 0;
	while ((dirp = readdir(dp)) != NULL) {
		tempfn = string(dirp->d_name);
		if( tempfn.find(".csv") != string::npos ) { // find number of csv files
			nCSV++;
		}
	}
	// rewind
	rewinddir(dp);

	cout << "[INFO] found " << nCSV << " dosepix files" << endl;

	// Prepare to retrieve all files
	string csv_string = ".csv";

	// Find an get dsc files now
	while ((dirp = readdir(dp)) != NULL) {

		tempfn = string(dirp->d_name);

		if( tempfn.find( csv_string.c_str() ) != string::npos ) { // found .dsc.dsc
			files.push_back( m_dirPath + "/" + tempfn );

		}
	}

	// These files are not in order.  Attempt to oder
	list_order(files);

	//vector<string>::iterator i = files.begin();
	//for ( ; i != files.end() ; i++ ) {
	//	cout << (*i) << endl;
	//}

	return nCSV;
}

int ListHandler::getListToLoop (TString /*tempScratchDir*/, vector<string> & files, vector<string> & dscfiles, vector<string> & idxfiles) {

	string tempfilename("lf_temp.txt");

	// If the user called the data files .dsc then the
	// dsc files will have extension .dsc.dsc
	// see if this happened
	DIR * dp;
	struct dirent *dirp;
	if((dp  = opendir(m_dirPath.c_str())) == NULL) {
		cout << "Error(" << errno << ") opening " << m_dirPath << endl;
		return errno;
	}

	string tempfn;
	bool doubledsc = false;
	while ((dirp = readdir(dp)) != NULL) {
		tempfn = string(dirp->d_name);
		if( tempfn.find(".dsc.dsc") != string::npos ) { // found .dsc.dsc
			doubledsc = true;
			break;
		}
	}
	// rewind
	rewinddir(dp);

	// See if there are .idx files (only present in multiframe)
	bool idxpresent = false;
	while ((dirp = readdir(dp)) != NULL) {
		tempfn = string(dirp->d_name);
		if( tempfn.find(".idx") != string::npos ) { // found .dsc.dsc
			idxpresent = true;
			break;
		}
	}
	// rewind
	rewinddir(dp);

	// Prepare to retrieve all files
	string dsc_string = ".dsc";
	if(doubledsc) {
		cout << "[WARNING] the dsc files have extension .dsc.dsc" << endl;
		dsc_string = ".dsc.dsc";
	}

	// Find an get dsc files now
	while ((dirp = readdir(dp)) != NULL) {

		tempfn = string(dirp->d_name);

		if( tempfn.find( dsc_string.c_str() ) != string::npos ) { // found .dsc.dsc
			dscfiles.push_back( m_dirPath + "/" + tempfn );

		}
	}

	// These files are not in order.  Attempt to oder
	list_order(dscfiles);

	// Using the dsc list build the rest
	vector<string>::iterator i = dscfiles.begin();
	string tempstr;
	size_t pos;
	for ( ; i != dscfiles.end() ; i++) {

		tempstr = *i;
		// data file
		pos = tempstr.find_last_not_of( dsc_string.c_str() );
		//cout << tempstr << " --> " << pos << endl;
		if(doubledsc) { // particular case
			tempstr = tempstr.substr(0, pos+2);
			tempstr += ".dsc";
		} else {
			tempstr = tempstr.substr(0, pos+1);
		}
		files.push_back( tempstr );


		// idx file
		if(idxpresent) {
			tempstr = *i;
			// data file
			pos = tempstr.find_last_not_of( dsc_string.c_str() );
			if(doubledsc) {
				tempstr = tempstr.substr(0, pos+2);
				tempstr += ".dsc"; // particular case
			} else {
				tempstr = tempstr.substr(0, pos+1);
			}
			tempstr += ".idx";
			idxfiles.push_back( tempstr );
		}

	}

	// close dir
	closedir(dp);

	return (int) files.size();
}



std::vector<std::string> ListHandler::getListToLoop(TString tempScratchDir, Int_t typeF){

	TString filename = "";
	if(tempScratchDir.Length() > 0){
		filename += tempScratchDir;
		filename += "/";
	}

	if (typeF==FRAME_FILES) filename += "listOfFiles.txt";
	if (typeF==DSC_FILES) filename += "listOfFiles.dsc.txt";
	if (typeF==IDX_FILES) filename += "listOfFiles.idx.txt";

	// Command to erase old listing files
	std::string command0 = "/bin/bash -c 'rm -f "; // listOfFiles.txt"   "'";
	command0 += filename.Data();
	command0 += "'";

	std::string command02 = "/bin/bash -c 'rm -f ";
	command02 += filename.Data();
	command02 += "'";

	std::string command03 = "/bin/bash -c 'rm -f ";
	command03 += filename.Data();
	command03 += "'";

	// Commands to find the number of related files
	std::string commandn = "/bin/bash -c \"ls ";
	std::string commandn2 = "/bin/bash -c 'ls ";
	std::string commandn3 = "/bin/bash -c 'ls ";

	// Commands to find the list of related files
	std::string command1 = "/bin/bash -c 'for a in ";
	std::string command2 = "/bin/bash -c 'for a in ";
	std::string command3 = "/bin/bash -c 'for a in ";

	// Number of frame files, i.e. not dsc or idx
	commandn += "'" + m_dirPath + "'" + " | awk '{ print !(/.dsc/||/.idx/) }' | grep 1 | wc -l >> "; //listOfFiles.txt'";
	commandn += filename.Data();
	commandn += " \"";

	commandn2 += "\"" + m_dirPath + "\"" + " | grep \".dsc\" | wc -l >> "; //listOfFiles.dsc.txt'";
	commandn2 += filename.Data();
	commandn2 += "'";

	commandn3 += "\"" + m_dirPath + "\"" + " | grep \".idx\" | wc -l >> "; //listOfFiles.idx.txt'";
	commandn3 += filename.Data();
	commandn3 += "'";

	// Get the files that are not dsc or idx, sometimes frame data comes
	// Without extension (.txt), or can be any other
	command1 += "\"" + m_dirPath + "\"" + "/* ; do [[ $a == *.dsc || $a == *.idx ]] || echo $a >> "; // listOfFiles.txt ; done'";
	command1 += filename.Data();
	command1 += " ; done'";

	// Get the files that are dsc
	command2 += "\"" + m_dirPath + "\"" + "/* ; do [[ $a == *.dsc ]] && echo $a >> "; // listOfFiles.dsc.txt ; done'";
	command2 += filename.Data();
	command2 += " ; done'";

	// Get the files that are idx
	command3 += "\"" + m_dirPath + "\"" + "/* ; do [[ $a == *.idx ]] && echo $a >> "; // listOfFiles.dsc.txt ; done'";
	command3 += filename.Data();
	command3 += " ; done'";


	/* get old files erased */
	if(typeF==FRAME_FILES) system(command0.c_str());
	if(typeF==DSC_FILES) system(command02.c_str());
	if(typeF==IDX_FILES) system(command03.c_str());

	if(typeF==FRAME_FILES) system(commandn.c_str());
	if(typeF==DSC_FILES) system(commandn2.c_str());
	if(typeF==IDX_FILES) system(commandn3.c_str());

	int rcm = 0;
	//rcm2 = 0, rcm3 = 0;
	if(typeF==FRAME_FILES) rcm = system(command1.c_str());
	//if(typeF==DSC_FILES) rcm2 = system(command2.c_str());
	//if(typeF==IDX_FILES) rcm3 = system(command3.c_str());

	bool throwerror = false;
	std::vector<std::string> listOfFiles; // list to return

	if (typeF==FRAME_FILES && rcm != 0) {
		throwerror = true;
	}
	//else if (typeF==DSC_FILES && rcm2 != 0){
	//throwerror = true;
	//}
	//std::cout << command2 << " : " << rcm << "," << rcm2 << std::endl;

	if(throwerror){
		std::cout << "[ERROR] unexpected error when listing files inside " << m_dirPath << std::endl;
		std::cout << "        please READ all the bash errors above. " << std::endl;
		std::cout << "        Other possible reasons: " << std::endl;
		std::cout << "        I need to make temporary files.  May be you " << std::endl;
		std::cout << "        don't have 'w' permissions in this directory ? " << std::endl;
		std::cout << "        Use the third (optional) parameter to this " << std::endl;
		std::cout << "        program to change the temp dir" << std::endl;
		return listOfFiles; // empty
	}

	fstream filestrList;
	filestrList.open(filename, fstream::in);

	std::string lineTemp;
	Int_t cntrInFile = 0;
	Int_t cntrInActualFiles = 0;
	Int_t nFiles = 0;

	while (filestrList.good()) {

		filestrList >> lineTemp;

		if(cntrInFile == 0)
		{
			nFiles = atoi(lineTemp.c_str());
		}
		if(cntrInFile > 0 && cntrInFile <= nFiles)
		{
			//std::cout << lineTemp << std::endl;
			listOfFiles.push_back(lineTemp);
			cntrInActualFiles++;
		}
		cntrInFile++;
	}

	if(nFiles != cntrInActualFiles){
		std::cout << "[ERROR] There should be " << nFiles << " frame (and .dsc) files in this directory and I found "
				<< cntrInFile << "... giving up.";
		exit(1);
	}

	return listOfFiles;
}
