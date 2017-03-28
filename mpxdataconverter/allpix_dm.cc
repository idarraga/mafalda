/**
 *  Author: John Idarraga <idarraga@cern.ch>
 *  Classes containing the frames
 */

#include <fstream>
#include <istream>

#include <string>
#include "TH2.h"

#include "MediPixWriteToEntuple.h"
#include "allpix_dm.h"

#include "TObject.h"

using namespace std;


ClassImp(FrameContainer)
ClassImp(FrameStruct)

FrameContainer::FrameContainer(){

	m_nEntriesPad = 0;
	m_nHitsInPad = 0;
	m_nChargeInPad = 0;
	m_isMCData = false;

}

void FrameContainer::FillOneElement(Int_t xi, Int_t yi, Int_t width, Int_t counts) {

	// X,Y,C --> X,C : yi*width + xi
	Int_t X = yi*width + xi;

	m_frameXC[X] += counts;        // TOT or count(binary detector)

	// If the pixel didn't exist this is an extra entry
	if ( m_frameXC.find(X) == m_frameXC.end() ) m_nEntriesPad++;
	// But always an extra hit
	m_nHitsInPad++;
	// Increase the total counts
	m_nChargeInPad += counts;

}

void FrameContainer::FillOneElement(Int_t xi, Int_t yi, Int_t width, Int_t counts, Int_t ToA, Int_t fastToA) {

	// X,Y,C --> X,C : yi*width + xi
	Int_t X = yi*width + xi;

	m_frameXC[X] += counts;        // TOT or count(binary detector)
	m_frameXToA[X] += ToA;        // TOT or count(binary detector)
	m_frameXFastToA[X] += fastToA;        // TOT or count(binary detector)

	// If the pixel didn't exist this is an extra entry
	if ( m_frameXC.find(X) == m_frameXC.end() ) m_nEntriesPad++;
	// But always an extra hit
	m_nHitsInPad++;
	// Increase the total counts
	m_nChargeInPad += counts;

}


void FrameContainer::FillOneElement(Int_t xi, Int_t counts) {

	m_frameXC[xi] += counts;        // TOT or count(binary detector)

	// If the pixel didn't exist this is an extra entry
	if(m_frameXC.find(xi) == m_frameXC.end()) m_nEntriesPad++;
	// But always an extra hit
	m_nHitsInPad++;
	// Increase the total counts
	m_nChargeInPad += counts;

}

void FrameContainer::FillOneElement(Int_t xi, Int_t yi, Int_t width, Int_t counts,
		Double_t truthE, Double_t E){

	// Fill pixel withouth MC info first
	FillOneElement(xi, yi, width, counts);

	// X,Y,C --> X,C : yi*width + xi
	Int_t X = yi*width + xi;

	m_frameXC_TruthE[X] += truthE; // Truth energy
	m_frameXC_E[X] += E;           // Energy with detector effects

}

void FrameContainer::SetLVL1(Int_t xi, Int_t yi, Int_t width, Int_t lvl1){

	// X,Y,C --> X,C : yi*width + xi
	Int_t X = yi*width + xi;

	m_lvl1[X] += lvl1;
}

void FrameContainer::ResetCountersPad(){

	m_nEntriesPad = 0;
	m_nHitsInPad = 0;
	m_nChargeInPad = 0;

}

FrameStruct::FrameStruct(TString dataset) 
: FrameContainer() {

	// clean matrix, every pixel to 0
	CleanUpMatrix();
	// reset counters m_nEntriesPad, m_nHitsInPad, m_nChargeInPad
	ResetCountersPad();
	// reset metadata
	RewindMetaDataValues();

	/* joint */
	fFrameId = -1;
	fMPXDataSetNumber = dataset;

}

void FrameStruct::RewindMetaDataValues(){

	/* head info */
	fFormat = 0;
	fWidth = 0;
	fHeight = 0;

	/* all metadata */
	fAcq_mode = 0;
	fAcq_time = 0.0;
	fApplied_filters = "";
	fAuto_erase_interval = 0.0;
	fAutoerase_interval_counter = 0;
	fBS_active = false;
	fChipboardID = "";
	fCoinc_live_time = 0.0;
	fCoincidence_delay = 0;
	fCoincidence_mode = 0;
	fCounters.clear();
	fDACs.clear();
	fHV = 0.0;
	fHw_timer = 0;
	fInterface = "";
	fMpx_clock = 0.0;
	fMpx_type = -1;
	fPolarity = -1;
	fStart_time = 0.0;
	fStart_timeS = "";
	fTimepix_clock = 0;
	fTrigger_time = 0.0;

	/* Joint data doesn't need to be rewinded
	 *  --> fFrameId
	 *  --> fMPXDataSetNumber
	 */
}

void FrameStruct::SetPrimaryVertex(Double_t vx, Double_t vy, Double_t vz){

	m_primaryVertex_x.push_back(vx);
	m_primaryVertex_y.push_back(vy);
	m_primaryVertex_z.push_back(vz);

}

void FrameStruct::FillMetaData(TString METAString, Int_t metaCode){

	TString tempEntry = "";

	// temporary values to be used
	int temp_i;
	double temp_d;

	switch (metaCode)
	{
	case ACQ_MODE_CODE:
		fAcq_mode = METAString.Atoi();
		break;
	case ACQ_TIME_CODE:
		fAcq_time = METAString.Atof();
		break;
	case CHIPBOARDID_CODE:
		fChipboardID = METAString;
		break;
	case DACS_CODE:
		for(Int_t i = 0 ; i < METAString.Length() ; i++)
		{
			if(METAString[i] != ' ')
			{
				tempEntry += METAString[i];
			}
			else if(METAString[i] == ' ' && tempEntry.Length() > 0)
			{
				//if(fDACsSize < MAX_DAQ_INT){
				//std::cout << tempEntry << std::endl;
				//fDACs[fDACsSize] = tempEntry.Atoi();
				//fDACsSize++;
				fDACs.push_back((Int_t)tempEntry.Atoi());
				//}
				tempEntry = "";
			}
		}
		tempEntry = "";
		break;
	case HW_TIMER_CODE:
		fHw_timer = METAString.Atoi();
		break;
	case INTERFACE_CODE:
		fInterface = METAString;
		break;
	case POLARITY_CODE:
		fPolarity = METAString.Atoi();
		break;
	case START_TIME_CODE:
		fStart_time = METAString.Atof();
		break;
	case START_TIME_S_CODE:
		fStart_timeS = METAString;
		break;
	case MPX_CLOCK_CODE:
		fMpx_clock = METAString.Atof();
		break;
	case HV_CODE:
		fHV = METAString.Atof();
		break;
	case TIMEPIX_CLOCK_CODE:
		temp_i = METAString.Atoi();
		temp_d = 10.0; // MHz
		switch (temp_i) {
		case 0:
			temp_d = 10.0;
			break;
		case 1:
			temp_d = 20.0;
			break;
		case 2:
			temp_d = 40.0;
			break;
		case 3:
			temp_d = 80.0;
			break;
		default:
			std::cout << "[WARNING] couldn't determine de timepix clock : " << endl;
			cout << METAString << endl;
			break;
		};
		fTimepix_clock = temp_d;
		break;
		case TIMEPIX_CLOCKINMHZ_CODE:
			fTimepix_clock = METAString.Atof();
			break;
		default:
			std::cout << "[WARNING] couldn't parse MetaData Info" << std::endl;
			break;
	}

	// FIXME ! .... have to decide how to give real dataSet numbers to data
	// PGP signature ?, md5 of the output file ? user decided ?
	// connect to db ?

	fMPXDataSetNumber = "1234ABCD";

}


FramesHandler::FramesHandler(TString dataset){

	//////////////////////////////////////////
	// Instanciate one frame
	// is going to be overwritten many times
	// but instanciated only once !
	m_aFrame = new FrameStruct(dataset);
	m_nFrames = 0;

	getAFrameMatrix_flag = false;
	getAFrameHist_flag = false;

	m_metaBit = -1;
	m_metaCode = -1;
	m_ParseAndHold = true;

	nFrames256x256 = 0;
	nFramesXYC = 0;

}

FramesHandler::~FramesHandler(){

	delete m_aFrame;
}

/*
TH2I * FramesHandler::getHistFrame(Int_t frameId, Int_t * frameMatrix){

	TString title = "frame ", hname = "frame_";
	title += frameId; hname += frameId;
	TH2I * histFrame = new TH2I(hname, title,
			MAX_FRAME_ROW, 0, MAX_FRAME_ROW-1,
			MAX_FRAME_COL, 0, MAX_FRAME_COL-1
	);

	Int_t iCntr = 0;
	Int_t jCntr = 0;
	for(Int_t i = 0 ; i < MAX_FRAME_ROW*MAX_FRAME_COL ; i++)
	{
		histFrame->Fill(jCntr, iCntr, frameMatrix[i]);
		iCntr++;
		if(iCntr >= MAX_FRAME_ROW)
		{
			iCntr = 0;
			jCntr++;
		}
		//std::cout << frameMatrix[i][j] << " ";
	}

	return histFrame;
}
 */

Int_t FramesHandler::IdentifyTypeOfInput2(TString fullDSCFileName, Int_t & width, Int_t & height, int & nFrames){

	// Find enconding
	int encoding = -1;

	fstream filestr;
	filestr.open(fullDSCFileName, fstream::in);

	Char_t temp[2048];
	string tempS = "";

	// Extract the very first line which is of the form
	// B000000X or A0000000X where B:=binary, A:=ASCII
	// and X is the number of frames associated to this
	// description file
	if(filestr.good()){
		filestr.getline(temp, 2048);
	}
	if(temp[0] == 'A') { // Ascii
		encoding = FSAVE_ASCII;
	}else if(temp[0] == 'B') { // Binary
		encoding = FSAVE_BINARY;
	}
	if(encoding == -1) return 0; // no data asociated to this dsc file (error)

	// Extract number of frames to process
	temp[0] = '0'; // replace the character and check the number of files
	nFrames = atoi(temp);
	//cout << "Number of frames associated to this dsc file = " << nFrames << endl;

	string infoS;
	string fSave = "";
	size_t found = std::string::npos;

	while (filestr.good())// && lineCntr < 1)
	{

		tempS = temp;

		// skip one line and extract the info in the third line
		filestr.getline(temp, 2048); // 2nd line
		filestr.getline(temp, 2048); // 3rd line
		tempS = temp;

		// At this point I need to analyze this type of string
		// Type=i16 [X,Y,C] width=1024 height=512
		infoS = WIDTH_STRING;
		found = tempS.find(infoS);
		// extract the width number
		int i = (int)found + (int)infoS.length();
		while ( tempS[i] != 0x20 && ( tempS[i] != 0xd && tempS[i] != '\0' ) )
		{
			if( tempS[i] != 0x20 )
			{
				fSave.append(1, tempS[i]);
			}
			i++;
		}

		// width ready
		width = atoi(fSave.c_str());
		m_width = width;

		fSave.clear();

		// extract the height number
		infoS = HEIGHT_STRING;
		found = tempS.find(infoS);
		//tempS[i] != 0xd &&
		i = (int) found + (int) infoS.length();
		while ( tempS[i] != 0x20 && ( tempS[i] != 0xd && tempS[i] != '\0' ) )
		{
			if ( tempS[i] != 0x20 )
			{
				fSave.append( 1, tempS[i] );
			}
			i++;
		}

		// height ready
		height = atoi(fSave.c_str());
		m_height = height;

		// Now extract the format
		// X,Y,C
		infoS = NEW_STRING_TYPEI16_XYC;
		if ( (found = tempS.find(infoS)) != std::string::npos) {
			return encoding | FSAVE_I16 | FSAVE_SPARSEXY;
		}

		infoS = NEW_STRING_TYPEU32_XYC;
		if ( (found = tempS.find(infoS)) != std::string::npos) {
			return encoding | FSAVE_U32 | FSAVE_SPARSEXY;
		}

		// Matrix
		infoS = NEW_STRING_TYPEI16_MATRIX;
		if ( (found = tempS.find(infoS)) != std::string::npos) {
			return encoding | FSAVE_I16;
		}

		infoS = NEW_STRING_TYPEU32_MATRIX;
		if ( (found = tempS.find(infoS)) != std::string::npos) {
			return encoding | FSAVE_U32;
		}

		// X,C
		infoS = NEW_STRING_TYPEI16_XC;
		if ( (found = tempS.find(infoS)) != std::string::npos) {
			return encoding | FSAVE_I16 | FSAVE_SPARSEX;
		}

		infoS = NEW_STRING_TYPEU32_XC;
		if ( (found = tempS.find(infoS)) != std::string::npos) {
			return encoding | FSAVE_I16 | FSAVE_SPARSEX;
		}

		break;

	}

	// otherwise, meaning number of frames to process.  0 is an error.
	return 0;
}

Int_t FramesHandler::IdentifyTypeOfInput(TString fullDSCFileName, Int_t & /*width*/, Int_t & /*height*/, int & /*nframes*/){

	fstream filestr;
	filestr.open(fullDSCFileName, fstream::in);

	Int_t lineCntr = 0;
	Char_t temp[2048];
	std::string tempS = "";
	std::string typeS = FRAME_TYPE_STRING;
	std::string binTypeS = BINARY_FRAME_TYPE_STRING;
	std::string fSave = "";
	size_t found = std::string::npos;

	while (filestr.good())// && lineCntr < 1)
	{

		filestr.getline(temp, 2048);
		tempS = temp;

		found = tempS.find(typeS);
		if (found != std::string::npos)
		{

			// found the format line, now pick up the format BITS
			for(int i = (int)typeS.length() ; i < (int)tempS.length() ; i++)
			{
				if(tempS[i] != ' ')
				{
					fSave.append(1,tempS[i]);
				}
				else
					break;
			}
		}

		/* check if binary */
		found = tempS.find(binTypeS);
		if (found != std::string::npos) {
			return FSAVE_BINARY;
		}

		lineCntr++;
	}

	filestr.close();

	return atoi(fSave.c_str());
}

void FrameContainer::CleanUpMatrix(){
	m_frameXC.clear();
	m_frameXToA.clear();
	m_frameXFastToA.clear();
	m_lvl1.clear();
	m_frameXC_TruthE.clear();
	m_frameXC_E.clear();
}

/* Rewind frame and metadata */
void FramesHandler::RewindAll(bool rewind_metadata){

	// clean matrix, every pixel to 0
	m_aFrame->CleanUpMatrix();
	// reset counters m_nEntriesPad, m_nHitsInPad, m_nChargeInPad
	m_aFrame->ResetCountersPad();
	// reset metadata if needed.  Not desired when processing multiframe
	if(rewind_metadata) m_aFrame->RewindMetaDataValues();

}

/**
 * Timepix3 specific
 */
int FramesHandler::readOneStampTimepix3(vector<unsigned int> stamp) {

	RewindAll();
	m_aFrame->IncreaseId(); // the first time it'll come to 0 (initialized at -1)
	m_nFrames++;

	// Identify the type of file
	Int_t width =  __timepix3_framewidth; // default, no meta-data available
	Int_t height = __timepix3_frameheight; // default, no meta-data available
	// Store this info in class members
	m_width = width;
	m_height = height;

	int ToA_MSB = 0;
	int ToA_LSB = 0;
	int ToA = 0;

	// For now fill 6 values
	// Take Id, [ x, y, ToT, ToA16, ToA14, fToA ] (6 values, from 1 to 6)
	unsigned int itr;
	unsigned int s_size = (unsigned int) stamp.size();
	for( itr = 0 ; itr < s_size ; itr += __timepix3_stampLength+1 ) { // Every 7 values (Throw the Id column)

		//cout << stamp[itr+0] << ", " << stamp[itr+1] << ", " << stamp[itr+2] << ", " << stamp[itr+3] << ", ";
		//cout << stamp[itr+4] << ", " << stamp[itr+5] << " : " << stamp[itr+6] << endl;

		// ToA is 30 bits
		// stamp[itr+4] : 16 bits (MSB) from the FPGA
		// stamp[itr+5] : 14 bits (LSB) from the chip
		ToA_MSB = (stamp[itr+4] << 14);
		ToA_LSB = stamp[itr+5];
		ToA = ToA_MSB | ToA_LSB;

		m_aFrame->FillOneElement ( 	stamp[itr+1],				// x
				stamp[itr+2],			// y
				__timepix3_framewidth,	// width
				stamp[itr+3], 			// ToT


				ToA,					// ToA ()
				stamp[itr+6]			// fastToA
		);

	}

	SetnX(width);
	SetnY(height);

	// processed one "frame".

	return 1;
}
/**
 * Dosepix specific
 */
int FramesHandler::readOneFrameDosepix(vector<int> frame, string /* timestamp */) {

	RewindAll();
	m_aFrame->IncreaseId(); // the first time it'll come to 0 (initialized at -1)
	m_nFrames++;

	// Identify the type of file
	Int_t width =  __dosepix_framewidth; // default, no meta-data available
	Int_t height = __dosepix_framewidth; // default, no meta-data available
	// Store this info in class members
	m_width = width;
	m_height = height;

	// First of all the string comes in reverse order the svc file
	// timestamp, 255, 254 .... 0
	reverse( frame.begin(), frame.end() );

	// Dosepix comes in this layout
	// 0 16 ........ 240
	// |
	// |
	// |
	// |
	// 15 .......... 255
	// --- Wirebonds ---

	// now read the map and fill the matrix
	vector<int>::iterator it = frame.begin();
	pair<int, int> pix;
	int x, y;
	int pos = 0;
	//cout << "----------------------------------------------" << endl;
	for ( ; it != frame.end(); it++ ) {

		pix = XtoXY(pos, __dosepix_framewidth); // this comes in the standard medipix layout

		y = pix.first;
		y = ( y - __dosepix_framewidth + 1 ) * -1;

		x = pix.second;

		//cout << pos << " | (" << pix.first << "," << pix.second << ") --> (" << x << "," << y << ")" << endl;

		// NO ZERO Suppression.  Only fill the element when the entry is different than 0.
		if ( (*it) > 0 ) {
			m_aFrame->FillOneElement ( x, y, __dosepix_framewidth, (*it) );
		}

		pos++;
	}
	//cout << " !!!!!!!!! frame " << pos << endl;
	SetnX(width);
	SetnY(height);

	// Metadata in the future ? // FIXME

	// processed one frame
	return 1;
}

pair<Int_t, Int_t> FramesHandler::XtoXY(Int_t X, Int_t dimX){
	return make_pair(X % dimX, X/dimX);
}

int FramesHandler::readOneFrame(TString fullFileName, TString fullDSCFileName, int * ftype){

	RewindAll();
	m_aFrame->IncreaseId(); // the first time it'll come to 0 (initialized at -1)
	m_nFrames++;

	// Identify the type of file, 256x256 or XYC
	Int_t width = 256; // default
	Int_t height = 256; // default

	int nFramesToRead = 1;
	//Int_t frameType = IdentifyTypeOfInput(fullDSCFileName, width, height, nFramesToRead);
	//if ( frameType == 0 ) {
	// try the newer layout
	int frameType = IdentifyTypeOfInput2(fullDSCFileName, width, height, nFramesToRead);
	//}
	*ftype = frameType;
	cout << "[INFO] the frame type is " << frameType << " and the number of frames associated is " << nFramesToRead << endl;

	if ( frameType == 0 ) { // if still 0
		std::cout << "[ERROR] could not determine the frame type --> " << fullDSCFileName << std::endl;
		return -1; // -1 frames to read
	}

	// Store this info in class members
	m_width = width;
	m_height = height;

	if(nFramesToRead > 1) { // multiframe
		return nFramesToRead;
	}

	/* ******************* */
	/* open the frame file */
	TString typeS = "";
	fstream filestr;

	Int_t temp = 0;
	Long_t totalCounts = 0;
	Long_t pixelHits = 0;
	Int_t cntri = 0;
	Int_t cntrj = 0;
	Int_t wholePadCntr = 0;

	if(frameType == (FSAVE_ASCII | FSAVE_I16) ||
			frameType == (FSAVE_ASCII | FSAVE_U32) ||
			frameType == (FSAVE_ASCII | FSAVE_DOUBLE) )
	{
		filestr.open(fullFileName, fstream::in);

		typeS = TYPE_256x256_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")" << std::endl;

		int maxocc = m_width*m_height;
		cout << "maxocc " << maxocc << " | " << filestr.good() << endl;
		while ( filestr.good() && wholePadCntr < maxocc ) {

			filestr >> temp;

			// See if the end of file has been reached
			if ( filestr.eof() ) break;

			if(temp > 0) { // only if this is a valid entry.  Since this is matrix even zeros are showing here
				totalCounts += temp;
				pixelHits++;
				m_aFrame->FillOneElement(cntri, cntrj, m_width, temp);
			}
			cntri++;

			if(cntri == m_width){
				cntrj++;
				cntri = 0;
			}
			wholePadCntr++;
		}
	}
	else if(frameType == (FSAVE_ASCII | FSAVE_I16 | FSAVE_SPARSEXY) ||
			frameType == (FSAVE_ASCII | FSAVE_U32 | FSAVE_SPARSEXY) ||
			frameType == (FSAVE_ASCII | FSAVE_DOUBLE | FSAVE_SPARSEXY) )
	{
		filestr.open(fullFileName, fstream::in);

		typeS = TYPE_XYC_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")" << std::endl;

		std::map<Int_t, Int_t> frameMap;
		for(Int_t itr = 0 ; itr < MAX_FRAME_OCC ; itr++)
		{
			frameMap[itr] = 0;
		}

		Int_t fillTime = 0;
		Int_t colr = 0, rowr = 0;
		while (filestr.good())
		{
			filestr >> temp;

			// See if the end of file has been reached
			if ( filestr.eof() ) break;

			if(fillTime%3 == 0) /* got X */ /* There is a harmles extra
					     read at this point.  No risk.*/
			{
				colr = temp;
				//std::cout << "colr: " << colr << std::endl;
				fillTime++;
			}
			else if(fillTime%3 == 1) /* got Y */
			{
				rowr = temp;
				//std::cout << "rowr: " << rowr << std::endl;
				fillTime++;
			}
			else if(fillTime%3 == 2) /* got Counts */
			{
				//std::cout << "val: " << temp << std::endl;

				totalCounts += temp;
				pixelHits++;
				m_aFrame->FillOneElement(colr, rowr, m_width, temp);

				fillTime = 0;
			}

			cntri++;

			if(cntri == m_width){
				cntrj++;
				cntri = 0;
			}

		}

		/* now read the dictionary and fill the matrix */
		// no need for this anymore
		/*
		std::map<Int_t, Int_t>::iterator it;
		cntri = 0;
		cntrj = 0;

		for ( it=frameMap.begin() ; it != frameMap.end(); it++ )
		{

			m_aFrame->FillOneElement(cntri, cntrj, m_width, (*it).second);
			cntri++;

			if(cntri == MAX_FRAME_COL){
				cntrj++;
				cntri = 0;
			}

		}
		 */

	}
	else if(frameType == (FSAVE_ASCII | FSAVE_I16 | FSAVE_SPARSEX) ||
			frameType == (FSAVE_ASCII | FSAVE_U32 | FSAVE_SPARSEX) ||
			frameType == (FSAVE_ASCII | FSAVE_DOUBLE | FSAVE_SPARSEX) )
	{
		filestr.open(fullFileName, fstream::in);

		typeS = TYPE_XC_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")"<< std::endl;

		std::map<Int_t, Int_t> frameMap;
		for(Int_t itr = 0 ; itr < MAX_FRAME_OCC ; itr++)
		{
			frameMap[itr] = 0;
		}

		Int_t fillTime = 0;
		Int_t colr = 0;//, rowr = 0;
		while (filestr.good())
		{
			filestr >> temp;

			// See if the end of file has been reached
			if ( filestr.eof() ) break;

			if(fillTime%2 == 0) /* got X */ /* There is a harmless extra
					     read at this point.  No risk.*/
			{
				colr = temp;
				//std::cout << "colr: " << colr << std::endl;
				fillTime++;
			}
			else if(fillTime%2 == 1) /* got Counts */
			{

				totalCounts += temp;
				pixelHits++;
				frameMap[colr] = temp;

				fillTime = 0;
			}

			cntri++;

			if(cntri == m_width){
				cntrj++;
				cntri = 0;
			}

		}

		// now read the map and fill the matrix
		std::map<Int_t, Int_t>::iterator it;
		for ( it=frameMap.begin() ; it != frameMap.end(); it++ ) {
			m_aFrame->FillOneElement((*it).first, (*it).second);
		}

	}
	//////////////////////////////////////////////////////////////////////////
	// Binaries
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	// Binary matrix missing TODO
	//////////////////////////////////////////////////////////////////////////

	else if (
			( frameType == (FSAVE_BINARY | FSAVE_I16 | FSAVE_SPARSEXY) )
	) {

		typeS = TYPE_XYC_BIN_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")"<< std::endl;

		//filestr.open(fullFileName, fstream::binary);
		filestr.open(fullFileName, istream::in);

		unsigned int bytesRead = 0;
		char tempByte[4];
		unsigned int x = 0x0, y = 0x0, counts = 0x0;

		while (filestr.good())
		{

			// read X, 32 bits
			x = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);

			// See if the end of file has been reached
			if ( filestr.eof() ) break;

			push_back_nbytes(&x, tempByte, 4);

			// read Y, 32 bits
			y = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);
			push_back_nbytes(&y, tempByte, 4);

			// Read Counts, 16 bits
			counts = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			push_back_nbytes(&counts, tempByte, 2);

			//printf("%x, %x --> %x\n", x, y, counts);
			//printf("%d, %d --> %d\n", x, y, counts);

			bytesRead += 10;

			m_aFrame->FillOneElement(x, y, m_width, counts);

		}

	}	else if (
			( frameType == (FSAVE_BINARY | FSAVE_U32 | FSAVE_SPARSEXY) )
	) {

		typeS = TYPE_XYC_BIN_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")"<< std::endl;

		//filestr.open(fullFileName, fstream::binary);
		filestr.open(fullFileName, istream::in);

		unsigned int bytesRead = 0;
		char tempByte[4];
		unsigned int x = 0x0, y = 0x0, counts = 0x0;

		while (filestr.good())
		{

			// read X, 32 bits
			x = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);

			// See if the end of file has been reached
			if ( filestr.eof() ) break;

			push_back_nbytes(&x, tempByte, 4);

			// read Y, 32 bits
			y = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);
			push_back_nbytes(&y, tempByte, 4);

			// Read Counts, 32 bits
			counts = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);
			push_back_nbytes(&counts, tempByte, 4);

			//printf("%x, %x --> %x\n", x, y, counts);
			//printf("%d, %d --> %d\n", x, y, counts);

			bytesRead += 12;

			m_aFrame->FillOneElement(x, y, m_width, counts);

		}

	} 	else if (
			( frameType == (FSAVE_BINARY | FSAVE_I16 | FSAVE_SPARSEX ) )
	) {

		typeS = TYPE_XYC_BIN_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")"<< std::endl;

		//filestr.open(fullFileName, fstream::binary);
		filestr.open(fullFileName, istream::in);

		unsigned int bytesRead = 0;
		char tempByte[4];
		unsigned int x = 0x0, counts = 0x0;

		while (filestr.good())
		{

			// read X, 32 bits
			x = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);

			// See if the end of file has been reached
			if ( filestr.eof() ) break;

			push_back_nbytes(&x, tempByte, 4);

			// Read Counts, 16 bits
			counts = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			push_back_nbytes(&counts, tempByte, 2);

			//printf("%x, %x --> %x\n", x, y, counts);
			//printf("%d, %d --> %d\n", x, y, counts);

			bytesRead += 6;

			m_aFrame->FillOneElement(x, counts);

		}

	} else if (
			( frameType == (FSAVE_BINARY | FSAVE_U32 | FSAVE_SPARSEX ) )
	) {

		typeS = TYPE_XYC_BIN_STRING;
		std::cout << "[INFO] opening: " << fullFileName << "  --|and|-->  " << fullDSCFileName
				<< " --> " << typeS << "(" << frameType << ")"<< std::endl;

		//filestr.open(fullFileName, fstream::binary);
		filestr.open(fullFileName, istream::in);

		unsigned int bytesRead = 0;
		char tempByte[4];
		unsigned int x = 0x0, counts = 0x0;

		while (filestr.good())
		{

			// read X, 32 bits
			x = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);

			// See if the end of file has been reached
			if ( filestr.eof() ) break;

			push_back_nbytes(&x, tempByte, 4);

			// Read Counts, 32 bits
			counts = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);
			push_back_nbytes(&counts, tempByte, 4);

			//printf("%x, %x --> %x\n", x, y, counts);
			//printf("%d, %d --> %d\n", x, y, counts);

			bytesRead += 8;

			m_aFrame->FillOneElement(x, counts);

		}

	} else {
		/* Unknown case ... giving up here */
	}

	// FIXME ! // get this from the dsc file
	SetnX(width);
	SetnY(height);

	filestr.close();

	// Now fetching MetaData
	fstream META_filestr;
	Int_t nLineMetaFile = 0;
	META_filestr.exceptions ( ifstream::eofbit | ifstream::failbit | ifstream::badbit );

	// read metadata
	try {

		META_filestr.open(fullDSCFileName, fstream::in);

		Char_t * METAtemp;
		METAtemp = new char[META_DATA_LINE_SIZE];

		while (META_filestr.good()) {

			META_filestr.getline(METAtemp, META_DATA_LINE_SIZE);
			if(m_ParseAndHold) parseMetaLine((TString)METAtemp, nLineMetaFile); // pass currentLine

			if(nLineMetaFile++ == m_metaBit){

				m_aFrame->FillMetaData((TString)METAtemp, m_metaCode);
				m_ParseAndHold = true;

			}

		}
	}
	catch (ifstream::failure e){
		//std::cout << "[ERROR] Exception opening/reading file: " << fullDSCFileName << std::endl;
		//std::cout << "        probably file doesn't not exist.  giving up." << std::endl;
		//exit(1);
	}
	META_filestr.close();

	return 1; // single frame
}

void FramesHandler::push_back_nbytes(unsigned int * val, char * bytes, int nbytes) {

	// indexes go like this
	// 0 --> lower byte  0x......XX
	// 3 --> higher byte 0xXX......

	*val &= 0x00000000;
	unsigned int tempVal;

	for(int idx = nbytes-1 ; idx >= 0 ; idx--) {

		// Get the byte
		tempVal = 0x0;
		tempVal ^= bytes[idx];
		// Clean up to have info for only one byte
		tempVal &= 0x000000ff;
		// Switch it to the right place
		for(int sw = 0 ; sw < idx ; sw++){
			tempVal = tempVal << 8;
		}
		// XOR the value
		*val ^= tempVal;

	}

}

void FramesHandler::push_back_nbytes(long long * val, char * bytes, int nbytes) {

	// indexes go like this
	// 0 --> lower byte  0x......XX
	// 3 --> higher byte 0xXX......

	*val &= 0x0000000000000000; // 64 bits
	unsigned int tempVal;

	for(int idx = nbytes-1 ; idx >= 0 ; idx--) {

		// Get the byte
		tempVal = 0x0;
		tempVal ^= bytes[idx];
		// Clean up to have info for only one byte
		tempVal &= 0x000000ff;
		// Switch it to the right place
		for(int sw = 0 ; sw < idx ; sw++){
			tempVal = tempVal << 8;
		}
		// XOR the value
		*val ^= tempVal;

	}

}

/* load a single frame pixel (X,Y,C) + truth info */
Bool_t FramesHandler::LoadFramePixel(Int_t col, Int_t row, Int_t counts, Double_t truthE, Double_t E){

	m_aFrame->FillOneElement(col, row, m_width, counts, truthE, E);

	return true;
}

/* load a single frame pixel (X,Y,C) */
Bool_t FramesHandler::LoadFramePixel(Int_t col, Int_t row, Int_t counts){

	m_aFrame->FillOneElement(col, row, m_width, counts);

	return true;
}


void FramesHandler::SetLVL1(Int_t col, Int_t row, Int_t lvl1){
	m_aFrame->SetLVL1(col, row, m_width, lvl1);
}

void FramesHandler::IncreaseCurrentFrameId(){
	m_aFrame->IncreaseId();
}

/* load a single frame MetaData entry (METAString, Int_t metaCode)*/
Bool_t FramesHandler::LoadFrameMetaData(TString metastring , Int_t metacode){

	m_aFrame->FillMetaData(metastring, metacode);

	return true;
}

void FramesHandler::parseMetaLine(TString aMetaLine, Int_t currentLine){

	// It is probably a good idea to change this to
	// "switch case" because it may be a bit faster
	if(aMetaLine.Contains(ACQ_MODE_STRING, TString::kExact)){
		//std::cout << "-> 1" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = ACQ_MODE_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(ACQ_TIME_STRING, TString::kExact)){
		//std::cout << "-> 2" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = ACQ_TIME_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(CHIPBOARDID_STRING, TString::kExact)){
		//std::cout << "-> 3" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = CHIPBOARDID_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(DACS_STRING, TString::kExact)){
		//std::cout << "-> 4" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = DACS_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(HW_TIMER_STRING, TString::kExact)){
		//std::cout << "-> 5" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = HW_TIMER_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(INTERFACE_STRING, TString::kExact)){
		//std::cout << "-> 6" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = INTERFACE_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(POLARITY_STRING, TString::kExact)){
		//std::cout << "-> 7" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = POLARITY_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(START_TIME_STRING, TString::kExact)){
		//std::cout << "-> 8" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = START_TIME_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(START_TIME_S_STRING, TString::kExact)){
		//std::cout << "-> 9" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = START_TIME_S_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(MPX_CLOCK_STRING, TString::kExact)){
		//std::cout << "-> 9" << std::endl;
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = MPX_CLOCK_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(HV_STRING, TString::kExact)){
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = HV_CODE;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(TIMEPIX_CLOCK_STRING, TString::kExact)
			&&
			! aMetaLine.Contains(TIMEPIX_CLOCKINMHZ_STRING, TString::kExact)
	){
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = TIMEPIX_CLOCK_CODE ;
		m_ParseAndHold = false;
	}
	else if(aMetaLine.Contains(TIMEPIX_CLOCKINMHZ_STRING, TString::kExact)){
		m_metaBit = currentLine + ACTUAL_VALUE_OFFSET;
		m_metaCode = TIMEPIX_CLOCKINMHZ_CODE ;
		m_ParseAndHold = false;
	}
	else{
		//metaBit = -1;
	}

}
/*
TH2I * FramesHandler::getAFrameHist(TString fullFileName, TString fullDSCFileName, TString){

	getAFrameMatrix_flag = false; getAFrameHist_flag = true;
	Bool_t frameTrigger = readOneFrame(fullFileName, fullDSCFileName);

	// check the integrity of the frame
	if(!frameSupervisor())
	{
		std::cout << "[ERROR] frame " << fullFileName << "\n";
		std::cout << "        does not seem to be a well formatted frame\n";
		std::cout << "        giving up." << std::endl;
		exit(1);
	}

	TH2I * frameHist = 0;
	TH2I * h_emptyFrame = 0;

	Int_t * matrix_i;
	matrix_i = m_aFrame->GetPixelsMatrix();
	if (matrix_i != 0)
	{
		frameHist = getHistFrame(m_aFrame->GetFrameId(), matrix_i);
	}
	else
	{
		std::cout << "[ERROR] frame matrix not initialized ... givin up" << std::endl;
		exit(1);
	}


	if(frameTrigger)
		return frameHist;
	else
		return h_emptyFrame;
}
 */
Int_t ** FramesHandler::getAFrameMatrix(TString fullFileName, TString /*fullDSCFileName*/){

	getAFrameMatrix_flag = true; getAFrameHist_flag = false;
	//Bool_t frameTrigger = readOneFrame(fullFileName, fullDSCFileName);

	// check the integrity of the frame
	if(!frameSupervisor())
	{
		std::cout << "[ERROR] frame " << fullFileName << "\n";
		std::cout << "        does not seem to be a well formatted frame\n";
		std::cout << "        giving up." << std::endl;
		exit(1);
	}
	/*
    if (m_aFrame->GetPixelsMatrix() != 0)
    return m_aFrame->GetPixelsMatrix();
	 */

	// not used !
	return 0x0;
}

Int_t FramesHandler::XYtoX(int x, int y, Int_t dimX){
	return y * dimX + x;
}

Bool_t FramesHandler::ProcessMultiframe(TString datafile, TString dscfile, TString idxfile, WriteToNtuple * wte, int ftype) {

	// I already know the type at this point

	fstream filestr;
	filestr.open(datafile, fstream::in);

	fstream filestr_idx;
	filestr_idx.open(idxfile, fstream::in);

	cout << "[INFO] opening    : " << datafile << "  --|and|-->  " << dscfile << " --> " << "(" << ftype << ")"<< endl;
	cout << "       multiframe : " << idxfile << endl;

	///////////////////////////////////////////////////////////////////////////////////////
	// Meta data has to be processed at the same time
	// Now fetching MetaData
	fstream META_filestr;
	Int_t nLineMetaFile = 0;
	//META_filestr.exceptions ( ifstream::eofbit | ifstream::failbit | ifstream::badbit );

	// Prepare to read metadata
	META_filestr.open(dscfile.Data(), fstream::in);
	Char_t * METAtemp;
	METAtemp = new char[META_DATA_LINE_SIZE];

	while ( META_filestr.good() ) {

		META_filestr.getline(METAtemp, META_DATA_LINE_SIZE);
		if(m_ParseAndHold) parseMetaLine((TString)METAtemp, nLineMetaFile); // pass currentLine

		if(nLineMetaFile++ == m_metaBit){

			//cout << METAtemp << endl;
			m_aFrame->FillMetaData((TString)METAtemp, m_metaCode);
			m_ParseAndHold = true;

		}

	}
	META_filestr.close();
	///////////////////////////////////////////////////////////////////////////////////////


	if(
			ftype == (FSAVE_ASCII | FSAVE_I16 | FSAVE_SPARSEX) ||
			ftype == (FSAVE_ASCII | FSAVE_U32 | FSAVE_SPARSEX)
	) {

		std::map<int, int> frameMap;
		std::map<int, int>::iterator it;

		int fillTime = 0;
		int xi = 0;
		char temps[256];
		bool newframe = false;
		int cntr = 0;

		while (filestr.good()) {

			// read data
			filestr.getline(temps, 256);
			// check EOF
			if(filestr.eof()) break;     // EOF

			// check next frame indicator
			if(temps[0] == '#') { // new frame
				newframe = true;
			}

			if(!newframe) {
				// extract both numbers
				string tempval;
				for(int i = 0 ; i < 256 ; i++) {
					//printf("%x\n",temps[i]);
					if(temps[i] != 0x9 && (temps[i] != 0xa && temps[i] != 0xd) ) { // 0x9 horizontal tab
						tempval.append(1, temps[i]);
					} else { // finish

						if(fillTime%2 == 0) {        // X
							xi = atoi( tempval.c_str() );
							tempval.clear();
							fillTime++;
						} else if(fillTime%2 == 1) { // C
							frameMap[xi] = atoi( tempval.c_str() );
							tempval.clear();
							fillTime = 0;
						}

					}
					if(temps[i] == 0xa || temps[i] == 0xd) break; // end of this line
				}

				//cout << xi << ", " << frameMap[xi] << endl;

			}

			if(newframe) {

				//cout << "Filling frame " << cntr << endl;

				// now read the dictionary and fill the matrix
				for ( it=frameMap.begin() ; it != frameMap.end(); it++ ) {
					// X,C
					//cout << (*it).first << ", " << (*it).second << endl;
					m_aFrame->FillOneElement((*it).first, (*it).second);
				}
				//
				SetnX(m_width);
				SetnY(m_height);
				m_aFrame->SetId(cntr);

				// clean map
				frameMap.clear();
				// finally fill vars
				wte->fillVars(this, false); // second parameters tells not to rewind metadata
				// ready to read next
				newframe = false;
				cntr++;
			}

		}

		filestr.close();

	} else if(
			ftype == (FSAVE_ASCII | FSAVE_I16 | FSAVE_SPARSEXY) ||
			ftype == (FSAVE_ASCII | FSAVE_U32 | FSAVE_SPARSEXY)
	) {

		std::map<int, int> frameMap;
		std::map<int, int>::iterator it;

		int fillTime = 0;
		int xi = 0, yi = 0;
		char temps[256];
		bool newframe = false;
		int cntr = 0;

		while (filestr.good()) {

			// read data
			filestr.getline(temps, 256);
			// check EOF
			if(filestr.eof()) break;     // EOF

			// check next frame indicator
			if(temps[0] == '#') { // new frame
				newframe = true;
			}

			if(!newframe) {
				// extract both numbers
				string tempval;
				for(int i = 0 ; i < 256 ; i++) {
					//printf("%x\n",temps[i]);
					if(temps[i] != 0x9 && (temps[i] != 0xa && temps[i] != 0xd) ) { // 0x9 horizontal tab
						tempval.append(1, temps[i]);
					} else { // finish

						if(fillTime%3 == 0) {        // X
							xi = atoi( tempval.c_str() );
							tempval.clear();
							fillTime++;
						} else if(fillTime%3 == 1){  // Y
							yi = atoi( tempval.c_str() );
							tempval.clear();
							fillTime++;
						} else if(fillTime%3 == 2) { // C

							frameMap[ XYtoX(xi, yi, m_width) ] = atoi( tempval.c_str() );
							tempval.clear();
							fillTime = 0;
						}


					}
					if(temps[i] == 0xa || temps[i] == 0xd) break; // end of this line
				}

				//cout << xi << ", " << yi << " : " << XYtoX(xi, yi, m_width) << endl;

			}

			if(newframe) {

				//cout << "Filling frame " << cntr << endl;

				// now read the dictionary and fill the matrix
				for ( it=frameMap.begin() ; it != frameMap.end(); it++ ) {
					// X,C
					//cout << (*it).first << ", " << (*it).second << endl;
					m_aFrame->FillOneElement((*it).first, (*it).second);
				}
				//
				SetnX(m_width);
				SetnY(m_height);
				m_aFrame->SetId(cntr);

				// clean map
				frameMap.clear();
				// finally fill vars
				wte->fillVars(this, false); // don't rewind metadata
				// ready to read next
				newframe = false;
				cntr++;
			}

		}
		filestr.close();

	} else if(
			ftype == (FSAVE_BINARY | FSAVE_I16 | FSAVE_SPARSEXY) ||
			ftype == (FSAVE_BINARY | FSAVE_U32 | FSAVE_SPARSEXY)
	) {

		// 2 := I16, 4 := U32
		int nCountBytes = 2;
		if ( ftype == (FSAVE_BINARY | FSAVE_I16 | FSAVE_SPARSEXY) ) nCountBytes = 2;
		if ( ftype == (FSAVE_BINARY | FSAVE_U32 | FSAVE_SPARSEXY) ) nCountBytes = 4;

		// Information from the idx file (indexation)
		// starting all at 0
		long long dscPos = 0, dataPos = 0, sfPos = 0;
		long long pixelAppendCntr = 0;

		// for data
		long long bytesRead = 0;
		char tempByte[4];
		unsigned int x = 0x0, y = 0x0, counts = 0x0;

		// Get Idx values.  This information indexes from second frame
		GetIdxValues(&filestr_idx, &dscPos, &dataPos, &sfPos);

		// Clean up
		//getFrameStructObject()->CleanUpMatrix();
		//getFrameStructObject()->ResetCountersPad();

		while ( filestr.good() ) {

			// read X, 32 bits
			x = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);

			// check EOF
			if(filestr.eof()) break;     // EOF

			push_back_nbytes(&x, tempByte, 4);

			// read Y, 32 bits
			y = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			filestr.get(tempByte[2]);
			filestr.get(tempByte[3]);
			push_back_nbytes(&y, tempByte, 4);

			// Read Counts, 16 bits or 32 bits
			counts = 0x0;
			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			if(nCountBytes == 4){ // U32
				filestr.get(tempByte[2]);
				filestr.get(tempByte[3]);
			}
			push_back_nbytes(&counts, tempByte, nCountBytes);

			//printf("%x, %x --> %x\n", x, y, counts);
			//printf("0x%x 0x%x, --> 0x%x | %d %d, --> %d\n", x, y, counts, x, y, counts);
			if(nCountBytes == 2) bytesRead += 10;
			if(nCountBytes == 4) bytesRead += 12;

			// Store this pixel info
			if (bytesRead < dataPos) {

				getFrameStructObject()->FillOneElement(x, y, 256, counts);
				pixelAppendCntr++;

			} else { // Change Frame

				getFrameStructObject()->IncreaseId();
				getFrameStructObject()->SetnX(256);
				getFrameStructObject()->SetnY(256);
				wte->fillVars(this, false); // don't rewind medatada
				// Get new Idx values.
				GetIdxValues(&filestr_idx, &dscPos, &dataPos, &sfPos);
				// Don't rewind the bytes counter.  dataPos matches this value
				// rewind pixel counter per frame
				pixelAppendCntr = 0;

			}

		}

	}

	//cout << "[INFO] Multiframe ... " << cntr << " frames processed" << endl;

	return true;
}

void FramesHandler::GetIdxValues(fstream * filestr_idx, long long * dscPos, long long * dataPos, long long * sfPos) {

	// for idx
	char tempByte64[8]; // 64 bits

	if(filestr_idx->good()) {
		// dsc pos
		for(int i = 0 ; i < 8 ; i++) {
			filestr_idx->get(tempByte64[i]);
		}
		push_back_nbytes(dscPos, tempByte64, 8);
		// data pos
		for(int i = 0 ; i < 8 ; i++) {
			filestr_idx->get(tempByte64[i]);
		}
		push_back_nbytes(dataPos, tempByte64, 8);
		// sf pos
		for(int i = 0 ; i < 8 ; i++) {
			filestr_idx->get(tempByte64[i]);
		}
		push_back_nbytes(sfPos, tempByte64, 8);

		//cout << "dscPos = " << *dscPos << ", dataPos = " << *dataPos << ", sfPos = " << *sfPos << endl;

	}

}

/*
void FramesHandler::push_back_nbytes(unsigned int * val, char * bytes, int nbytes) {

	// indexes go like this
	// 0 --> lower byte  0x......XX
	// 3 --> higher byte 0xXX......

 *val &= 0x00000000;
	unsigned int tempVal;

	for(int idx = nbytes-1 ; idx >= 0 ; idx--) {

		// Get the byte
		tempVal = 0x0;
		tempVal ^= bytes[idx];
		// Clean up to have info for only one byte
		tempVal &= 0x000000ff;
		// Switch it to the right place
		for(int sw = 0 ; sw < idx ; sw++){
			tempVal = tempVal << 8;
		}
		// XOR the value
 *val ^= tempVal;

	}

}
 */

Bool_t FramesHandler::frameSupervisor(){

	/** check to perform when the matrix for one frame
	 *  has been completely loaded.
	 */

	if(m_aFrame->GetEntriesPad() != MAX_FRAME_OCC)
	{
		return false;
	}
	else
	{
		return true;
	}

}
