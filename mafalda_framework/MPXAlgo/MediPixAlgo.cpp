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

#ifndef MediPixAlgo_cxx
#define MediPixAlgo_cxx

#include "MediPixAlgo.h"
#include "MPXStoreGate/CandidateContainer.h"
#include "MAFTools/MAFTools.h"
#include "MPXAlgo/Highlighter.h"
#include "BlobsFinder/BlobsFinder.h"

#include "TLatex.h"

ClassImp(MediPixAlgo)

MediPixAlgo::MediPixAlgo(){

	m_globalCounter = 0;
	mySignals.nextFrame = 0;
	m_hasConfiguration = false;
	m_readConfiguration = false;

	// Set the Log service
	Log.setAlgoName("MediPixAlgo");
	Log.OutputLevel = MSG::INFO;

	// Extended mask
	m_inputExtendedMask = false;
	m_extendedMaskFile = "";

}

/**
 * Read configuration from filename
 *
 */
void MediPixAlgo::ReadConfiguration(const Char_t * fileN){

	m_configurationFilename = fileN;
	m_readConfiguration = true;
}
/**
 * Read configuration from default filename
 *
 */
void MediPixAlgo::ReadConfiguration(){
	m_readConfiguration = true;
}

void MediPixAlgo::PreInitMPXAlgo(){

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Deal with the extended mask here

	if (m_inputExtendedMask) { // In this case the code will produce the extended mask

		// If not filename given
		if( m_extendedMaskFile == TString("") ) {
			TString oname = TString(GetOutputBaseName().c_str());
			oname += "_ExtendedMask.txt";
			m_extendedMaskFile = oname;
		}

		Log << MSG::INFO << "Attempting to load extended mask file : " << m_extendedMaskFile << endreq;
		infile_extendedmask = new ifstream(m_extendedMaskFile.Data());
		if( ! infile_extendedmask->good() ) {
			Log << MSG::ERROR << "Coudn't load extended mask file : " << m_extendedMaskFile << ", giving up." << endreq;
			exit(1);
		}
		LoadExtendedMask();
		DumpExtendedMask();
		Log << MSG::INFO << "Extended mask successfully loaded. " << GetAlgoName() << " will use it." << endreq;

	}

}

void MediPixAlgo::InitMPXAlgo(TString algoName_i, MediPixAnalysisCore * myCore){

	m_MPXdata = myCore;
	algoName = algoName_i;
	Log.setAlgoName(algoName);

	// some internal config
	m_MarkerFontSizeDrawingSeparateWindow = 0.05;
}

/**
 * Some task pending after calling daughter Init member but still in Initialization stage
 */
void MediPixAlgo::PostInitMPXAlgo(){


	////////////////////////////////////////////////////////////////////////////////////////////////////
	// Deal with the configuration

	if(!m_hasConfiguration || !m_readConfiguration)
		return;

	if(m_configurationFilename.empty())
		m_configurationFilename = "configuration_algos.txt";

	Log << MSG::ALWAYS << "Reading configuration from \""
			<< m_configurationFilename.c_str() << "\" config file"<< endreq;

	ifstream ifs(m_configurationFilename.c_str(), ifstream::in);
	if(!ifs.is_open()){
		Log << MSG::ALWAYS << "Could not find configuration file.  Using default values." << endreq;
		return;
	}

	// see how many configuration objects the algos have registered
	Int_t nObjs = GetNumberOfSpecialObjectsLessEqualThan(MPXDefs::CONF);
	Log << MSG::ALWAYS << "There are " << nObjs << " configuration objects registered in the" << endreq;
	Log << MSG::ALWAYS << "  StoreGate up to this point.  Working on these values." << endreq;

	if(nObjs == 0){
		Log << MSG::ALWAYS << "There is a configuration file and there seem to be values" << endreq;
		Log << MSG::ALWAYS << "  registered but I can not find them in the StoreGate." << endreq;
		Log << MSG::ALWAYS << "  This should not be hapenning, something's wrong !" << endreq;
		return;
	}

	Char_t temp[MAX_OCCUPANCY_CONF_LINE];
	string tempS;
	Int_t npos;
	string t_type, t_value, t_algo, t_varname;
	while(ifs.good())
	{
		ifs.getline(temp, MAX_OCCUPANCY_CONF_LINE);
		tempS.clear();
		tempS = temp;

		// type
		npos = tempS.find_last_of(" ");
		t_type = tempS.substr(npos+1, tempS.length()-1);
		tempS = tempS.substr(0, npos);

		// value
		npos = tempS.find_last_of(" ");
		t_value = tempS.substr(npos+1, tempS.length()-1);
		tempS = tempS.substr(0, npos);

		// var
		npos = tempS.find_last_of(" ");
		t_varname = tempS.substr(npos+1, tempS.length()-1);
		tempS = tempS.substr(0, npos);

		// algo, the only part left
		t_algo = tempS;

		// now go touch the configuration
		TouchConfiguraton(t_type, t_value, t_varname, t_algo);
	}


}

void MediPixAlgo::LoadExtendedMask(){

	int tempv;
	int cntr = 0;
	while( infile_extendedmask->good() ) {
		*infile_extendedmask >> tempv;
		if(infile_extendedmask->eof()) break; // check EOF
		m_extendedMask.push_back( tempv );
		cntr++;
	}

	if( cntr == 0 ) {
		Log << MSG::ERROR << "The mask file exist and was read but no information was found.  Giving up ..." << endreq;
		exit(1);
	}

	Log << Log.OutputLevel << "Number of pixels in the mask = " << cntr << endreq;

}

void MediPixAlgo::DumpExtendedMask(){

	int cntr = 0;
	int maxocc = (int) m_extendedMask.size();
	Log << MSG::DEBUG << "-- Dump extended mask ----------------------------------------------------------------" << endreq;
	for(int i = 0 ; i < maxocc ; i++) {

		if( cntr%256 == 0) { Log << MSG::DEBUG << endreq; }
		else { Log << MSG::DEBUG << ' '; }

		Log << MSG::DEBUG << m_extendedMask[i];

		cntr++;
	}
	Log << MSG::DEBUG << endreq;
	Log << MSG::DEBUG << "-- End of extended mask ----------------------------------------------------------------" << endreq;

}

bool MediPixAlgo::PixelInExtendedMask (pair<int, int> pix, int sizex) {

	if ( std::find(m_extendedMask.begin(), m_extendedMask.end(), MAFTools::XYtoX( pix, sizex )) != m_extendedMask.end() ) {

		return true;
	}

	return false;
}

bool MediPixAlgo::PixelInExtendedMask (int x, int y, int sizex) {

	if ( std::find(m_extendedMask.begin(), m_extendedMask.end(), MAFTools::XYtoX( x, y, sizex )) != m_extendedMask.end() ) {
		return true;
	}

	return false;
}


bool MediPixAlgo::PixelInExtendedMaskBelongsToCluster (blob cl, int width, int height) {
	return MAFTools::PixelInListIsAClusterConstituent(cl, m_extendedMask, width, height);
}

void MediPixAlgo::TouchConfiguraton(string type, string val, string varname, string algo)
{

	MPXDefs::SpecialObjs typeFetch;
	if(type == "int")
		typeFetch = MPXDefs::CONF_INT;
	else if(type == "float")
		typeFetch = MPXDefs::CONF_FLOAT;
	else if(type == "double")
		typeFetch = MPXDefs::CONF_DOUBLE;
	else if(type == "bool")
		typeFetch = MPXDefs::CONF_BOOL;

	std::vector<CandidateContainer *> conf = GetObjectsSpecialWithAuthor(typeFetch, GetAlgoName());
	std::vector<CandidateContainer *>::iterator itr = conf.begin();

	if(type == "int"){
		for( ; itr != conf.end() ; itr++){
			ConfigurationValue<Int_t> * p = (ConfigurationValue<Int_t> *)(*itr);
			if(p->GetConfigName() ==  varname && algoName == algo){
				Log << MSG::ALWAYS << "  found " << p->GetConfigName() << " | setting value to " << val << endreq;
				p->SetConfigValue( atoi(val.c_str()) );
			}
		}
	}
	else if(type == "float"){
		for( ; itr != conf.end() ; itr++){
			ConfigurationValue<Float_t> * p = (ConfigurationValue<Float_t> *)(*itr);
			if(p->GetConfigName() ==  varname && algoName == algo){
				Log << MSG::ALWAYS << "  found " << p->GetConfigName() << " | setting value to " << val << endreq;
				p->SetConfigValue( atof(val.c_str()) );
			}
		}
	}
	else if(type == "double"){
		for( ; itr != conf.end() ; itr++){
			ConfigurationValue<Double_t> * p = (ConfigurationValue<Double_t> *)(*itr);
			if(p->GetConfigName() ==  varname && algoName == algo){
				Log << MSG::ALWAYS << "  found " << p->GetConfigName() << " | setting value to " << val << endreq;
				p->SetConfigValue( atof(val.c_str()) );
			}
		}
	}
	else if(type == "bool"){
		for( ; itr != conf.end() ; itr++){
			ConfigurationValue<Bool_t> * p = (ConfigurationValue<Bool_t> *)(*itr);
			if(p->GetConfigName() ==  varname && algoName == algo){
				Log << MSG::ALWAYS << "  found " << p->GetConfigName() << " | setting value to " << val << endreq;
				p->SetConfigValue( atof(val.c_str()) );
			}
		}
	}

}

void MediPixAlgo::SetConfigurationFilename(const Char_t * cnf){

	m_configurationFilename = cnf;

}

TTree * MediPixAlgo::getMyTree(){

	return m_MPXdata->outputNtuple->getTree(algoName);
}

TFile * MediPixAlgo::getMyROOTFile(){
	return m_MPXdata->outputNtuple->getFile();
}

void MediPixAlgo::ListenToTheManager(AnalysisManager * myManager_i)
{
	// get a handle in the manager ptr
	m_myManager = myManager_i;

	// information from manager
	m_inputFilesMap = m_myManager->GetInputFilesMap();
}


/** ***********************************************
 * 
 * DM Interface ! 
 *
 */

/**
 * DisconnectMeFromRun informs the manager that this algorithm should
 * be disconnected starting from the next event.
 */
void MediPixAlgo::DisconnectMeFromRun () {
	m_myManager->DisconnectAlgo(algoName);
}

std::map<int,int> MediPixAlgo::GetTOTMap(){
	return m_myManager->GetTOTMap();
}

Int_t MediPixAlgo::GetMatrixElement(Int_t col, Int_t row){ 
	return m_myManager->GetMatrixElement(col, row);
}

Int_t MediPixAlgo::GetToA(Int_t col, Int_t row){
	return m_myManager->GetToA(col, row);
}
Int_t MediPixAlgo::GetFastToA(Int_t col, Int_t row){
	return m_myManager->GetFastToA(col, row);
}

Int_t MediPixAlgo::GetMatrixElement(std::pair<Int_t, Int_t> coordinatesPair){
	return m_myManager->GetMatrixElement(coordinatesPair.first, coordinatesPair.second);
}

Int_t MediPixAlgo::GetToA(std::pair<Int_t, Int_t> coordinatesPair){
	return m_myManager->GetToA(coordinatesPair.first, coordinatesPair.second);
}

Int_t MediPixAlgo::GetFastToA(std::pair<Int_t, Int_t> coordinatesPair){
	return m_myManager->GetFastToA(coordinatesPair.first, coordinatesPair.second);
}

Int_t MediPixAlgo::GetCalibEnergy (Int_t col, Int_t row) {
	return m_myManager->GetCalibEnergy(MAFTools::XYtoC(col, row, GetWidth()));
}

// This method allows the user to register the calibrated energy in the data model
// That info doesn't really belong to the data as it doesn't come normally with it
// but in the next algo the user will have access to it through MediPixAlgo::GetCalibEnergy
void MediPixAlgo::SetCalibEnergy(std::pair<Int_t, Int_t> pix, Double_t e){
	SetCalibEnergy(pix.first, pix.second, e);
}
void MediPixAlgo::SetCalibEnergy(Int_t col, Int_t row, Double_t e){
	m_myManager->SetCalibEnergy(MAFTools::XYtoC(col, row, GetWidth()), e);
}

Int_t MediPixAlgo::GetLVL1(Int_t col, Int_t row){
	return m_myManager->GetToA(col, row);
	//return m_myManager->GetLVL1(col, row);
}

Int_t MediPixAlgo::GetLVL1(std::pair<Int_t, Int_t> coordinatesPair){
	return m_myManager->GetLVL1(coordinatesPair.first, coordinatesPair.second);
}

Double_t MediPixAlgo::GetMatrixElementMCEdep(Int_t col, Int_t row){
	return m_myManager->GetMatrixElementMCEdep(col, row);
}

Double_t MediPixAlgo::GetMatrixElementMCEdep(std::pair<Int_t, Int_t> coordinatesPair){
	return m_myManager->GetMatrixElementMCEdep(coordinatesPair.first, coordinatesPair.second);
}

Int_t MediPixAlgo::GetMatrixXdim(){
	return m_myManager->GetfWidth();
}
Int_t MediPixAlgo::GetMatrixYdim(){
	return m_myManager->GetfHeight();
}
/* ****************** */

UInt_t MediPixAlgo::GetMatrixCropAsBitWord(Int_t col1, Int_t col2, Int_t row1, Int_t row2){ 

	UInt_t cropMatrixBitWord = 0x0;
	UInt_t mask = 0x01;
	Int_t rowItr, colItr;

	for(rowItr = row1 ; rowItr < row2+1 ; rowItr++) /* exclude first, last rows and colums */
	{
		for(colItr = col1 ; colItr < col2+1 ; colItr++)
		{

			cropMatrixBitWord |= mask * ( (Bool_t) GetMatrixElement(colItr, rowItr));

			//std::cout << "--> In pixel: "<< GetMatrixElement(colItr, rowItr) << std::endl;
			//std::cout << "    * mask    : " << mask << std::endl;
			//std::cout << "    row, col: " << rowItr << ", " << colItr << std::endl;
			//std::cout << "    * bitword : " << cropMatrixBitWord << std::endl;

			mask = mask << 1;
		}
	}

	return cropMatrixBitWord;
}

std::vector<Int_t> MediPixAlgo::GetRowAsVector(Int_t row){

	std::vector<Int_t> rowV;
	Int_t xDim = GetMatrixXdim();

	for(Int_t colItr = 0; colItr < xDim ; colItr++)
	{
		rowV.push_back(GetMatrixElement(colItr, row));
	}

	return rowV;
}

void MediPixAlgo::DumpDACs(MSG::Level l) {

	vector<int> d = GetDAQs();
	vector<int>::iterator i = d.begin();

	Log << l << "DACs : ";
	for ( ; i != d.end() ; i++ ) {
		Log << l << (*i) << " ";
	}
	Log << l << endreq;

}

std::vector<Int_t> MediPixAlgo::GetColAsVector(Int_t col){

	std::vector<Int_t> colV;
	Int_t yDim = GetMatrixYdim();

	for(Int_t rowItr = 0; rowItr < yDim ; rowItr++)
	{
		colV.push_back(GetMatrixElement(col, rowItr));
	}

	return colV;
}

Int_t MediPixAlgo::GetEntriesPad(){
	return m_myManager->GetnEntriesPad();
}

Int_t MediPixAlgo::GetHitsInPad(){
	return m_myManager->GetnHitsInPad();
}
Int_t MediPixAlgo::GetChargeInPad(){ // alias of GetTotalTOT
	return m_myManager->GetnChargeInPad();
}
Int_t MediPixAlgo::GetTotalTOT(){
	return GetChargeInPad();
}
Bool_t MediPixAlgo::GetIsMCData(){
	return m_myManager->GetisMCData();
}

Bool_t MediPixAlgo::IsMCData(){
	return m_myManager->GetisMCData();
}

Int_t MediPixAlgo::GetFormat(){
	return m_myManager->GetfFormat();
}

Int_t MediPixAlgo::GetFrameWidth(){
	return m_myManager->GetfWidth();
}

Int_t MediPixAlgo::GetFrameHeight(){
	return m_myManager->GetfHeight();
}

Int_t MediPixAlgo::GetAcqMode(){
	return m_myManager->GetfAcq_mode();
}

Double_t MediPixAlgo::GetAcqTime(){
	return m_myManager->GetfAcq_time();
}

TString MediPixAlgo::GetAppliedFilters(){
	return m_myManager->GetfApplied_filters();
}

Double_t MediPixAlgo::GetAutoEraseInterval(){
	return m_myManager->GetfAuto_erase_interval();
}

Int_t MediPixAlgo::GetAutoeraseIntervalCounter(){
	return m_myManager->GetfAutoerase_interval_counter();
}

Bool_t MediPixAlgo::GetBSActive(){
	return m_myManager->GetfBS_active();
}

TString MediPixAlgo::GetChipboardID(){
	return m_myManager->GetfChipboardID();
}

Double_t MediPixAlgo::GetCoincLiveTime(){
	return m_myManager->GetfCoinc_live_time();
}

//UChar_t --> Byte_t
UChar_t MediPixAlgo::GetCoincidenceDelay(){
	return m_myManager->GetfCoincidence_delay();
}

//UChar_t --> Byte_t
UChar_t MediPixAlgo::GetCoincidenceMode(){
	return m_myManager->GetfCoincidence_mode();
}

std::vector<Int_t> MediPixAlgo::GetCounters(){
	return m_myManager->GetfCounters();
}

std::vector<Int_t> MediPixAlgo::GetDAQs(){
	return m_myManager->GetfDACs();
}

Int_t MediPixAlgo::GetTHL(){
	return m_myManager->GetTHL();
}

Double_t MediPixAlgo::GetHV(){ // same as GetBiasVoltage
	return m_myManager->GetfHV();
}

Double_t MediPixAlgo::GetBiasVoltage(){ // same as GetHV
	return m_myManager->GetfHV();
}

Int_t MediPixAlgo::GetHwTimer(){
	return m_myManager->GetfHw_timer();
}

TString MediPixAlgo::GetInterface(){
	return m_myManager->GetfInterface();
}

Double_t MediPixAlgo::GetMpxClock(){
	return m_myManager->GetfMpx_clock();
}

Int_t MediPixAlgo::GetMpxType(){
	return m_myManager->GetfMpx_type();
}

Int_t MediPixAlgo::GetPolarity(){
	return m_myManager->GetfPolarity();
}

Double_t MediPixAlgo::GetStartTime(){
	return m_myManager->GetfStart_time();
}

TString MediPixAlgo::GetStartTimeS(){
	return m_myManager->GetfStart_timeS();
}

Double_t MediPixAlgo::GetTimepixClock(){
	return m_myManager->GetfTimepix_clock();
}

Double_t MediPixAlgo::GetTriggerTime(){
	return m_myManager->GetfTrigger_time();
}

Long_t MediPixAlgo::GetFrameId(){
	return m_myManager->GetfFrameId();
}

TString MediPixAlgo::GetDataSetNumber(){
	return m_myManager->GetfMPXDataSetNumber();
}

Int_t MediPixAlgo::GetnEntriesChain(){ // same as GetNFrames
	return m_myManager->GetnEntriesChain();
}
Int_t MediPixAlgo::GetNFrames(){ // same as GetnEntriesChain
	return m_myManager->GetnEntriesChain();
}

//std::vector<TVector3> MediPixAlgo::GetPrimaryMCVertex(){
//return m_myManager->GetPrimaryMCVertex();
//};

// In the future I should write an error handler --> TODO
void MediPixAlgo::__FATALMC_ERROR(){
	cout << "[FATAL] Request of MC related information on real data ... giving up" << endl;
	exit(1);
}

Int_t MediPixAlgo::GetPrimaryMCVertex_N(){
	if( !IsMCData() ) { __FATALMC_ERROR(); }
	return m_myManager->GetPrimaryMCVertex_N();
}
Double_t MediPixAlgo::GetPrimaryMCVertex_X(int i){
	if( !IsMCData() ) { __FATALMC_ERROR(); }
	return m_myManager->GetPrimaryMCVertex_X(i);
}
Double_t MediPixAlgo::GetPrimaryMCVertex_Y(int i){
	if( !IsMCData() ) { __FATALMC_ERROR(); }
	return m_myManager->GetPrimaryMCVertex_Y(i);
}
Double_t MediPixAlgo::GetPrimaryMCVertex_Z(int i){
	if( !IsMCData() ) { __FATALMC_ERROR(); }
	return m_myManager->GetPrimaryMCVertex_Z(i);
}
/* ****************************************** */

TApplication * MediPixAlgo::GetApplication(){
	return m_myManager->GetApplication();
}

/* signals to the steer */
AlgoSignalsHandler MediPixAlgo::SignalToSteering(){
	return mySignals;
}

void MediPixAlgo::JumpToFrame(Long64_t jumpTo, Int_t direction_i){
	mySignals.nextFrame = jumpTo;
	mySignals.direction = direction_i;
}

void MediPixAlgo::ResetSignals(AlgoSignalsHandler newSignals){
	mySignals.nextFrame = newSignals.nextFrame;
}

///////////////////////////////////////////////////////////////
// Connection to StoreGate
Bool_t MediPixAlgo::PushToStoreGate(CandidateContainer * candidate, MPXDefs::Flags flg){ // fixing inconvenient member names
	return PullToStoreGateAccess(candidate, flg);
}
Bool_t MediPixAlgo::PullToStoreGateAccess(CandidateContainer * candidate, MPXDefs::Flags){
	m_myManager->GetStoreGate()->SaveObject(candidate);
	return true;
}

void MediPixAlgo::FillValuesForDisplay(Highlighter * hl, blob bl){

	hl->UploadDisplayValue("Number of Pixels     ", bl.bP.nPixels);
	hl->UploadDisplayValue("Inner Pixels         ", bl.bP.nInnerPixels);
	hl->UploadDisplayValue("Volume (TOT)         ", bl.bP.clusterTOT);
	hl->UploadDisplayValue("Cluster Energy       ", bl.bP.clusterEnergy);
	hl->UploadDisplayValue("Width x              ", bl.bP.width_x);
	hl->UploadDisplayValue("Width y              ", bl.bP.width_y);
	hl->UploadDisplayValue("Geo center x         ", bl.bP.geoCenter_x);
	hl->UploadDisplayValue("Geo center y         ", bl.bP.geoCenter_y);
	hl->UploadDisplayValue("Weighted center x    ", bl.bP.weightedCenter_x);
	hl->UploadDisplayValue("Weighted center y    ", bl.bP.weightedCenter_y);
	hl->UploadDisplayValue("Box area             ", bl.bP.boxArea);
	hl->UploadDisplayValue("Circle area          ", bl.bP.circleArea);
	hl->UploadDisplayValue("Ellipse area         ", bl.bP.ellipseArea);
	hl->UploadDisplayValue("  circle/elipse      ", bl.bP.circleArea/bl.bP.ellipseArea);
	hl->UploadDisplayValue("Elipse A             ", bl.bP.ellipseA);
	hl->UploadDisplayValue("Elipse B             ", bl.bP.ellipseB);
	hl->UploadDisplayValue("Rotation angle [deg] ", bl.bP.rotAngle * 180 / TMath::Pi());
	hl->UploadDisplayValue("Chisquare/dof        ", bl.bP.chisquare_OverDof);
	hl->UploadDisplayValue("Fit slope            ", bl.bP.fitSlope);
	hl->UploadDisplayValue("Fit cut              ", bl.bP.fitCut);
	hl->UploadDisplayValue("Balance to minimum   ", bl.bP.balanceToMin);
	hl->UploadDisplayValue("Balance to maximum   ", bl.bP.balanceToMin);
	hl->UploadDisplayValue("Min dist to center   ", bl.bP.minToGeoCenter);
	hl->UploadDisplayValue("Max dist to center   ", bl.bP.maxToGeoCenter);

	// upload the level1
	hl->UploadDisplayValue("lvl1                 ", bl.bP.lvl1);

	//hl->UploadDisplayValue("type of blob (string)", bl.GetTypeAsString());

}

// This is a User oriented function.  No consideration here for special objects
Int_t MediPixAlgo::GetNumberOfObjectsWithAuthor(TString authorName){
	// !!! WARNING, if there are no hits in the pad, Here we return 0 ! even if there
	// are service objects !!!
	if(GetHitsInPad() == 0 ) return 0;
	return m_myManager->GetStoreGate()->GetNObjWithAuthor(authorName);
}

Int_t MediPixAlgo::GetNumberOfSpecialObjects(MPXDefs::SpecialObjs objectType){
	return m_myManager->GetStoreGate()->GetNObjWithType(objectType);
}

Int_t MediPixAlgo::GetNumberOfSpecialObjectsLessEqualThan(MPXDefs::SpecialObjs objectType){
	return m_myManager->GetStoreGate()->GetNObjWithTypeLessEqualThan(objectType);
}

std::vector<CandidateContainer *> MediPixAlgo::GetObjectsSpecial(MPXDefs::SpecialObjs objectType){  
	return m_myManager->GetStoreGate()->GetObjsSpecial(objectType);
}

std::vector<CandidateContainer *> MediPixAlgo::GetObjectsSpecialWithAuthor(MPXDefs::SpecialObjs objectType
		, TString author){
	return m_myManager->GetStoreGate()->GetObjsSpecialWithAuthor(objectType, author);
}

CandidateContainer * MediPixAlgo::GetObjectFromAuthor(TString authorName, Int_t objIndex){
	return m_myManager->GetStoreGate()->GetObjFrom(authorName, objIndex);
}

///////////////////////////////////////////////////////////////
// Configuration values Viewer Connection Overloaded ... Class
//  ConfigurationValue is templated !
//template <class T> T MediPixAlgo::hola(T a){
//  return a;
//}
//template <class T> void MediPixAlgo::RegisterConfigurationValue(T * val, Char_t * valname){
// pull a 'ConfigurationValue' object to the store gate
//confValue_int = new ConfigurationValue<T>(this, val, valname);
//PullToStoreGateAccess(confValue_int, MPXDefs::DO_NOT_SERIALIZE_ME);
//}

void MediPixAlgo::BadConfigValue(const Char_t * valname){
	Log << MSG::ALWAYS << "Configuration value name can not have spaces ..." << endreq;
	Log << MSG::ALWAYS << "  I can't register this value : "<< valname << endreq;
}

void MediPixAlgo::RegisterConfigurationValue(Int_t * val, const Char_t * valname){
	// pull a 'ConfigurationValue' object to the store gate
	if( ((string)valname).find(" ") != string::npos){
		BadConfigValue(valname);
		return;
	}
	confValue_int = new ConfigurationValue<Int_t>(this, val, valname, MPXDefs::CONF_INT);
	PullToStoreGateAccess(confValue_int, MPXDefs::DO_NOT_SERIALIZE_ME);

	m_hasConfiguration = true;
}
void MediPixAlgo::RegisterConfigurationValue(Float_t * val, const Char_t * valname){
	// pull a 'ConfigurationValue' object to the store gate
	if( ((string)valname).find(" ") != string::npos ){
		BadConfigValue(valname);
		return;
	}
	confValue_float = new ConfigurationValue<Float_t>(this, val, valname, MPXDefs::CONF_FLOAT);
	PullToStoreGateAccess(confValue_float, MPXDefs::DO_NOT_SERIALIZE_ME);

	m_hasConfiguration = true;
}
void MediPixAlgo::RegisterConfigurationValue(Double_t * val, const Char_t * valname){
	// pull a 'ConfigurationValue' object to the store gate
	if( ((string)valname).find(" ") != string::npos ){
		BadConfigValue(valname);
		return;
	}
	confValue_double = new ConfigurationValue<Double_t>(this, val, valname, MPXDefs::CONF_DOUBLE);
	PullToStoreGateAccess(confValue_double, MPXDefs::DO_NOT_SERIALIZE_ME);

	m_hasConfiguration = true;
}

void MediPixAlgo::RegisterConfigurationValue(Bool_t * val, const Char_t * valname){
	// pull a 'ConfigurationValue' object to the store gate
	if( ((string)valname).find(" ") != string::npos ){
		BadConfigValue(valname);
		return;
	}
	confValue_bool = new ConfigurationValue<Bool_t>(this, val, valname, MPXDefs::CONF_BOOL);
	PullToStoreGateAccess(confValue_double, MPXDefs::DO_NOT_SERIALIZE_ME);

	m_hasConfiguration = true;
}

string MediPixAlgo::GetOutputBaseName ( ) {
	if(!m_myManager) {
		Log << MSG::ERROR << "You are calling GetOutputBaseName() at the wrong state (probably the constructor of your algo).  Move it to the Init()." << endreq;
		return string("temp");
	}
	return m_myManager->GetOutputBaseName();
}

TCanvas * MediPixAlgo::DrawInSeparateWindow(vector<TH1 *> g , MSG::Level vl){

	if(g.empty()) {
		Log << vl << "Nothing to draw ... return NULL pointer" << endreq;
		return 0x0;
	}
	TString cname; // take the cname from the first histo
	if(vl == MSG::DEBUG) {
		vector<TH1 *>::iterator i = g.begin();
		Log << vl << "Attempting to draw in a separate canvas the following objects : " << endreq;
		for( ; i != g.end() ; i++){
			if(!cname.Length()) {cname = (*i)->GetName();}
			Log << vl << (*i)->GetName() << " with " << (*i)->GetEntries() << " entries "<< endreq;
		}
	}
	///////////////////////////////////////////////////////////////////
	TCanvas * c11 = new TCanvas("TH1","TH1");
	int ngraphs = (int)g.size();
	int xdiv = floor( TMath::Sqrt(ngraphs) );
	int ydiv = ceil( (double)ngraphs/(double)xdiv );
	c11->Divide(xdiv, ydiv);

	vector<TH1 *>::iterator itr = g.begin();
	int cntr2 = 0;
	double posx = 0.;
	double posy = 0.;
	TString drawOpt = "";
	for ( ; itr != g.end(); itr++) {
		c11->cd(cntr2+1);
		(*itr)->Draw(drawOpt);
		(*itr)->GetXaxis()->SetLabelSize(0.1);
		(*itr)->GetYaxis()->SetLabelSize(0.1);

		posx = (*itr)->GetBinCenter( 1 );
		posy = (*itr)->GetMaximum() * 0.9;

		TLatex l1;
		l1.SetTextSize(0.12);
		l1.DrawLatex(posx, posy, (*itr)->GetTitle() );

		cntr2++;
		drawOpt = "same";
	}
	return c11;
}

TCanvas * MediPixAlgo::DrawInSeparateWindow(vector<TH2 *> g , MSG::Level vl){

	if(g.empty()) {
		Log << vl << "Nothing to draw ... return NULL pointer" << endreq;
		return 0x0;
	}
	TString cname; // take the cname from the first histo
	if(vl == MSG::DEBUG) {
		vector<TH2 *>::iterator i = g.begin();
		Log << vl << "Attempting to draw in a separate canvas the following objects : " << endreq;
		for( ; i != g.end() ; i++){
			if(!cname.Length()) {cname = (*i)->GetName();}
			Log << vl << (*i)->GetName() << " with " << (*i)->GetEntries() << " entries "<< endreq;
		}
	}
	///////////////////////////////////////////////////////////////////
	TCanvas * c11 = new TCanvas("TH2","TH2");
	int ngraphs = (int)g.size();
	int xdiv = floor( TMath::Sqrt(ngraphs) );
	int ydiv = ceil( (double)ngraphs/(double)xdiv );
	c11->Divide(xdiv, ydiv);

	vector<TH2 *>::iterator itr = g.begin();
	int cntr2 = 0;
	double posx = 0.;
	double posy = 0.;
	for ( ; itr != g.end(); itr++) {

		c11->cd(cntr2+1);
		(*itr)->Draw("colz");
		//(*itr)->SetT
		(*itr)->GetXaxis()->SetLabelSize(0.05);
		(*itr)->GetYaxis()->SetLabelSize(0.05);

		posx = (*itr)->GetBinCenter( 1 );
		posy = (*itr)->GetMaximum() * 0.9;

		TLatex l1;
		l1.SetTextSize(0.05);
		l1.DrawLatex(posx, posy, (*itr)->GetName() );

		cntr2++;
	}
	return c11;
}

TCanvas * MediPixAlgo::DrawInSeparateWindow(vector<TF2 *> g , MSG::Level vl){

	if(g.empty()) {
		Log << vl << "Nothing to draw ... return NULL pointer" << endreq;
		return 0x0;
	}

	if(vl == MSG::DEBUG) {
		vector<TF2 *>::iterator i = g.begin();
		Log << vl << "Attemping to draw in a separate canvas the following objects : " << endreq;
		for( ; i != g.end() ; i++){
			Log << vl << (*i)->GetName() << endreq;
		}
	}

	///////////////////////////////////////////////////////////////////
	TCanvas * c11 = new TCanvas("TF2","TF2");
	int ngraphs2D = (int)g.size();
	int xdiv = floor( TMath::Sqrt(ngraphs2D) );
	int ydiv = ceil( (double)ngraphs2D/(double)xdiv );
	c11->Divide(xdiv, ydiv);

	vector<TF2 *>::iterator itr = g.begin();
	int cntr2 = 0;
	for ( ; itr != g.end(); itr++) {
		c11->cd(cntr2+1);
		(*itr)->Draw("colz");

		cntr2++;
	}

	return c11;
}

TCanvas * MediPixAlgo::DrawInSeparateWindow(vector<TGraph2D *> g, MSG::Level vl, TString gopt, TString zaxis_opt){

	if(g.empty()) {
		Log << vl << "Nothing to draw ... return NULL pointer" << endreq;
		return 0x0;
	}

	if(vl == MSG::DEBUG) {
		vector<TGraph2D *>::iterator i = g.begin();
		Log << vl << "Attemping to draw in a separate canvas the following objects : " << endreq;
		for( ; i != g.end() ; i++){
			Log << vl << (*i)->GetName() << endreq;
		}
	}

	///////////////////////////////////////////////////////////////////
	TCanvas * c11 = new TCanvas("TGraph2D","TGraph2D");
	int ngraphs2D = (int)g.size();
	int xdiv = floor( TMath::Sqrt(ngraphs2D) );
	int ydiv = ceil( (double)ngraphs2D/(double)xdiv );
	c11->Divide(xdiv, ydiv);

	vector<TGraph2D *>::iterator itr = g.begin();
	int cntr2 = 0;
	for ( ; itr != g.end(); itr++) {
		c11->cd(cntr2+1);
		//(*itr)->Draw("tri1");
		(*itr)->Draw(gopt); // "tri1" by default
		(*itr)->GetXaxis()->SetLabelSize(0.1);
		(*itr)->GetYaxis()->SetLabelSize(0.1);
		(*itr)->GetZaxis()->SetLabelSize(0.08);
		if(zaxis_opt == "log" || zaxis_opt == "log") c11->GetPad(cntr2+1)->SetLogz();

		TLatex l1;
		l1.SetTextSize(0.12);
		l1.DrawLatex(0,0, (*itr)->GetName() );

		cntr2++;
	}


	return c11;
}

TCanvas * MediPixAlgo::DrawInSeparateWindow(vector<TGraph *> g , TString cname, MSG::Level vl){

	if(g.empty()) {
		Log << vl << "Nothing to draw ... return NULL pointer" << endreq;
		return 0x0;
	}

	if(vl == MSG::DEBUG) {
		vector<TGraph *>::iterator i = g.begin();
		Log << vl << "Attemping to draw in a separate canvas the following objects : " << endreq;
		for( ; i != g.end() ; i++){
			Log << vl << (*i)->GetName() << endreq;
		}
	}

	///////////////////////////////////////////////////////////////////
	TString canvasname = "TGraph2D_";
	canvasname += cname;
	TCanvas * c11 = new TCanvas(canvasname, canvasname);
	int ngraphs2D = (int)g.size();
	int xdiv = floor( TMath::Sqrt(ngraphs2D) );
	int ydiv = ceil( (double)ngraphs2D/(double)xdiv );
	c11->Divide(xdiv, ydiv);

	vector<TGraph *>::iterator itr = g.begin();
	int cntr2 = 0;
	for ( ; itr != g.end(); itr++) {
		c11->cd(cntr2+1);

		(*itr)->Draw("A*");
		(*itr)->SetMarkerStyle(20);
		//(*itr)->Draw("triw");
		(*itr)->GetXaxis()->SetLabelSize(m_MarkerFontSizeDrawingSeparateWindow);
		(*itr)->GetYaxis()->SetLabelSize(m_MarkerFontSizeDrawingSeparateWindow);

		TLatex l1;
		l1.SetTextSize(0.1);
		l1.DrawLatex(0,0, (*itr)->GetName() );

		cntr2++;
	}


	return c11;
}

TCanvas * MediPixAlgo::DrawInSeparateWindow(vector<vector< TGraph * > > g , TString cname, MSG::Level vl){

	if(g.empty()) {
		Log << vl << "Nothing to draw ... return NULL pointer" << endreq;
		return 0x0;
	}



	///////////////////////////////////////////////////////////////////
	TString canvasname = "TGraph2D_";
	canvasname += cname;
	TCanvas * c11 = new TCanvas(canvasname, canvasname);
	int ngraphs2D = (int)g.size();

	Log << MSG::ALWAYS << "NGRAPHS = " << ngraphs2D << endreq;

	int xdiv = floor( TMath::Sqrt(ngraphs2D) );
	int ydiv = ceil( (double)ngraphs2D/(double)xdiv );
	c11->Divide(xdiv, ydiv);

	vector<vector< TGraph * > >::iterator itr = g.begin();
	vector< TGraph * >::iterator i;

	int cntr2 = 0;
	TString drawOpt = "";
	int color = 50;
	for ( ; itr != g.end(); itr++) {

		c11->cd(cntr2+1);
		drawOpt = "A*l";
		color = 50;

		Log << MSG::ALWAYS << "Contours = " << (int)(*itr).size() << endreq;

		for( i = (*itr).begin() ; i != (*itr).end() ; i++ ) {

			(*i)->Draw(drawOpt); drawOpt = "*l";
			(*i)->SetMarkerStyle(7);
			(*i)->SetMarkerColor(color++);
			(*i)->GetXaxis()->SetLabelSize(m_MarkerFontSizeDrawingSeparateWindow);
			(*i)->GetYaxis()->SetLabelSize(m_MarkerFontSizeDrawingSeparateWindow);

		}

		//TLatex l1;
		//l1.SetTextSize(0.1);
		//l1.DrawLatex( 0, 0, (*(*itr).begin())->GetName() );

		cntr2++;
	}


	return c11;
}


TCanvas * MediPixAlgo::DrawInSeparateWindow(vector<TGraph *> g , vector<double> val, TString cname, MSG::Level vl){

	if ( g.empty() ) {
		Log << vl << "Nothing to draw ... return NULL pointer" << endreq;
		return 0x0;
	}

	if( g.size() != val.size() ) {
		Log << vl << "[ERROR] both vectors should have the same size" << endreq;
		return 0x0;
	}

	if(vl == MSG::DEBUG) {
		vector<TGraph *>::iterator i = g.begin();
		Log << vl << "Attemping to draw in a separate canvas the following objects : " << endreq;
		for( ; i != g.end() ; i++){
			Log << vl << (*i)->GetName() << endreq;
		}
	}

	///////////////////////////////////////////////////////////////////
	TString canvasname = "TGraph2D_";
	canvasname += cname;
	TCanvas * c11 = new TCanvas(canvasname, canvasname);
	int ngraphs2D = (int)g.size();
	int xdiv = floor( TMath::Sqrt(ngraphs2D) );
	int ydiv = ceil( (double)ngraphs2D/(double)xdiv );
	c11->Divide(xdiv, ydiv);

	vector<TGraph *>::iterator itr = g.begin();
	vector<double>::iterator itrVal = val.begin();

	int cntr2 = 0;
	for ( ; itr != g.end(); itr++) {
		c11->cd(cntr2+1);

		(*itr)->Draw("A*");
		(*itr)->SetMarkerStyle(20);
		//(*itr)->Draw("triw");
		(*itr)->GetXaxis()->SetLabelSize(m_MarkerFontSizeDrawingSeparateWindow);
		(*itr)->GetYaxis()->SetLabelSize(m_MarkerFontSizeDrawingSeparateWindow);

		TLatex l1;
		l1.SetTextSize(0.2);
		TString valS = "";
		valS += TString::Format("%.2f",*itrVal);
		l1.DrawLatex((*itr)->GetXaxis()->GetXmax()*0.2, TMath::Abs( (*itr)->GetYaxis()->GetXmax() )*0.8, (*itr)->GetName() );
		l1.DrawLatex((*itr)->GetXaxis()->GetXmax()*0.7, TMath::Abs( (*itr)->GetYaxis()->GetXmax() )*0.7, valS );

		itrVal++;

		cntr2++;
	}


	return c11;
}


/**
 * Both vector need to have the same length
 */
TCanvas * MediPixAlgo::DrawInSeparateWindow(vector<TGraph *> g, vector<TF1 *> f, TString cname, MSG::Level vl){

	if(g.size() != f.size()) {
		Log << MSG::INFO << "The vector don't have the same size.  Giving up.  Return null pointer of TCanvas." << endreq;
		return 0x0;
	}

	if(g.empty()) {
		Log << vl << "Nothing to draw ... return NULL pointer" << endreq;
		return 0x0;
	}

	if(vl == MSG::DEBUG) {
		vector<TGraph *>::iterator i = g.begin();
		Log << vl << "Attemping to draw in a separate canvas the following objects : " << endreq;
		for( ; i != g.end() ; i++){
			Log << vl << (*i)->GetName() << endreq;
		}
	}

	///////////////////////////////////////////////////////////////////
	TString canvasname = "TGraph2D_";
	canvasname += cname;
	TCanvas * c11 = new TCanvas(canvasname, canvasname);
	int ngraphs2D = (int)g.size();
	int xdiv = floor( TMath::Sqrt(ngraphs2D) );
	int ydiv = ceil( (double)ngraphs2D/(double)xdiv );
	c11->Divide(xdiv, ydiv);

	vector<TGraph *>::iterator itr = g.begin();
	vector<TF1 *>::iterator itr_f = f.begin();

	int cntr2 = 0;
	for ( ; itr != g.end(); itr++) {
		c11->cd(cntr2+1);

		(*itr)->Draw("A*");
		(*itr)->SetMarkerStyle(20);
		//(*itr)->Draw("triw");
		(*itr)->GetXaxis()->SetLabelSize(0.08);
		(*itr)->GetYaxis()->SetLabelSize(0.08);

		TLatex l1;
		l1.SetTextSize(0.1);
		l1.DrawLatex(0,0, (*itr)->GetName() );

		// TF1 on top
		(*itr_f)->Draw("same");

		(*itr)->GetYaxis()->SetRangeUser( 0, (*itr)->GetYaxis()->GetXmax() );

		itr_f++;
		cntr2++;
	}


	return c11;
}

/**
 * Both vector need to have the same length
 */
TCanvas * MediPixAlgo::DrawInSeparateWindow(vector<TH1 *> g, vector<TF1 *> f, TString cname, MSG::Level vl){

	if (g.size() != f.size()) {
		Log << MSG::INFO << "The vector don't have the same size.  Giving up.  Return null pointer of TCanvas." << endreq;
		return 0x0;
	}

	if (g.empty()) {
		Log << vl << "Nothing to draw ... return NULL pointer" << endreq;
		return 0x0;
	}

	if (vl == MSG::DEBUG) {
		vector<TH1 *>::iterator i = g.begin();
		Log << vl << "Attemping to draw in a separate canvas the following objects : " << endreq;
		for( ; i != g.end() ; i++){
			Log << vl << (*i)->GetName() << endreq;
		}
	}

	///////////////////////////////////////////////////////////////////
	TString canvasname = "TGraph2D_";
	canvasname += cname;
	TCanvas * c11 = new TCanvas(canvasname, canvasname);
	int ngraphs2D = (int)g.size();
	int xdiv = floor( TMath::Sqrt(ngraphs2D) );
	int ydiv = ceil( (double)ngraphs2D/(double)xdiv );
	c11->Divide(xdiv, ydiv);

	vector<TH1 *>::iterator itr = g.begin();
	vector<TF1 *>::iterator itr_f = f.begin();

	int cntr2 = 0;
	for ( ; itr != g.end(); itr++) {
		c11->cd(cntr2+1);

		(*itr)->Draw("");
		(*itr)->SetMarkerStyle(20);
		//(*itr)->Draw("triw");
		(*itr)->GetXaxis()->SetLabelSize(0.08);
		(*itr)->GetYaxis()->SetLabelSize(0.08);

		TLatex l1;
		l1.SetTextSize(0.1);
		l1.DrawLatex(0,0, (*itr)->GetName() );

		// TF1 on top
		(*itr_f)->Draw("same");

		//(*itr)->GetYaxis()->SetRangeUser( 0, (*itr)->GetYaxis()->GetXmax() );

		itr_f++;
		cntr2++;
	}


	return c11;
}


#endif
