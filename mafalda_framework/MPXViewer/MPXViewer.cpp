/*
 * 	Copyright 2009 John Idarraga
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

#ifndef MPXViewer_cpp
#define MPXViewer_cpp

#include "MPXViewer.h"

class MPXViewerControl;

using namespace MSG;

ClassImp(MPXViewer)

MPXViewer::MPXViewer(){

	//g_theApp = new TApplication("Output", 0, NULL);

	Saturation = -1;
	Maxlightness = -1;
	Minlightness = -1;
	MaxHue = -1; // from 0 to 360
	MinHue = -1;
	Lightness = -1;

	setHistoPalette(VIS_MOD_JET);
	//setHistoPalette(VIS_MOD_GRAY);

	nHitsInPadCut = 1;
	nChargeInPadCut = 1;

	/* holds metadata to show in Viewer */
	updateInfo = new ControlUpdateInfo;
	updateInfo->frameId = -1;
	updateInfo->frameIdLocal = -1;
	updateInfo->currentFile = "";
	updateInfo->currentFileId = -1;
	updateInfo->nEntries = -1;
	updateInfo->nEntriesLocal = -1;
	updateInfo->acqTime = 0;
	updateInfo->nHitsInPad = -1;
	updateInfo->nChargeInPad = -1;
	updateInfo->vTHL = -1;
	updateInfo->startTimeS = "";
	updateInfo->startTimeD = 0.;

	feedBackInfo_New = new ControlFeedbackInfo;
	feedBackInfo_New->maxHisto = 0;
	feedBackInfo_New->minHisto = 0;
	feedBackInfo_New->autoAdjust = true;
	feedBackInfo_New->m_plotTOT = true;
	feedBackInfo_New->m_plotCalib = false;
	feedBackInfo_New->m_plotLvl1 = false;
	feedBackInfo_New->m_plotMC = false;
	feedBackInfo_New->m_plotTruthMC = false;

	feedBackInfo_Old = new ControlFeedbackInfo;
	feedBackInfo_Old->maxHisto = 0;
	feedBackInfo_Old->minHisto = 0;
	feedBackInfo_Old->autoAdjust = true;
	feedBackInfo_Old->m_plotTOT = true;
	feedBackInfo_Old->m_plotCalib = false;
	feedBackInfo_Old->m_plotLvl1 = false;
	feedBackInfo_Old->m_plotMC = false;
	feedBackInfo_Old->m_plotTruthMC = false;

	/* steering */
	vSteer = new ViewerSteer;
	vSteer->direction = SEEK_FORWARD;
	vSteer->oldDirection = SEEK_FORWARD;

	/* signals */
	signalGenerated = false;

	/* a title */
	m_viewerTitle = "";
}

void MPXViewer::SetCuts(Int_t nHitsInPadCut_i, Int_t nChargeInPadCut_i){
	nHitsInPadCut = nHitsInPadCut_i;
	nChargeInPadCut = nChargeInPadCut_i;
}

void MPXViewer::Init(){

	c1 = new TCanvas("c", "MAFalda Viewer", 700, 600);
	//b1 = new TBrowser();
	//c1->Browse(b1);
	//ConnectTBrowser(b1);

	c1->Show();

	///////////////////////////////////////////////////////////////////////////////
	// Fetch from the StoreGate The configuration values
	//  Those objects are in the StoreGate with special type MPXDefs::SpecialObjs
	std::vector<CandidateContainer *> specialObjs_conf_int = GetObjectsSpecial(MPXDefs::CONF_INT);
	std::vector<CandidateContainer *>::iterator cVItr_int = specialObjs_conf_int.begin();

	std::vector<CandidateContainer *> specialObjs_conf_float = GetObjectsSpecial(MPXDefs::CONF_FLOAT);
	std::vector<CandidateContainer *>::iterator cVItr_float = specialObjs_conf_float.begin();

	std::vector<CandidateContainer *> specialObjs_conf_double = GetObjectsSpecial(MPXDefs::CONF_DOUBLE);
	std::vector<CandidateContainer *>::iterator cVItr_double = specialObjs_conf_double.begin();

	if(!specialObjs_conf_int.empty() || !specialObjs_conf_float.empty() || !specialObjs_conf_double.empty())
	{
		Log << MSG::DEBUG << "Registered Algorithms' Configuration : " << endreq;
		for( ; cVItr_double != specialObjs_conf_double.end() ; cVItr_double++)
			Log << MSG::DEBUG << "  " << ((ConfigurationValue<Double_t> *)(*cVItr_double))->GetAuthor()
			<< " : " << ((ConfigurationValue<Double_t> *)(*cVItr_double))->GetConfigName().c_str()
			<< "<double> = " << ((ConfigurationValue<Double_t> *)(*cVItr_double))->GetConfigValue() << endreq;


		for( ; cVItr_float != specialObjs_conf_float.end() ; cVItr_float++)
			Log << MSG::DEBUG << "  " << ((ConfigurationValue<Float_t> *)(*cVItr_float))->GetAuthor()
			<< " : " << ((ConfigurationValue<Float_t> *)(*cVItr_float))->GetConfigName().c_str()
			<< "<float> = " << ((ConfigurationValue<Float_t> *)(*cVItr_float))->GetConfigValue() << endreq;


		for( ; cVItr_int != specialObjs_conf_int.end() ; cVItr_int++)
			Log << MSG::DEBUG << "  " << ((ConfigurationValue<Int_t> *)(*cVItr_int))->GetAuthor()
			<< " : " << ((ConfigurationValue<Int_t> *)(*cVItr_int))->GetConfigName().c_str()
			<< "<int> = " << ((ConfigurationValue<Int_t> *)(*cVItr_int))->GetConfigValue() << endreq;

	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////


	std::vector<CandidateContainer *> specialObjs_conf[3] = {specialObjs_conf_int, specialObjs_conf_float,
			specialObjs_conf_double };

	int sizex = GetMatrixXdim();
	int sizey = GetMatrixYdim();
	pMenu = new MPXViewerControl(GetApplication(), vSteer, specialObjs_conf, sizex, sizey);

}

Bool_t MPXViewer::SearchHotPixel()
{

	///////////// TODO, TO BE REMOVED !!!

	Int_t nHotPixels = 0;
	for(Int_t rowItr = 0 ; rowItr < GetMatrixYdim() ; rowItr++)
	{
		for(Int_t colItr = 0 ; colItr < GetMatrixXdim() ; colItr++)
		{
			if(GetMatrixElement(colItr, rowItr) > 600)
			{
				nHotPixels++;
			}
			if(nHotPixels > 4)
				return false;
		}
	}

	return true;
}

void MPXViewer::Execute(){

	///////////////////////////////////////////////////////////////////////////////
	// Get extra object to draw, things of top of the matrix
	//  Those objects are in the StoreGate with special type MPXDefs::SpecialObjs
	std::vector<CandidateContainer *> specialObjs = GetObjectsSpecial(MPXDefs::VIS);
	Log << MSG::DEBUG << "Drawable objects --> " << (Int_t)specialObjs.size() << endreq;

	std::vector<CandidateContainer *> specialObjsSkip = GetObjectsSpecial(MPXDefs::VIS_SKIP);
	Log << MSG::DEBUG << "Objects signaling to skip --> " << (Int_t)specialObjsSkip.size() << endreq;

	// Extra drawing.  Special objects
	bool skipFromAlgo = false;
	if ( (int)specialObjsSkip.size() > 0 ) skipFromAlgo = true;

	TH2I * h_frame = 0;

	Bool_t overflowFlg = false;
	if(GetHitsInPad() < 0){
		Log << MSG::WARNING << "HitsInPad is less than zero ! ... probably variable overflow ocurred !" << endreq;
		overflowFlg = true;
	}

	if(GetChargeInPad() < 0){
		Log << MSG::WARNING << "ChargeInPad is less than zero ! ... probably variable overflow ocurred !" << endreq;
		overflowFlg = true;
	}

	if ( !overflowFlg && (GetHitsInPad() < nHitsInPadCut || GetChargeInPad() < nChargeInPadCut) ) {
		Log << MSG::DEBUG << "Too few active pixels in the frame ... skipping" << endreq;
		Log << MSG::DEBUG << "Cut set at (" << nHitsInPadCut << ", " <<
				nChargeInPadCut << ") (hits, TOT/counts) " << endreq;

	} else if ( skipFromAlgo ) { // Signal sent from an algorithm indicating to skip for any other reason

		Log << MSG::DEBUG << "Signal sent from : " << " | signaling to skip frame " << endreq;

	} else {

		/* get histogram to display */
		h_frame = getFrameHist();

		/* sending info to the menu */
		updateInfo->frameId = GetGlobalCounters();
		updateInfo->frameIdLocal = GetFrameId();
		updateInfo->currentFile = FindCurrentFile();
		updateInfo->currentFileId = FindCurrentFileId(updateInfo->currentFile);
		updateInfo->nEntries = GetNFrames(); // total !
		updateInfo->nEntriesLocal = GetNEntriesForFile(updateInfo->currentFile); // per file
		updateInfo->acqTime = GetAcqTime();
		updateInfo->nHitsInPad = GetHitsInPad();
		updateInfo->nChargeInPad = GetChargeInPad();
		updateInfo->vTHL = GetTHL();
		updateInfo->startTimeS = GetStartTimeS();
		updateInfo->startTimeD = GetStartTime();

		/* draw the histogram */
		c1->cd();
		// Draw !!
		h_frame->Draw("colz");
		// Interesting.  The object 'palette' of the class TPaletteAxis is created only
		// when I update the canvas, not before.  Thanks gdb[backtrace] ;)
		c1->Update();

		// apply properties
		histoProperties(h_frame, c1);

		// Custom Title
		TString title = "frame ";
		title += GetFrameId();

		// if the frame is MPXGeant4 data
		if(GetIsMCData()){
			title += " | Geant4 data";
		}

		TLatex * l1 = new TLatex();
		l1->SetTextSize(0.03);
		UInt_t yposT = GetMatrixYdim();
		UInt_t xposT = GetMatrixXdim();

		l1->DrawLatex(0, yposT + 10, title);
		if(m_viewerTitle.Length()){
			l1->DrawLatex(xposT/2., yposT + 10, m_viewerTitle);
		}

		if ( feedBackInfo_New->m_plotTOT ) {
			l1->DrawLatex(xposT, yposT*1.05, "TOT");
		} else if ( feedBackInfo_New->m_plotCalib ) {
			l1->DrawLatex(xposT, yposT*1.05, "Energy");
			l1->DrawLatex(xposT, yposT*1.01, "[keV]");
		} else if ( feedBackInfo_New->m_plotLvl1 ) {
			l1->DrawLatex(xposT, yposT*1.05, "LVL1");
		} else if ( feedBackInfo_New->m_plotMC ){
			l1->DrawLatex(xposT, yposT*1.05, "MC");
			l1->DrawLatex(xposT, yposT*1.01, "[keV]");
		} else if ( feedBackInfo_New->m_plotTruthMC ){
			l1->DrawLatex(xposT, yposT*1.05, "TruthMC");
			l1->DrawLatex(xposT, yposT*1.01, "[keV]");
		}

		// Extra drawing.  Special objects
		std::vector<CandidateContainer *>::iterator highItr = specialObjs.begin();
		for ( ; highItr != specialObjs.end() ; highItr++) {
			((Highlighter *)(*highItr))->GetDrawable()->Draw();
		}

		// Interesting.  The object 'palette' of the class TPaletteAxis is created only
		// when I update the canvas, not before.  Thanks gdb !!! :)
		c1->Update();

		Log << MSG::INFO << "Global frame Id: "
				<< updateInfo->frameId << "/"
				<< updateInfo->nEntries
				<< "  [ Indexes start at zero (frameIndex/Total) ]"
				<< endreq;
		Log << MSG::INFO << "Current file ["
				<< updateInfo->currentFileId
				<< "/" << GetTotalNumberOfFiles()
				<< "] \""
				<< updateInfo->currentFile
				<< "\" | Local Id : "
				<< updateInfo->frameIdLocal << "/"
				<< updateInfo->nEntriesLocal
				<< endreq;
		Log << MSG::INFO << "|------------------- Displaying frame: "
				<< GetFrameId()
				<< " -------------------|"
				<< endreq;

		/* update properties */
		feedBackInfo_New = pMenu->ControlUpdate(updateInfo, feedBackInfo_Old);

		GetApplication()->Run(kTRUE);

		/* if properties changed, update ! */
		if(feedBackInfo_Old != feedBackInfo_New){
			feedBackInfo_Old = feedBackInfo_New;
		}

	}

	// Check for a request of disconnection through the ViewerControl
	if ( feedBackInfo_Old->m_disconnect ) {
		DisconnectMeFromRun();
	}

	/* steer jump */
	if(vSteer->direction == JUMP)
	{
		/* send a signal*/
		signalGenerated = true;
		/* set the change (nextFrame, direction) */
		JumpToFrame(feedBackInfo_Old->jumpTo, vSteer->oldDirection);
		vSteer->direction = vSteer->oldDirection;
	}
	/* steer back or forward */
	else if(vSteer->oldDirection == SEEK_FORWARD &&
			vSteer->direction == SEEK_BACKWARDS)
	{
		/* send a signal*/
		signalGenerated = true;
		/* set the change (nextFrame, direction) */
		JumpToFrame(DONT_JUMP, SEEK_BACKWARDS);
		vSteer->oldDirection = vSteer->direction;
	}
	else if(vSteer->oldDirection == SEEK_BACKWARDS &&
			vSteer->direction == SEEK_FORWARD)
	{
		/* send a signal*/
		signalGenerated = true;
		/* set the change (nextFrame, direction) */
		JumpToFrame(DONT_JUMP, SEEK_FORWARD);
		vSteer->oldDirection = vSteer->direction;
	}
	else if(vSteer->oldDirection == SEEK_FORWARD &&
			vSteer->direction == SEEK_FORWARD)
	{
		/* send a signal*/
		signalGenerated = true;
		/* set the change (nextFrame, direction) */
		JumpToFrame(DONT_JUMP, SEEK_FORWARD);
		vSteer->oldDirection = vSteer->direction;
	}
	else if(vSteer->oldDirection == SEEK_BACKWARDS &&
			vSteer->direction == SEEK_BACKWARDS)
	{
		/* send a signal*/
		signalGenerated = true;
		/* set the change (nextFrame, direction) */
		JumpToFrame(DONT_JUMP, SEEK_BACKWARDS);
		vSteer->oldDirection = vSteer->direction;
	}
	else if(vSteer->direction == DONT_MOVE)
	{
		/* send a signal*/
		signalGenerated = true;
		/* set the change (nextFrame, direction) */
		JumpToFrame(DONT_JUMP, DONT_MOVE);
	}

	if(h_frame)
		delete h_frame;
}

/**  
 *  This function is extremely inneficient but we are only calling
 *  this when the the viewer is ON.  This never executes on a batch
 *  run
 */

TString MPXViewer::FindCurrentFile(){

	Int_t globalCntr = GetGlobalCounters();
	map<TString, Int_t> ifM = GetInputFilesMap();

	TString f;
	map<TString, Int_t>::iterator itr = ifM.begin();

	Int_t filesCntr = 0;
	for( ; itr != ifM.end() ; itr++)
	{
		filesCntr += (*itr).second;
		if( globalCntr >= filesCntr )
		{
			continue;
		}
		else
		{
			return (*itr).first;
		}
	}

	return TString("");
}

Int_t MPXViewer::FindCurrentFileId(TString f){

	map<TString, Int_t> ifM = GetInputFilesMap();
	map<TString, Int_t>::iterator itr = ifM.begin();

	Int_t filesCntr = 0;
	for( ; itr != ifM.end() ; itr++){

		if(f == (*itr).first)
			return filesCntr;

		filesCntr++;
	}

	return -1; // didn't find the file ... shouldn't happen.
}

void MPXViewer::Finalize(){


}

TH2I * MPXViewer::getFrameHist(){

	TString title = "frame ", hname = "frame_";
	title += GetFrameId(); hname += GetFrameId();

	if(GetIsMCData()){ // if the frame is MPXGeant4 data
		title += " | MPXGeant4";
		hname += "_MPXGeant4";
	}

	UInt_t ysizeT = GetMatrixYdim();
	UInt_t xsizeT = GetMatrixXdim();

	TH2I * histFrame = new TH2I(hname, title,
			xsizeT, 0, xsizeT,
			ysizeT, 0, ysizeT);

	// don't run throught the whole matrix
	// bring back the map<int,int> and fill only what's needed !!!
	int ei;
	for(UInt_t rowItr = 0 ; rowItr < ysizeT ; rowItr++)
	{
		for(UInt_t colItr = 0 ; colItr < xsizeT ; colItr++)
		{

			if ( feedBackInfo_New->m_plotTOT ) {  // plot TOT or counts

				histFrame->Fill(colItr, rowItr, GetMatrixElement(colItr, rowItr));

			}  else if ( feedBackInfo_New->m_plotCalib ) {  // plot Energy through Calib

				// Get the energy from calibration, only if available in the run
				ei = GetCalibEnergy(colItr, rowItr);

				if( ei <= 0 ){

					histFrame->Fill(colItr, rowItr, 0);

				} else {

					//Log << MSG::INFO << "Energy = " << ei << endreq;
					histFrame->Fill(colItr, rowItr, ei);

				}

			} else if ( feedBackInfo_New->m_plotLvl1 ) {  // plot LVL1 info

				// lvl1 case, lvl1 is put to -1 if there is no lvl1 info available
				if(GetLVL1(colItr, rowItr) <= 0){
					histFrame->Fill(colItr, rowItr, 0);
				}else{
					histFrame->Fill( colItr, rowItr, GetLVL1(colItr, rowItr) );
				}

			} else if ( feedBackInfo_New->m_plotTruthMC ) { // plot MC info

				histFrame->Fill(colItr, rowItr, GetMatrixElementMCEdep(colItr, rowItr));

			}



		}
	}

	return histFrame;
}

void MPXViewer::setHistoPalette(Int_t visMode) {

	// Palettes prepared by Olivier.
	// That's the part I (Olivier) am adding.

	switch (visMode) {

	case VIS_MOD_GRAY:
		Saturation = GRAY_SATURATION_VAL;
		Maxlightness = GRAY_MAXLIGHTNESS_VAL;
		Minlightness = GRAY_MINLIGHTNESS_VAL;
		MaxHue = GRAY_MAXHUE_VAL;
		MinHue = GRAY_MINHUE_VAL;
		break;

	case VIS_MOD_JET:
		Saturation = JET_SATURATION_VAL;
		Maxlightness = JET_MAXLIGHTNESS_VAL;
		Minlightness = JET_MINLIGHTNESS_VAL;
		MaxHue = JET_MAXHUE_VAL;
		MinHue = JET_MINHUE_VAL;
		break;

	case VIS_MOD_HOT:
		Saturation = HOT_SATURATION_VAL;
		Maxlightness = HOT_MAXLIGHTNESS_VAL;
		Minlightness = HOT_MINLIGHTNESS_VAL;
		MaxHue = HOT_MAXHUE_VAL;
		MinHue = HOT_MINHUE_VAL;
		break;

	case VIS_MOD_COOL:
		Saturation = COOL_SATURATION_VAL;
		Maxlightness = COOL_MAXLIGHTNESS_VAL;
		Minlightness = COOL_MINLIGHTNESS_VAL;
		MaxHue = COOL_MAXHUE_VAL;
		MinHue = COOL_MINHUE_VAL;
		break;

	default:

		break;
	}

	int          index;
	float        saturation, lightness, hue, r, g, b; // rv, gv, bv;
	TColor     * color;

	for (Int_t i = 0 ; i < MAX_COLORS ; i++)
	{
		index = palette[i] = MAX_COLORS+1+i;

		// Check if the color already exists
		color = gROOT->GetColor(index);
		if(!color)
			color = new TColor(index, 0, 0, 0);

		lightness = Maxlightness-(i+1)*((Maxlightness-Minlightness)/MAX_COLORS);
		hue = MaxHue-(i+1)*((MaxHue-MinHue)/MAX_COLORS);
		saturation = Saturation;

		color->HLStoRGB(hue, lightness, saturation, r, g, b);
		color->SetRGB(r, g, b);

	}

	gStyle->SetPalette(MAX_COLORS, palette);

}

void MPXViewer::histoProperties(TH2I * h1, TCanvas * c){

	// no status
	h1->SetStats(false);

	// seting min max in auto-adjust or manual mode
	if(feedBackInfo_Old->autoAdjust)
	{
		Log << MSG::DEBUG << "Frame min and max values are set to auto-adjust mode" << endreq;
		feedBackInfo_Old->maxHisto = (Int_t)h1->GetMaximum();
		feedBackInfo_Old->minHisto = (Int_t)h1->GetMinimum();

		if(feedBackInfo_Old->minHisto == feedBackInfo_Old->maxHisto){
			feedBackInfo_Old->minHisto--;
		}

	}
	else
	{
		Log << MSG::DEBUG << "Frame min and max values are set to manual mode" << endreq;

		// Handle bad values.  They can't be negative, guaranteed by
		// TGNumberEntry object.  Handle only equal values
		if(feedBackInfo_Old->minHisto == feedBackInfo_Old->maxHisto){
			h1->SetMinimum(feedBackInfo_Old->minHisto);
			feedBackInfo_Old->maxHisto = feedBackInfo_Old->minHisto + 1;
			Log << MSG::WARNING << "[WARNING] Same values for min and max in histogram ! trying to fix (max=min+1)" << endreq;
		}

	}

	h1->SetMaximum(feedBackInfo_Old->maxHisto);
	h1->SetMinimum(feedBackInfo_Old->minHisto);

	c->ToggleEventStatus();
	c->SetHighLightColor(2);
	c->SetBorderMode(0);

	c->SetFrameBorderMode(0);

	// Background of the graph
	//c->SetFrameFillColor(palette[0]);
	c->SetFillColor(kWhite);

	c->SetBorderMode(0);

	// fix palette
	TPaletteAxis * thePalette = (TPaletteAxis *)h1->GetListOfFunctions()->FindObject("palette");
	if(thePalette){
		thePalette->SetLabelSize(0.03);
		double posx2 = thePalette->GetX2NDC();
		double posx1 = thePalette->GetX1NDC();
		thePalette->SetX2NDC((posx2 + posx1)/2.);
		// The label size can be done from the style styles/MedipixStyle.C
		//thePalette->SetLabelSize(0.03);
	}
	c->Update();

}

AlgoSignalsHandler MPXViewer::SignalToSteering(){
	Log << MSG::LOOP_DEBUG << "Sending signal to steering !" << endreq;
	return GetSignals();
}

Bool_t MPXViewer::SignalFlag(){

	if(!signalGenerated)
		return signalGenerated;
	else
		signalGenerated = false;
	return true;

}

#endif
