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

#ifndef ViewerControl_cpp
#define ViewerControl_cpp

#include "ViewerControl.h"
#include "MPXViewer.h"

#include <stdio.h>
#include <time.h>

ClassImp(MPXViewerControl)

MPXViewerControl::MPXViewerControl(TApplication * g_theApp_i, ViewerSteer * vSteer_i, std::vector<CandidateContainer *> * objs_conf, int sizex, int sizey)
: TGMainFrame(gClient->GetRoot(), 1, 2, kHorizontalFrame)
{

	g_theApp = g_theApp_i;

	// Set the Log service
	Log.setAlgoName("MPXViewerControl");
	Log.OutputLevel = MSG::INFO;

	/////////////////////////////
	m_objs_conf_int = objs_conf[0];
	m_objs_conf_float = objs_conf[1];
	m_objs_conf_double = objs_conf[2];
	//m_objs_conf_bool = 0; // objs_conf[3];
	m_objs_conf_itr = 0; // manage position on vector m_objs_conf

	//ShowResults();
	Build();
	m_p = 0;

	feedBackInfoControl = new ControlFeedbackInfo;
	feedBackInfoControl->maxHisto = 1;
	feedBackInfoControl->minHisto = 0;
	feedBackInfoControl->autoAdjust = true;
	feedBackInfoControl->jumpTo = 0;
	feedBackInfoControl->m_plotTOT = true;
	feedBackInfoControl->m_plotLvl1 = false;
	feedBackInfoControl->m_plotMC = false;
	feedBackInfoControl->m_plotTruthMC = false;
	feedBackInfoControl->m_disconnect = false;

	vSteerControl = vSteer_i;

	m_matrixSizeX = sizex;
	m_matrixSizeY = sizey;

}

MPXViewerControl::~MPXViewerControl(){

}

void MPXViewerControl::Build(){

	SetCleanup(kDeepCleanup);

	// Are there any parameters coming from MPXAlgo's ?
	Bool_t thereAreAlgoParameters = ( !m_objs_conf_double.empty()
			|| !m_objs_conf_float.empty()
			|| !m_objs_conf_int.empty() );
	// Find out how many algorithms are registered
	if(thereAreAlgoParameters) GenerateAlgorithmMap();

	///////////////////////////////////////////////////////////////////////////////
	// Frames

	// Algorithms
	TGVerticalFrame * contentsAlgo = new TGVerticalFrame(this);
	AddFrame(contentsAlgo, new TGLayoutHints(kLHintsRight | kLHintsExpandY, 5, 5, 5, 5));

	if(thereAreAlgoParameters){
		// Separator 1 ////////
		TGVertical3DLine *separator1 = new TGVertical3DLine(this);
		AddFrame(separator1, new TGLayoutHints(kLHintsRight | kLHintsExpandY));  // Separator 1
	}

	// MetaData
	TGVerticalFrame *contentsMeta = new TGVerticalFrame(this);
	AddFrame(contentsMeta, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 5, 5));

	///////////////////////////////////
	// Algos

	// number of TGGroupFrame objects = number of algorithms with configuration
	Int_t nGroups = (Int_t)m_algorithmMap.size();
	TGGroupFrame * algoGroup[nGroups];

	// Number of configuration values needing a TextMargin object
	Int_t nTextMargin = (Int_t)m_objs_conf_int.size()
    				+ (Int_t)m_objs_conf_float.size()
    				+ (Int_t)m_objs_conf_double.size();
	fconfBoxes = new TextMargin*[nTextMargin];

	std::map<TString, Int_t>::iterator aItr = m_algorithmMap.begin();
	Int_t cntrAlgos = 0;
	Int_t cntrBoxes = 0;
	TString label = "";

	// check the order of filling here
	// 1) float
	// 2) double
	// 3) int

	Int_t nConf = 0;
	for( ; aItr != m_algorithmMap.end() ; aItr++){

		algoGroup[cntrAlgos] = new TGGroupFrame(contentsAlgo, (*aItr).first);
		algoGroup[cntrAlgos]->SetTitlePos(TGGroupFrame::kCenter);

		// (*aItr).second, is the number of configuration values in this algorithm
		nConf += (*aItr).second;

		//////////////////////////////////////////////////////
		/* Composite frame for algo config <Float_t> information, only if requested */
		if(!m_objs_conf_float.empty()) {

			for(int i = 0 ; i < (Int_t)m_objs_conf_float.size() ; i++) {

				// only if it belongs to this algorith
				if(((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetAuthor() == (*aItr).first)
				{
					label = ((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetConfigName();
					// create instance
					fconfBoxes[cntrBoxes] = new TextMargin(algoGroup[cntrAlgos], label);
					// build unique identifier
					TString mappingName = BuildKeyNameTextMargin((*aItr).first, label);
					// push into map for further reference
					m_confMapping[mappingName] = fconfBoxes[cntrBoxes];

					algoGroup[cntrAlgos]->AddFrame(fconfBoxes[cntrBoxes], new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));

					fconfBoxes[cntrBoxes]->GetTGNumberEntry()->Connect("ValueSet(Long_t)", "MPXViewerControl",
							(MPXViewerControl *)this, "SetConfigAlgoValue()");
					fconfBoxes[cntrBoxes]->GetEntry()->Connect("ReturnPressed()", "MPXViewerControl",
							(MPXViewerControl *)this, "SetConfigAlgoValue()");

					cntrBoxes++;
				}
			}

		}

		//////////////////////////////////////////////////////
		/* Composite frame for algo config <Double_t> information, only if requested */
		if(!m_objs_conf_double.empty()) {

			for(int i = 0 ; i < (Int_t)m_objs_conf_double.size() ; i++) {

				// only if it belongs to this algorith
				if(((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetAuthor() == (*aItr).first)
				{
					label = ((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetConfigName();
					// create instance
					fconfBoxes[cntrBoxes] = new TextMargin(algoGroup[cntrAlgos], label);
					// build unique identifier
					TString mappingName = BuildKeyNameTextMargin((*aItr).first, label);
					// push into map for further reference
					m_confMapping[mappingName] = fconfBoxes[cntrBoxes];

					algoGroup[cntrAlgos]->AddFrame(fconfBoxes[cntrBoxes], new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));
					fconfBoxes[cntrBoxes]->GetTGNumberEntry()->Connect("ValueSet(Long_t)", "MPXViewerControl",
							this, "SetConfigAlgoValue()");
					fconfBoxes[cntrBoxes]->GetEntry()->Connect("ReturnPressed()", "MPXViewerControl",
							this, "SetConfigAlgoValue()");

					cntrBoxes++;
				}
			}

		}


		//////////////////////////////////////////////////////
		/* Composite frame for algo config <Int_t> information, only if requested */
		if(!m_objs_conf_int.empty()) {

			for(int i = 0 ; i < (Int_t)m_objs_conf_int.size() ; i++) {

				// only if it belongs to this algorithm
				if(((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetAuthor() == (*aItr).first)
				{

					label = ((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetConfigName();

					// create instance
					fconfBoxes[cntrBoxes] = new TextMargin(algoGroup[cntrAlgos], label);
					// build unique identifier
					TString mappingName = BuildKeyNameTextMargin((*aItr).first, label);
					// push into map for further reference
					m_confMapping[mappingName] = fconfBoxes[cntrBoxes];

					algoGroup[cntrAlgos]->AddFrame(fconfBoxes[cntrBoxes], new TGLayoutHints(kLHintsExpandX, 0, 0, 2, 2));
					fconfBoxes[cntrBoxes]->GetTGNumberEntry()->Connect("ValueSet(Long_t)", "MPXViewerControl",
							this, "SetConfigAlgoValue()");
					fconfBoxes[cntrBoxes]->GetEntry()->Connect("ReturnPressed()", "MPXViewerControl",
							this, "SetConfigAlgoValue()");

					cntrBoxes++;
				}
			}

		}

		contentsAlgo->AddFrame(algoGroup[cntrAlgos], new TGLayoutHints(kLHintsExpandX));

		cntrAlgos++;
	}

	// checkpoint
	if(nConf == cntrBoxes){
		Log << MSG::INFO << "All configuration values have been asigned a slot." << endreq;
	}


	//////////////////////////////////////////////////////////////////////////////
	// Meta data
	TGHorizontalFrame * contentsMetaList = new TGHorizontalFrame(contentsMeta);
	contentsMeta->AddFrame(contentsMetaList, new TGLayoutHints(kLHintsLeft | kLHintsTop |
			kLHintsExpandX | kLHintsExpandY, 0, 0, 5, 5));

	// Create the table and add it to the TableTest that is a TGMainFrame.

	fListBox = new TGListBox(contentsMetaList);
	fListBox->AddEntry("Frame Id (global): ", _frameId);
	fListBox->AddEntry("Frame Id: ", _frameIdLocal);
	fListBox->AddEntry("File : ", _fileName);
	fListBox->AddEntry("Acq Time: ", _acqTime);
	fListBox->AddEntry("N. Hits: ", _nHits);
	fListBox->AddEntry("N. Counts: ", _nCounts);
	fListBox->AddEntry("THL: ", _Thl);
	fListBox->AddEntry("Start time: ", _StS);
	fListBox->Resize(50,100);
	contentsMetaList->AddFrame(fListBox, new TGLayoutHints(kLHintsTop | kLHintsLeft |
			kLHintsExpandX | kLHintsExpandY,
			5, 5, 5, 5));

	///////////////////////////////////////////////
	// controls

	// Steering controls // put a TGHorizontalFrame inside contentsMeta
	TGVerticalFrame * contentsSteering = new TGVerticalFrame(contentsMeta);
	contentsMeta->AddFrame(contentsSteering, new TGLayoutHints(kLHintsBottom | kLHintsExpandY, 0, 0, 2, 2));

	///////////////////////////////////////////////////////////////////////////
	// FB
	TGHorizontalFrame * contentsSteeringFB = new TGHorizontalFrame(contentsSteering);
	contentsSteering->AddFrame(contentsSteeringFB, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));

	/* seek back button */
	fButtonSeekBack = new TGTextButton(contentsSteeringFB, "        <<       ", 450);
	fButtonSeekBack->Connect("Clicked()", "MPXViewerControl", this, "seekBack()");
	fButtonSeekBack->SetToolTipText("seek back");
	contentsSteeringFB->AddFrame(fButtonSeekBack, new TGLayoutHints(kLHintsCenterX, 10, 2, 2, 2));
	/* seek forward button */
	fButtonSeekForward = new TGTextButton(contentsSteeringFB, "      seek >>     ", 550);
	fButtonSeekForward->Connect("Clicked()", "MPXViewerControl", this, "seekForward()");
	fButtonSeekForward->SetToolTipText("seek forward (next interesting frame)");
	contentsSteeringFB->AddFrame(fButtonSeekForward, new TGLayoutHints(kLHintsCenterX, 2, 10, 2, 2));

	// Separator
	//TGHorizontal3DLine * separator2 = new TGHorizontal3DLine(contentsMeta);
	//contentsMeta->AddFrame(separator2, new TGLayoutHints(kLHintsRight | kLHintsExpandX));

	///////////////////////////////////////////////////////////////////////////
	// Jumps
	TGHorizontalFrame * contentsSteeringJump = new TGHorizontalFrame(contentsSteering);
	contentsSteering->AddFrame(contentsSteeringJump, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));

	// set Jump to
	fJumpToFrame = new TGNumberEntry(contentsSteeringJump, 0, 5, 600, TGNumberFormat::kNESInteger,
			TGNumberFormat::kNEANonNegative,
			TGNumberFormat::kNELLimitMinMax,
			0, 9999999);
	fJumpToFrame->Connect("ValueSet(Long_t)", "MPXViewerControl", this, "SetJumpTo()");
	fJumpToFrame->GetNumberEntry()->Connect("ReturnPressed()", "MPXViewerControl", this, "SetJumpTo()");
	contentsSteeringJump->AddFrame(fJumpToFrame, new TGLayoutHints(kLHintsCenterX, 0, 0, 2, 2));

	// Jump to button
	fButtonJumpToFrame = new TGTextButton(contentsSteeringJump, "  jump to  ", 650);
	fButtonJumpToFrame->Connect("Clicked()", "MPXViewerControl", this, "DoJumpTo()");
	fButtonJumpToFrame->SetToolTipText("jump to frame");
	contentsSteeringJump->AddFrame(fButtonJumpToFrame, new TGLayoutHints(kLHintsCenterX, 0, 0, 2, 2));

	///////////////////////////////////////////////////////////////////////////
	// min max Frame
	TGHorizontalFrame * contentsSteeringMinMax = new TGHorizontalFrame(contentsSteering);
	contentsSteering->AddFrame(contentsSteeringMinMax, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));

	// set mininum NumberEntry
	fHistoMin = new TGNumberEntry(contentsSteeringMinMax, 0, 5, 750, TGNumberFormat::kNESInteger,
			TGNumberFormat::kNEANonNegative,
			TGNumberFormat::kNELLimitMinMax,
			0, 99999);
	fHistoMin->Connect("ValueSet(Long_t)", "MPXViewerControl", this, "DoSetMin()");
	fHistoMin->GetNumberEntry()->Connect("ReturnPressed()", "MPXViewerControl", this, "DoSetMin()");
	contentsSteeringMinMax->AddFrame(fHistoMin, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0, 1, 2, 2));

	// set maximum NumberEntry
	fHistoMax = new TGNumberEntry(contentsSteeringMinMax, 0, 5, 800, TGNumberFormat::kNESInteger,
			TGNumberFormat::kNEANonNegative,
			TGNumberFormat::kNELLimitMinMax,
			0, 99999);
	fHistoMax->Connect("ValueSet(Long_t)", "MPXViewerControl", this, "DoSetMax()");
	fHistoMax->GetNumberEntry()->Connect("ReturnPressed()", "MPXViewerControl", this, "DoSetMax()");
	contentsSteeringMinMax->AddFrame(fHistoMax, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 1, 0, 2, 2));

	// set auto-adjust enable disable
	fdisableAutoAdjust = new TGCheckButton(contentsSteeringMinMax, "auto-adjust\nON/OFF");
	fdisableAutoAdjust->SetOn();
	fdisableAutoAdjust->Connect("Toggled(Bool_t)", "MPXViewerControl", this, "ToogleAutoAdjust()");
	contentsSteeringMinMax->AddFrame(fdisableAutoAdjust, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 2, 2));

	///////////////////////////////////////////////////////////////////////////
	// Switch Plot
	TGHorizontalFrame * contentsSwitchPlot = new TGHorizontalFrame(contentsSteering);
	contentsSteering->AddFrame(contentsSwitchPlot, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));

	// Plot switch: TOT, lvl1
	m_plotSwitchTOT = new TGRadioButton(contentsSwitchPlot, "TOT/counts");
	m_plotSwitchTOT->SetOn(true, true); // tot by default
	m_plotSwitchTOT->Connect("Toggled(Bool_t)", "MPXViewerControl", this, "TooglePlotValue_TOT()");
	contentsSwitchPlot->AddFrame(m_plotSwitchTOT, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 2, 2));

	// Plot switch: Energy through calib
	m_plotSwitchCalib = new TGRadioButton(contentsSwitchPlot, "Enegy(calib)");
	m_plotSwitchCalib->SetOn(false, false); // tot by default
	m_plotSwitchCalib->Connect("Toggled(Bool_t)", "MPXViewerControl", this, "TooglePlotValue_Calib()");
	contentsSwitchPlot->AddFrame(m_plotSwitchCalib, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 2, 2));

	m_plotSwitchLvl1 = new TGRadioButton(contentsSwitchPlot, "lvl1");
	m_plotSwitchLvl1->SetOn(false, false);
	m_plotSwitchLvl1->Connect("Toggled(Bool_t)", "MPXViewerControl", this, "TooglePlotValue_lvl1()");
	contentsSwitchPlot->AddFrame(m_plotSwitchLvl1, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 2, 2));

	m_plotSwitchMC = new TGRadioButton(contentsSwitchPlot, "MC_E");
	m_plotSwitchMC->SetOn(false, false);
	m_plotSwitchMC->Connect("Toggled(Bool_t)", "MPXViewerControl", this, "TooglePlotValue_MC()");
	contentsSwitchPlot->AddFrame(m_plotSwitchMC, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 2, 2));

	m_plotSwitchTruthMC = new TGRadioButton(contentsSwitchPlot, "MC_TruthE");
	m_plotSwitchTruthMC->SetOn(false, false);
	m_plotSwitchTruthMC->Connect("Toggled(Bool_t)", "MPXViewerControl", this, "TooglePlotValue_TMC()");
	contentsSwitchPlot->AddFrame(m_plotSwitchTruthMC, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 2, 2));

	/////////////////////////////////////////////////////////////////////////////////
	// Reprocess Frame
	TGHorizontalFrame * contentsSteeringRep = new TGHorizontalFrame(contentsSteering);
	contentsSteering->AddFrame(contentsSteeringRep, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));


	if(thereAreAlgoParameters)
	{
		// Reprocess
		fButtonReprocessEvent = new TGTextButton(contentsSteeringRep, "Reprocess", 850);
		fButtonReprocessEvent->Connect("Clicked()", "MPXViewerControl", this, "ReprocessCurrentFrame()");
		fButtonReprocessEvent->SetToolTipText("Reprocess current frame with new configuration values");
		contentsSteeringRep->AddFrame(fButtonReprocessEvent, new TGLayoutHints(kLHintsCenterX, 5, 5, 2, 2));

		// SaveConf
		fSaveConfiguration = new TGTextButton(contentsSteeringRep, "Save config", 900);
		fSaveConfiguration->Connect("Clicked()", "MPXViewerControl", this, "SaveConfiguration()");
		fSaveConfiguration->SetToolTipText("Save current list of configuration values to a plain text file");
		contentsSteeringRep->AddFrame(fSaveConfiguration, new TGLayoutHints(kLHintsCenterX, 5, 5, 2, 2));
	}

	// Quit button
	fButtonQuit = new TGTextButton(contentsSteeringRep, " Quit ", 950);
	fButtonQuit->Connect("Clicked()", "MPXViewerControl", this, "hardQuit()");
	fButtonQuit->SetToolTipText("Quit MAFalda");
	contentsSteeringRep->AddFrame(fButtonQuit, new TGLayoutHints(kLHintsCenterX, 5, 5, 2, 2));

	/////////////////////////////////////////////////////////////////////////////////
	// Special Frame
	TGHorizontalFrame * contentsSteeringSpecial = new TGHorizontalFrame(contentsSteering);
	contentsSteering->AddFrame(contentsSteeringSpecial, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));

	TGTextButton * fButtonDisconnect = new TGTextButton(contentsSteeringSpecial, "Disconnect", 850);
	fButtonDisconnect->Connect("Clicked()", "MPXViewerControl", this, "DisconnectAndRun()");
	fButtonDisconnect->SetToolTipText("Disconnect Viewer algorithm and complete run");
	contentsSteeringSpecial->AddFrame(fButtonDisconnect, new TGLayoutHints(kLHintsCenterX, 5, 5, 2, 2));

	/*
  ///////////////////////////////////////////////////////////////////////////////
  // credits
  TGHorizontalFrame * contentsCredits = new TGHorizontalFrame(contentsSteering);
  contentsSteering->AddFrame(contentsCredits, new TGLayoutHints(kLHintsRight | kLHintsExpandX, 0, 0, 0, 0));

  TGText * creditsJohn = new TGText(this, "John Idarraga <idarraga@cern.ch>");
  contentsCredits->AddFrame(creditsJohn, new TGLayoutHints(kLHintsRight, 5, 5, 2, 2));
	 */

	// Separator ////////
	TGHorizontal3DLine *separator2 = new TGHorizontal3DLine(contentsSteering);
	contentsSteering->AddFrame(separator2, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

	///////////////////////////////////////////////////////////////////////
	Connect("CloseWindow()", "TApplication", g_theApp, "Terminate()");
	DontCallClose();

	MapSubwindows();
	//Resize();

	SetWMSizeHints(GetDefaultWidth(), GetDefaultHeight(), 1000, 1000, 0 ,0);
	SetWindowName("MAFalda Viewer Control | J. Idarraga <idarraga@cern.ch>");
	MapRaised();

}

void MPXViewerControl::SaveConfiguration(){

	string fn = "configuration_algos.txt";
	ofstream outfile (fn.c_str());

	Log << MSG::INFO << "Saving configuration to \"" << fn.c_str() << "\"" << endreq;

	for(int i = 0 ; i < (Int_t)m_objs_conf_double.size() ; i++){
		outfile << ((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetAuthor()
	    				<< " " << ((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetConfigName()
	    				<< " " << ((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetConfigValue()
	    				<< " double" << std::endl;
	}

	for(int i = 0 ; i < (Int_t)m_objs_conf_float.size() ; i++){
		outfile << ((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetAuthor()
	    				<< " " << ((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetConfigName()
	    				<< " " << ((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetConfigValue()
	    				<< " float" << std::endl;
	}

	for(int i = 0 ; i < (Int_t)m_objs_conf_int.size() ; i++){
		outfile << ((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetAuthor()
	    				  << " " << ((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetConfigName()
	    				  << " " << ((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetConfigValue()
	    				  << " int" << std::endl;
	}

	for(int i = 0 ; i < (Int_t)m_objs_conf_bool.size() ; i++){
		outfile << ((ConfigurationValue<Bool_t> *)m_objs_conf_bool[i])->GetAuthor()
	    				  << " " << ((ConfigurationValue<Bool_t> *)m_objs_conf_bool[i])->GetConfigName()
	    				  << " " << ((ConfigurationValue<Bool_t> *)m_objs_conf_bool[i])->GetConfigValue()
	    				  << " bool" << std::endl;
	}

	outfile.close();

}

void MPXViewerControl::SetConfigAlgoValue(){

	// Each time SetConfigAlgoValueFloat is called I need to
	// sweep all slots, because there is no way I can find the id
	// of the one calling and match it to the id of m_objs_conf_XXX.
	// There are only a few, not a performance issue.

	// config values
	std::map<TString, Int_t>::iterator aItr = m_algorithmMap.begin();
	TString label = "";

	for( ; aItr != m_algorithmMap.end() ; aItr++) {

		for(int i = 0 ; i < (Int_t)m_objs_conf_int.size() ; i++) {

			// only if it belongs to this algorithm
			if(((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetAuthor() == (*aItr).first)
			{
				label = ((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetConfigName();

				// build unique identifier
				TString mappingName = BuildKeyNameTextMargin((*aItr).first, label);

				int vald = atoi(Form("%ld", m_confMapping[mappingName]->GetTGNumberEntry()->GetNumberEntry()->GetIntNumber()));
				((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->SetConfigValue(vald);

			}
		}


		for(int i = 0 ; i < (Float_t)m_objs_conf_float.size() ; i++) {

			// only if it belongs to this algorithm
			if(((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetAuthor() == (*aItr).first)
			{
				label = ((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetConfigName();

				// build unique identifier
				TString mappingName = BuildKeyNameTextMargin((*aItr).first, label);

				float val = atof(Form("%f", m_confMapping[mappingName]->GetTGNumberEntry()->GetNumberEntry()->GetNumber()));
				((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->SetConfigValue(val);

			}
		}


		for(int i = 0 ; i < (Double_t)m_objs_conf_double.size() ; i++) {

			// only if it belongs to this algorithm
			if(((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetAuthor() == (*aItr).first)
			{
				label = ((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetConfigName();

				// build unique identifier
				TString mappingName = BuildKeyNameTextMargin((*aItr).first, label);

				double val = atof(Form("%f", m_confMapping[mappingName]->GetTGNumberEntry()->GetNumberEntry()->GetNumber()));
				((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->SetConfigValue(val);

			}
		}

	}

}

void MPXViewerControl::SetJumpTo(){

	// Slot method connected to the ValueSet(Long_t) signal.
	feedBackInfoControl->jumpTo = (Long64_t)atoi(Form("%ld",fJumpToFrame->GetNumberEntry()->GetIntNumber()));

}
void MPXViewerControl::DoJumpTo(){

	feedBackInfoControl->jumpTo = (Long64_t)atoi(Form("%ld",fJumpToFrame->GetNumberEntry()->GetIntNumber()));
	g_theApp->Terminate();
	vSteerControl->direction = JUMP;
}

void MPXViewerControl::DoSetMin()
{
	// Slot method connected to the ValueSet(Long_t) signal.
	feedBackInfoControl->minHisto = atoi(Form("%ld",fHistoMin->GetNumberEntry()->GetIntNumber()));
	g_theApp->Terminate();
	vSteerControl->direction = DONT_MOVE;

}

void MPXViewerControl::DoSetMax()
{
	// Slot method connected to the ValueSet(Long_t) signal.
	feedBackInfoControl->maxHisto = atoi(Form("%ld",fHistoMax->GetNumberEntry()->GetIntNumber()));
	g_theApp->Terminate();
	vSteerControl->direction = DONT_MOVE;
}

void MPXViewerControl::ToogleAutoAdjust()
{

	// toogle autoAdjust flag
	if(feedBackInfoControl->autoAdjust){
		feedBackInfoControl->autoAdjust = kFALSE;
	}else{
		feedBackInfoControl->autoAdjust = kTRUE;
	}

	// Make sure the information has been taken from the NumberEntry
	// objects since no call to g_theApp->Terminate() is made here.  The
	// slots might not have been read yet (it happens when being in the
	// first frame).
	feedBackInfoControl->minHisto = atoi(Form("%ld",fHistoMin->GetNumberEntry()->GetIntNumber()));
	feedBackInfoControl->maxHisto = atoi(Form("%ld",fHistoMax->GetNumberEntry()->GetIntNumber()));

}

void MPXViewerControl::TooglePlotValue_TOT(){
	m_plotSwitchTOT->SetOn     (true, true);
	m_plotSwitchCalib->SetOn   (false, false);
	m_plotSwitchLvl1->SetOn    (false, false);
	m_plotSwitchMC->SetOn      (false, false);
	m_plotSwitchTruthMC->SetOn (false, false);
	feedBackInfoControl->SetPlotExclusive(__TOT);
}
void MPXViewerControl::TooglePlotValue_Calib(){
	m_plotSwitchTOT->SetOn     (false, false);
	m_plotSwitchCalib->SetOn     (true, true);
	m_plotSwitchLvl1->SetOn    (false, false);
	m_plotSwitchMC->SetOn      (false, false);
	m_plotSwitchTruthMC->SetOn (false, false);
	feedBackInfoControl->SetPlotExclusive(__Calib);
}

void MPXViewerControl::TooglePlotValue_lvl1(){
	m_plotSwitchTOT->SetOn     (false, false);
	m_plotSwitchCalib->SetOn   (false, false);
	m_plotSwitchLvl1->SetOn    (true, true);
	m_plotSwitchMC->SetOn      (false, false);
	m_plotSwitchTruthMC->SetOn (false, false);
	feedBackInfoControl->SetPlotExclusive(__LVL1);
}
void MPXViewerControl::TooglePlotValue_MC(){
	m_plotSwitchTOT->SetOn     (false, false);
	m_plotSwitchCalib->SetOn   (false, false);
	m_plotSwitchLvl1->SetOn    (false, false);
	m_plotSwitchMC->SetOn      (true, true);
	m_plotSwitchTruthMC->SetOn (false, false);
	feedBackInfoControl->SetPlotExclusive(__MC);
}
void MPXViewerControl::TooglePlotValue_TMC(){
	m_plotSwitchTOT->SetOn     (false, false);
	m_plotSwitchCalib->SetOn   (false, false);
	m_plotSwitchLvl1->SetOn    (false, false);
	m_plotSwitchMC->SetOn      (false, false);
	m_plotSwitchTruthMC->SetOn (true, true);
	feedBackInfoControl->SetPlotExclusive(__TRUTH_MC);
}


ControlFeedbackInfo * MPXViewerControl::ControlUpdate(ControlUpdateInfo * updateInfo_i, ControlFeedbackInfo * feedBackInfoControl_i){

	/* histos min max */
	fHistoMin->GetNumberEntry()->SetIntNumber(feedBackInfoControl_i->minHisto);
	fHistoMax->GetNumberEntry()->SetIntNumber(feedBackInfoControl_i->maxHisto);

	/* Meta data to Show, update on each frame */
	TGString frameIdS = "Global frame Id : ";
	frameIdS += updateInfo_i->frameId;
	frameIdS += " / ";
	frameIdS += updateInfo_i->nEntries;
	TGTextLBEntry * le = (TGTextLBEntry *)fListBox->GetEntry(_frameId);
	le->SetText(new TGString(frameIdS));
	//le->SetBackgroundColor(...);

	TGString frameIdLocalS = "Frame Id in file : ";
	frameIdLocalS += updateInfo_i->frameIdLocal;
	frameIdLocalS += " / ";
	frameIdLocalS += updateInfo_i->nEntriesLocal;
	le = (TGTextLBEntry *) fListBox->GetEntry(_frameIdLocal);
	le->SetText(new TGString(frameIdLocalS));

	TGString frameFileS = "File [";
	frameFileS += updateInfo_i->currentFileId;
	frameFileS += "] : ";
	frameFileS += updateInfo_i->currentFile;
	le = (TGTextLBEntry *) fListBox->GetEntry(_fileName);
	le->SetText(new TGString(frameFileS));

	TGString acqTimeS = "";
	acqTimeS += acqTimeS.Format("Acq. time: %.6f",updateInfo_i->acqTime); // N_SIGNIFICANT_FIGURES
	acqTimeS += " [s]";
	le = (TGTextLBEntry *)fListBox->GetEntry(_acqTime);
	le->SetText(new TGString(acqTimeS));

	double occupancy = 0.;
	if(m_matrixSizeX != 0 && m_matrixSizeY != 0) {
		occupancy = (double)updateInfo_i->nHitsInPad / (m_matrixSizeX * m_matrixSizeY);
	}
	TGString nHitsS = "N. Hits: ";
	nHitsS += updateInfo_i->nHitsInPad;
	if(occupancy > 0.01) { // show the occupancy only if it is very high > 1.%
		nHitsS += " | Occupancy: ";
		nHitsS += TString::Format("%1.f", occupancy).Data();
	}
	le = (TGTextLBEntry *)fListBox->GetEntry(_nHits);
	le->SetText(new TGString(nHitsS));

	TGString nCountsS = "N. Counts: ";
	if(updateInfo_i->nChargeInPad < 0){
		nCountsS += "overflow or not available !";
	}else{
		nCountsS += updateInfo_i->nChargeInPad;
	}
	le = (TGTextLBEntry *)fListBox->GetEntry(_nCounts);
	le->SetText(new TGString(nCountsS));

	TGString THLS = "THL: ";
	THLS += updateInfo_i->vTHL;
	le = (TGTextLBEntry *)fListBox->GetEntry(_Thl);
	le->SetText(new TGString(THLS));

	// THIS IS LOCAL TIME !!! WARNING !!!
	TGString StS = "Time(local): ";
	if(updateInfo_i->startTimeD == 0.){
		StS += "not available";
	}else{
		StS += CreateStartTimeString(updateInfo_i->startTimeD);
	}
	le = (TGTextLBEntry *)fListBox->GetEntry(_StS);
	le->SetText(new TGString(StS));

	// config values
	std::map<TString, Int_t>::iterator aItr = m_algorithmMap.begin();
	TString label = "";

	for( ; aItr != m_algorithmMap.end() ; aItr++) {

		for(int i = 0 ; i < (Int_t)m_objs_conf_int.size() ; i++) {

			// only if it belongs to this algorithm
			if(((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetAuthor() == (*aItr).first)
			{
				label = ((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetConfigName();
				// build unique identifier
				TString mappingName = BuildKeyNameTextMargin((*aItr).first, label);

				m_confMapping[mappingName]->
				GetTGNumberEntry()->
				SetNumber(
						((ConfigurationValue<Int_t> *)m_objs_conf_int[i])->GetConfigValue()
				);

			}
		}

		for(int i = 0 ; i < (Int_t)m_objs_conf_float.size() ; i++) {

			// only if it belongs to this algorithm
			if(((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetAuthor() == (*aItr).first)
			{
				label = ((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetConfigName();
				// build unique identifier
				TString mappingName = BuildKeyNameTextMargin((*aItr).first, label);

				m_confMapping[mappingName]->
				GetTGNumberEntry()->
				SetNumber(
						((ConfigurationValue<Float_t> *)m_objs_conf_float[i])->GetConfigValue()
				);

			}
		}


		for(int i = 0 ; i < (Int_t)m_objs_conf_double.size() ; i++) {

			// only if it belongs to this algorithm
			if(((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetAuthor() == (*aItr).first)
			{
				label = ((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetConfigName();
				// build unique identifier
				TString mappingName = BuildKeyNameTextMargin((*aItr).first, label);

				m_confMapping[mappingName]->
				GetTGNumberEntry()->
				SetNumber(
						((ConfigurationValue<Double_t> *)m_objs_conf_double[i])->GetConfigValue()
				);

			}
		}

	}

	MapSubwindows();
	Resize();
	MapWindow();

	/* Info that comes back to the Viewer */
	/* Set values in the TGNumberEntry instances as they come from MPXViewer */
	/* ---
     if(feedBackInfoControl_i->autoAdjust){

     fHistoMin->SetIntNumber((Int_t)feedBackInfoControl_i->minHisto);
     fHistoMax->SetIntNumber((Int_t)feedBackInfoControl_i->maxHisto);
     feedBackInfoControl->minHisto = feedBackInfoControl_i->minHisto;
     feedBackInfoControl->maxHisto = feedBackInfoControl_i->maxHisto;
     }
	 */

	return feedBackInfoControl;
}

TString MPXViewerControl::CreateStartTimeString(Double_t startTime){

	// THIS IS LOCAL TIME !!! WARNING !!!
	time_t localTime = (int)floor(startTime);

	return TString(ctime(&localTime));
}

void MPXViewerControl::PrintEventStats(){
	printf("In PrintEventStats()\n");
}

void MPXViewerControl::jump100F(){
	//g_theApp->Terminate();
	//g_direction = 51;
}
void MPXViewerControl::jump100B(){
	//g_theApp->Terminate();
	//g_direction = 50;
}
void MPXViewerControl::seekForward( ){

	g_theApp->Terminate();
	vSteerControl->direction = SEEK_FORWARD;
	//g_direction = 1;

}

void MPXViewerControl::seekBack(){

	g_theApp->Terminate();
	vSteerControl->direction = SEEK_BACKWARDS;
	//g_direction = 0;

}

void MPXViewerControl::softQuit(){

}

void MPXViewerControl::hardQuit(){

	g_theApp->Terminate();
	exit(1);
	//g_direction = 0;

}

void MPXViewerControl::ReprocessCurrentFrame(){

	// first pick up information in the configuration TGBoxes if any
	SetConfigAlgoValue();

	// now go ahead
	g_theApp->Terminate();
	vSteerControl->direction = DONT_MOVE;

}

void MPXViewerControl::DisconnectAndRun(){

	// Optional save config
	SaveConfiguration();

	// Through the feedBack MPXViewer will Schedule itself for disconnection
	feedBackInfoControl->m_disconnect = true;

	// now go ahead
	g_theApp->Terminate();
	vSteerControl->direction = SEEK_FORWARD;

}


void MPXViewerControl::SetHighlightersEnabled(Bool_t){

}

//______________________________________________________________________________
void MPXViewerControl::DoLeftMargin(char *)
{
	// Set left text margin.
	//fButton->SetLeftMargin(atoi(val));
	//gClient->NeedRedraw(fButton);
}  

void MPXViewerControl::ShowResults(){

	//Log << MSG::DEBUG << "Building the Control Bar" << endreq;
	Build();

}

void MPXViewerControl::GenerateAlgorithmMap(){

	// get a list of authors and the number of configuration values per author
	std::vector<CandidateContainer *>::iterator itr;

	itr = m_objs_conf_float.begin();
	for( ; itr != m_objs_conf_float.end() ; itr++){
		TString aut = ((ConfigurationValue<Float_t> *)(*itr))->GetAuthor();
		m_algorithmMap[aut]++;
	}

	itr = m_objs_conf_int.begin();
	for( ; itr != m_objs_conf_int.end() ; itr++){
		TString aut = ((ConfigurationValue<Int_t> *)(*itr))->GetAuthor();
		m_algorithmMap[aut]++;
	}

	itr = m_objs_conf_double.begin();
	for( ; itr != m_objs_conf_double.end() ; itr++){
		TString aut = ((ConfigurationValue<Double_t> *)(*itr))->GetAuthor();
		m_algorithmMap[aut]++;
	}

}

//////////////////////////////////////////////////////////////
// ControlFeedbackInfo implementation

Bool_t ControlFeedbackInfo::operator==(const ControlFeedbackInfo & oldFeedback) const {

	if(maxHisto == oldFeedback.maxHisto &&
			minHisto == oldFeedback.minHisto &&
			autoAdjust == oldFeedback.autoAdjust &&
			jumpTo == oldFeedback.jumpTo)
		return true;
	else
		return false;
}

Bool_t ControlFeedbackInfo::operator!=(const ControlFeedbackInfo & oldFeedback) const {
	return !(*this == oldFeedback);
}

#endif
