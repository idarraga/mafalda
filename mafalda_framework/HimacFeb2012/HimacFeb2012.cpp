/*
 * 	Copyright 2011 John Idarraga
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

#ifndef __HimacFeb2012_cpp
#define __HimacFeb2012_cpp

#include "HimacFeb2012.h"
#include "BlobsFinder/BlobsFinder.h"
#include "MAFTools/MAFTools.h"
#include "MAFTools/ConvexHull.h"
#include "MPXAlgo/Highlighter.h"

#include "TGraphDelaunay.h"
#include "TString.h"
#include "TCut.h"

#include "TGraph2D.h"
#include "TF2.h"

using namespace MSG;

ClassImp(HimacFeb2012)

// Gaus 2D
Double_t gaus2(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	Double_t r2 = Double_t((x[1]-par[3])/par[4]);
	return par[0]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}
Double_t fun2(Double_t *x, Double_t *par) {
	Double_t *p1 = &par[0];
	Double_t result = gaus2(x,p1);
	return result;
}


HimacFeb2012::HimacFeb2012() : MediPixAlgo(), CalibrationLoader(this) {

	// taken from data
	m_timepixClock = 0.;

	// used for basic selection
	m_sel.minInner = 10;
	m_sel.maxCircleOverElipse = 2.0;
	m_sel.borderExclusion = 1;
	m_sel.meshDiv = 0;

	// special selection
	m_sel.maxInner = 1000;
	m_sel.minClusterSize = 1;
	m_sel.maxClusterSize = 500;
	m_sel.minClusterTOT  = 1;
	m_sel.maxClusterTOT  = 5000;
	m_sel.minCircleOverElipse = 0.8; // Round selection, very tight
	m_sel.minHullCord = 0.8;
	m_sel.maxHullCord = 2.3;
	m_sel.bleeding_derivative_thl = 800.;
	m_sel.bleeding_derivative_counts = 10;
	m_sel.multipeak_divider = 3.;
	m_sel.bleed_min = 0.;
	m_sel.bleed_max = 100000.;
	m_sel.bleedToClusterE_cut = 0.7;
	m_sel.bleedConcentricDiffLimit = 0.05;
	m_sel.ndivitionsTGraph2D = 0;
	m_sel.distanceMultiplierClosingPoly = 2.0;

	m_acceptBorder = true;

	// An extra TCanvas *
	m_extraCanvas.clear();

	m_sizex = 256;
	m_sizey = 256;
	m_extraDrawing = 0;

	// Summatory of the elements in this queue at all times
	m_nContoursQueueSum = 0;
	m_nContoursQueueState = __nContour_reached_plateau;
	m_nContoursQueueState = __nContour_initState;

}

void HimacFeb2012::Init(){

	// to output ntuple
	getMyTree()->Branch("frameId", &m_frameId);
	getMyTree()->Branch("nClusters", &m_nClusters);
	getMyTree()->Branch("innerPixels", &m_innerPixels);
	getMyTree()->Branch("clusterSize", &m_clusterSize);
	getMyTree()->Branch("clusterTOT", &m_clusterTOT);
	getMyTree()->Branch("clusterEnergy", &m_clusterEnergy);
	getMyTree()->Branch("clusterEnergyReg", &m_clusterEnergyReg);

	getMyTree()->Branch("clusterGeoX", &m_clusterGeoX);
	getMyTree()->Branch("clusterGeoY", &m_clusterGeoY);
	getMyTree()->Branch("pixelTOT", &m_pixelTOT);
	getMyTree()->Branch("circleArea", &m_circleArea);
	getMyTree()->Branch("ellipseArea", &m_ellipseArea);
	getMyTree()->Branch("circleOverEllipse", &m_circleOverEllipse);
	getMyTree()->Branch("hullCordsFraction", &m_hullCordsFraction);

	getMyTree()->Branch("ebleed", &m_ebleed);
	getMyTree()->Branch("ebleed2", &m_ebleed2);
	getMyTree()->Branch("ehot", &m_ehot);
	getMyTree()->Branch("ehot2", &m_ehot2);

	getMyTree()->Branch("volumeRemainsLevel0", &m_volumeRemainsLevel0);
	getMyTree()->Branch("ebleedOverCluster", &m_ebleedOverCluster);

	getMyTree()->Branch("areaSkirtTotalRatioFactors", &m_areaSkirtTotalRatioFactors);
	getMyTree()->Branch("clusterLength", &m_clusterLength);

	// A configuration value that can be tuned from the Viewer
	//RegisterConfigurationValue(&m_minNClusters, "minNClusters");

	RegisterConfigurationValue(&m_sel.minInner, "minInner");
	RegisterConfigurationValue(&m_sel.maxInner, "maxInner");
	RegisterConfigurationValue(&m_sel.minClusterSize, "minClusterSize");
	RegisterConfigurationValue(&m_sel.maxClusterSize, "maxClusterSize");
	RegisterConfigurationValue(&m_sel.maxCircleOverElipse, "maxCircleOverElipse");
	RegisterConfigurationValue(&m_sel.bleeding_derivative_thl, "bleedDerTHL");
	RegisterConfigurationValue(&m_sel.bleeding_derivative_counts, "bleedDerCounts");
	RegisterConfigurationValue(&m_sel.multipeak_divider, "multipeakDivider");
	RegisterConfigurationValue(&m_sel.bleed_min, "bleed_min");
	RegisterConfigurationValue(&m_sel.bleed_max, "bleed_max");
	RegisterConfigurationValue(&m_sel.bleedConcentricDiffLimit, "bleedConcentricDiff");

	//RegisterConfigurationValue(&m_sel.minClusterTOT, "minClusterTOT");
	//RegisterConfigurationValue(&m_sel.maxClusterTOT, "maxClusterTOT");
	//RegisterConfigurationValue(&m_sel.minCircleOverElipse, "minCircleOverElipse");
	//RegisterConfigurationValue(&m_sel.maxCircleOverElipse, "maxCircleOverElipse");
	RegisterConfigurationValue(&m_sel.minHullCord, "minHullCord");
	RegisterConfigurationValue(&m_sel.maxHullCord, "maxHullCord");
	RegisterConfigurationValue(&m_sel.meshDiv, "meshDiv");
	RegisterConfigurationValue(&m_extraDrawing, "extraDrawing");
	RegisterConfigurationValue(&m_sel.bleedToClusterE_cut, "bleedToClusterE_cut");
	RegisterConfigurationValue(&m_sel.ndivitionsTGraph2D, "ndivitionsTGraph2D");
	RegisterConfigurationValue(&m_sel.distanceMultiplierClosingPoly, "distanceMultiplierClosingPoly");

	// FIXME  ... not handling correctly booleans
	//RegisterConfigurationValue(&m_acceptBorder, "acceptBorder");

	// Clustersize map
	m_clusterSizeMap = new TH2F("clusterSizeMap","clusterSizeMap", 256, 0, 256, 256, 0, 256);
	m_clusterSizeEntries = new TH2I("clusterSizeEntries","clusterSizeEntries", 256, 0, 256, 256, 0, 256);
	m_clusterTOTMap = new TH2F("clusterTOTMap", "clusterTOTMap", 256, 0, 256, 256, 0, 256);
	m_clusterEMap = new TH2F("clusterEMap", "clusterEMap", 256, 0, 256, 256, 0, 256);

	// for now this is only for G04 - Larry's 300um, no Al layer on top
	if(m_inputMask) {
		LoadMask();
	}
}

void HimacFeb2012::LoadMask(){

	int maxocc = m_sizex * m_sizey;
	m_Mask = new int[maxocc];
	for(int i = 0 ; i < maxocc ; i++){
		m_Mask[i] = 0;
	}

	// Detector G04
	// define the little square x: 164 --> 176, y: 60 --> 100
	for(int j = 0 ; j < m_sizey ; j++) {
		for(int i = 0 ; i < m_sizex ; i++) {
			if( (i >= 164 && i <= 176) && (j >= 60 && j <= 100) ) {
				m_Mask[ MAFTools::XYtoC(i,j,m_sizex) ] = 1;
			}
		}
	}

	/*
	int tempv;
	int cntr = 0;
	while( infile_extendedmask->good() && cntr < maxocc) {
	 *infile_extendedmask >> tempv;
		if(infile_extendedmask->eof()) break; // check EOF
		m_extendedMask[cntr] = tempv;
		cntr++;
	}
	 */

}

vector<pair<double, double> > HimacFeb2012::FindLimitPixels_MaxMinXY(blob theBlob) {

	// result vector
	vector<pair<double, double> > limits;
	limits.push_back( make_pair( 0 , 0) );
	limits.push_back( make_pair( 0 , 0) );

	set< pair<Int_t, Int_t> > contents = theBlob.GetContentSet();
	set< pair<Int_t, Int_t> >::iterator i = contents.begin();

	int minx = 0x0fffffff, maxx = 0, miny = 0x0fffffff, maxy = 0;

	for( ; i != contents.end() ; i++) {

		if((*i).first  < minx) minx = (*i).first;

		if((*i).second < miny) miny = (*i).second;

		if((*i).first  > maxx) maxx = (*i).first;

		if((*i).second > maxy) maxy = (*i).second;

	}

	limits[0] = make_pair( (double)minx, (double)miny );
	limits[1] = make_pair( (double)maxx, (double)maxy );

	return limits;
}

void HimacFeb2012::Execute(){

	// See the dacs
	//DumpDACs(MSG::INFO);

	// Figure out the timepix clock for this frame
	m_timepixClock = GetTimepixClock();
	// In the case of USBLite this will give zero and I have to use the other clock
	if( m_timepixClock == 0. ) m_timepixClock = GetMpxClock();

	// clean up stuff for the extra window display
	vector<TCanvas *>::iterator ican = m_extraCanvas.begin();
	for( ; ican != m_extraCanvas.end() ; ican++ ) {
		delete (*ican);
	}
	m_extraCanvas.clear();
	m_graph2DVector.clear();
	m_graphVector.clear();
	m_tf2Vector.clear();
	m_th2Vector.clear();

	// clean up all slices
	vector<vector<TGraph *> >::iterator ialls = m_allslices.begin();
	for( ; ialls != m_allslices.end() ; ialls++ ){
		(*ialls).clear();
	}
	vector<vector<TF1 *> >::iterator iallk = m_allkernels.begin();
	for( ; iallk != m_allkernels.end() ; iallk++ ){
		(*iallk).clear();
	}
	vector<vector<TGraph * > >::iterator iallq = m_clusterCountours.begin();
	for( ; iallq != m_clusterCountours.end() ; iallq++){
		(*iallq).clear();
	}
	m_clusterCountours.clear();

	// Ask the store gate if the previous algorithm (BlobsFinder --> reponsible for clustering)
	//  sent any objects to the StoreGate.
	Int_t lastObject = GetNumberOfObjectsWithAuthor("BlobsFinder");
	if(lastObject == 0)
		return;

	// If so, get the pointer to the last object.  BlobsFinder, as it comes out of the box, sends
	//  a single object containing all clusters.
	m_aB = (AllBlobsContainer *) GetObjectFromAuthor("BlobsFinder", lastObject-1);

	// AllBlobsContainer is a box full of Blobs(Clusters). Now we can iterate over all clusters inside.
	vector<blob> blobsVector = m_aB->GetBlobsVector();
	m_frameId = GetFrameId();
	m_nClusters = (int) blobsVector.size();
	Log << MSG::INFO << "Number of blobs from clustering = " << m_nClusters << endreq;
	vector<blob>::iterator blobsItr = blobsVector.begin(); //allBlobs.begin();

	// a pixel
	pair<int, int> pix;
	int clusterIndx = 0;
	// just to store a hot blob for display
	TGraph2D * g = 0x0;
	bool skipFlag = false;

	for ( ; blobsItr != blobsVector.end() ; blobsItr++ ) {

		// personal pronoun ;)
		blob b = *blobsItr;

		// Minimum selection
		if ( ! m_acceptBorder ) {
			// Do not accept clusters in the border
			if( ClusterInBorder(*blobsItr) ) continue;
		}

		if(b.bP.nInnerPixels <= m_sel.minInner) continue;

		/////////////////////////////////////////////////////////////////////////////////////////////
		// Selection !!!
		// see if any of the pixels is in the mask
		if(m_inputMask) {
			if( MAFTools::PixelInListIsAClusterConstituent(b, m_Mask, m_sizex, m_sizey) ) continue;
		}

		if( ! ParticleTagging(*blobsItr) ) continue;

		vector< pair<double, double> > limitPixels = FindLimitPixels_MaxMinXY(*blobsItr);
		double clusterLength = MAFTools::CalcDistance( limitPixels[0], limitPixels[1] );
		m_clusterLength.push_back( clusterLength );

		// start storing info
		//m_hullCordsFraction.push_back( hullfraction );

		// cluster Size
		m_clusterSize.push_back( (*blobsItr).bP.nPixels );

		// inner pixels
		m_innerPixels.push_back( (*blobsItr).bP.nInnerPixels );

		// Save total charge in TOT
		m_clusterTOT.push_back( (*blobsItr).bP.totalCharge );

		// Save a few geometrical properties
		m_circleArea.push_back( (*blobsItr).bP.circleArea );
		m_ellipseArea.push_back( (*blobsItr).bP.ellipseArea );
		m_circleOverEllipse.push_back( (*blobsItr).bP.circleArea / (*blobsItr).bP.ellipseArea );

		// Calibration //////////////////////////////////////////////////////////////////////////////////////////////
		// Pairs [(X,Y), counts]
		list< pair < pair<int, int>, int > >::iterator listItr;
		double clusterEnergy = 0.;
		int clusterTOT = 0;
		double cE = 0.;
		bool badCalibInCluster = false;
		double clockFactor = 1.;



		for ( listItr = (*blobsItr).blobContent.begin() ; listItr != (*blobsItr).blobContent.end() ; listItr++ ) {

			if(! IsMCData() ) {


				// pixel corrdinates
				pix = (*listItr).first;
				int tot = GetMatrixElement(pix);
				clusterTOT += tot;

				// clock
				clockFactor = m_timepixClock / CalibrationLoader::GetCalibClk();
				Log << MSG::LOOP_DEBUG << "Timepix clock = " << m_timepixClock
						<< " | calib clock " << CalibrationLoader::GetCalibClk()
				<< " | factor = " << clockFactor
				<< endreq;

				tot *= clockFactor; // clock factor

				if( CalibIsOK() ) {

					cE = CalculateAndGetCalibEnergy(pix, tot);

					if (cE < 0.) { // bad pixel or bad calib
						badCalibInCluster = true;
					}

				}

				// Store the calib energy for this pixel in the data model
				SetCalibEnergy(pix, cE);
				// And in the cluster
				b.SetPixelEnergy(cE, pix);

				clusterEnergy += cE;

				Log << MSG::LOOP_DEBUG << "pix(" << pix.first << " , " << pix.second << ") : "
						<< "TOT = " << tot << " | E(calib) = " << cE << endreq;

			} else {  // MC data

				pix = (*listItr).first;
				cE = GetMatrixElementMCEdep(pix); // Truth energy --> LET
				clusterEnergy += cE;

			}

		}

		// Set cluster energy
		b.SetClusterEnergy( clusterEnergy );

		// do not process this cluster
		if(badCalibInCluster) continue;

		/*
		// Storing the energy of each pixel as the energy of the cluster
		for ( listItr = (*blobsItr).blobContent.begin() ; listItr != (*blobsItr).blobContent.end() ; listItr++ ) {
			pix = (*listItr).first;
			SetCalibEnergy(pix, clusterEnergy);
		}
		 */

		// Some maps done only on selected objects
		m_clusterSizeMap->Fill( b.bP.geoCenter_x , b.bP.geoCenter_y, b.bP.clusterSize);
		m_clusterSizeEntries->Fill( b.bP.geoCenter_x , b.bP.geoCenter_y);

		m_clusterTOTMap->Fill( b.bP.geoCenter_x , b.bP.geoCenter_y, b.bP.clusterTOT );
		m_clusterEMap->Fill( b.bP.geoCenter_x , b.bP.geoCenter_y, clusterEnergy);

		// write to file
		//g2->Write(funcname);

		// Shape study
		double FittedVolumeFactor = 0.;
		//if(!g) {

		// Before anything let's have a finner grid
		TH2F * h = MAFTools::CropCluster(b, GetFrameId(), "E", "one");
		double meanEnergy = b.bP.clusterEnergy / (double)b.bP.clusterSize;
		m_th2Vector.push_back( h );

		// Get profiles for all slices in this cluster
		vector<TGraph * > slices = MAFTools::GetAllProfiles(b, GetFrameId(), "X", "E");
		Log << MSG::INFO << "Mean energy = " << meanEnergy << endreq;
		// Before sending them for drawing select only the slices containing points above the mean
		vector<TGraph * > slicesAboveMean = SelectSlicesAboveMean(slices, meanEnergy);
		m_allslices.push_back( slicesAboveMean );
		// Create kernel density functions for all this important slices
		// I need vector<double> for each histogram
		vector<TF1 * > kerneldensity = GetKernelDensityFunctions(clusterIndx, slicesAboveMean);
		m_allkernels.push_back( kerneldensity );
		// Find critial points
		//FindAllCritialPoints(slicesAboveMean, kerneldensity);

		/*
		// Fits to Gaus2D
		// Fit
		const int npar = 5;
		// Amplitude, mean, sigma, mean, sigma
		int hotest_E = MAFTools::GetEnergyOfHotestPixel(b);
		Double_t f2params[npar] = { hotest_E, b.bP.width_x/2. , 1, b.bP.width_y/2. , 1};
		TString f2_name = "f2_";
		f2_name += b.bP.geoCenter_x;
		f2_name += "_";
		f2_name += b.bP.geoCenter_y;
		TF2 * f2 = new TF2(f2_name, fun2, 0, b.bP.width_x, 0, b.bP.width_y, npar);
		f2->SetParameters(f2params);
		//f2->FixParameter(0, hotest_E);
		//TString fitPars = "RQN";
		TString fitPars = "RN";
		h->Fit(f2,fitPars,"");
		m_tf2Vector.push_back(f2);

		// See this volume
		double calc = 0.,vol = 0.;
		for (double j = 0 ; j < b.bP.width_y ; j++) {
			for (double i = 0 ; i < b.bP.width_x ; i++) {
				calc = f2->Eval( i + 0.5, j + 0.5 );
				if(calc > 1) {
					vol += f2->Eval(i,j);
					//Log << MSG::INFO << i << " , " << j << " = " << f2->Eval(i,j) << endreq;
				}
			}
		}

		Log << MSG::INFO << "ClusterEnergy before = " << TString::Format("%.2f", clusterEnergy) << " keV" << endreq;

		// TODO ! Changing the calibrated eneryg vlaue !!!!
		FittedVolumeFactor = vol/ b.bP.clusterTOT;
		// get a simple factor
		//clusterEnergy /= FittedVolumeFactor;

		Log << MSG::INFO << "Fitted volume " << TString::Format("%.1f",vol) <<
				" | fine integral = " << TString::Format("%.1f", f2->Integral(0, b.bP.width_x, 0, b.bP.width_y, 0.001)) << endreq;

		Log << MSG::INFO << "Factor to measured volume = " << TString::Format("%.2f",FittedVolumeFactor) << endreq;
		double pval = TMath::Prob( f2->GetChisquare()/f2->GetNDF(), f2->GetNDF() );

		Log << MSG::INFO << "P - value = " << TString::Format("%.2f", pval) << endreq;

		Log << MSG::INFO << "ClusterEnergy = " << TString::Format("%.2f", clusterEnergy) << " keV" << endreq;
		Log << MSG::INFO << "Geo      Center (" << b.bP.geoCenter_x << "," <<  b.bP.geoCenter_y << ")" << endreq;
		Log << MSG::INFO << "Weighted Center (" << b.bP.weightedCenter_x << "," <<  b.bP.weightedCenter_y << ")" << endreq;
		 */

		// Change the cluster information.  Add the energy

		// Store the cluster energy
		m_aB->SetClusterEnergy(clusterEnergy, clusterIndx);
		(*blobsItr).bP.clusterEnergy = clusterEnergy; // keV

		// Save total charge in keV
		m_clusterEnergy.push_back( clusterEnergy );
		// store the geo center of the cluster
		m_clusterGeoX.push_back( (*blobsItr).bP.geoCenter_x );
		m_clusterGeoY.push_back( (*blobsItr).bP.geoCenter_y );

		///////////////////////////////////////////////////////////////////////
		// Cluster bleeding ///////////////////////////////////////////////////
		// ClusterBleeding(b, f2, h);
		double bleeding = ClusterBleeding(b, clusterIndx, 0x0, h);



		if ( bleeding > m_sel.bleedToClusterE_cut ) {

			Highlighter * hl = new Highlighter(b.bP.weightedCenter_x,
					b.bP.weightedCenter_y,
					"arrow", this);
			hl->SetLineWidth(1);
			FillValuesForDisplay(hl, b);

			//PullToStoreGateAccess(hl, MPXDefs::DO_NOT_SERIALIZE_ME);

			// Skip flag
			skipFlag = true;

		}

		//if(clusterEnergy < 5000) {
		//	skipFlag = true;
		//}

		// a blob index integer counter
		clusterIndx++;

		// FIXME !
		break;

	}


	if( skipFlag ) {
		Signals * sig = new Signals(this, Signals::__SIGNAL_SKIP);
		//PullToStoreGateAccess(sig, MPXDefs::DO_NOT_SERIALIZE_ME);
	}

	if( GetHitsInPad() / double(256*256) < 0.05 ) {
		Signals * sig = new Signals(this, Signals::__SIGNAL_SKIP);
		//PullToStoreGateAccess(sig, MPXDefs::DO_NOT_SERIALIZE_ME);
	}

	/*
	// look at the contents of bleeding and decide if I want to skip
	vector<double>::iterator bleed_itr = m_ebleed.begin();
	bool skip = true;
	for ( ; bleed_itr != m_ebleed.end() ; bleed_itr++ ) {
		if( (*bleed_itr) > m_sel.bleed_min && (*bleed_itr) < m_sel.bleed_max ) {
			skip = false;
		}
	}
	if ( skip ) {
		Signals * sig = new Signals(this, Signals::__SIGNAL_SKIP);
		PullToStoreGateAccess(sig, MPXDefs::DO_NOT_SERIALIZE_ME);
	}
	 */

	// paint extra stuff
	if( m_extraDrawing ) {

		if(!m_tf2Vector.empty()) {
			m_extraCanvas.push_back( DrawInSeparateWindow(m_tf2Vector, MSG::DEBUG) );
		}

		if(!m_th2Vector.empty()) {
			m_extraCanvas.push_back( DrawInSeparateWindow(m_th2Vector, MSG::DEBUG) );
		}

		if( !m_graphVector.empty() ) {
			//m_extraCanvas.push_back( DrawInSeparateWindow(m_graphVector, m_ebleed, MSG::DEBUG) );
			m_extraCanvas.push_back( DrawInSeparateWindow(m_graphVector, m_ebleedOverCluster, MSG::DEBUG) );
		}

		if( !m_clusterCountours.empty() ) {
			m_extraCanvas.push_back( DrawInSeparateWindow(m_clusterCountours, "contour", MSG::DEBUG ) );
		}

		int cntrcn = 0;
		iallk = m_allkernels.begin();
		for(ialls = m_allslices.begin() ; ialls != m_allslices.end() ; ialls++ ) {
			TString cnamecn;
			cnamecn += cntrcn;
			//m_extraCanvas.push_back( DrawInSeparateWindow((*ialls), (*iallk), cnamecn, MSG::DEBUG ) );

			iallk++; // functions iterator
			cntrcn++;
		}

		if( ! m_graph2DVector.empty() ) {
			m_extraCanvas.push_back( DrawInSeparateWindow( m_graph2DVector, MSG::DEBUG, "tri1", "log" ) );
		}


	}

	// Fill the output tree of this algorithm
	getMyTree()->Fill();

	// WARNING ! don't forget to clean up your variables for the next TTree::Fill call
	m_nClusters = 0;
	m_innerPixels.clear();
	m_clusterSize.clear();
	m_clusterTOT.clear();
	m_clusterEnergy.clear();
	m_clusterEnergyReg.clear();

	m_clusterGeoX.clear();
	m_clusterGeoY.clear();

	m_pixelTOT.clear();

	m_circleArea.clear();
	m_ellipseArea.clear();
	m_circleOverEllipse.clear();

	m_hullCordsFraction.clear();

	m_ebleed.clear();
	m_ebleed2.clear();
	m_ehot.clear();
	m_ehot2.clear();

	m_volumeRemainsLevel0.clear();
	m_ebleedOverCluster.clear();
	m_volcanoBit.clear();
	m_clusterLength.clear();

	m_areaSkirtTotalRatioFactors.clear();

}

double HimacFeb2012::IsoVolumeAnalysisBleeding(cluster c, TH2 * h, double sigmas) {

	// Take the originally defined ellipse (from Clustering) and calculate the amount of energy out of that ellipse
	// Should be very little or zero
	double a = c.bP.ellipseA;
	double b = c.bP.ellipseB;
	double a_last = a, b_last = b;

	// weighted center
	int wc_tr_x = TMath::Floor( c.bP.width_x/2. + (c.bP.weightedCenter_x - c.bP.geoCenter_x) );
	int wc_tr_y = TMath::Floor( c.bP.width_y/2. + (c.bP.weightedCenter_y - c.bP.geoCenter_y) );
	Log << MSG::INFO << "Translated weighted center : (" << wc_tr_x << "," << wc_tr_y << ")" << endreq;

	//list< pair < pair<int, int>, int > > cdes = c.GetClusterDescription();
	//list< pair < pair<int, int>, int > >::iterator i = cdes.begin();
	int xbins = h->GetNbinsX();
	int ybins = h->GetNbinsY();
	double x, y, theta, point_r, r2;
	int n_inside = 0, e_inside = 0;
	double fraction_pixels, fraction_e, der;
	vector<double> x_vec;
	vector<double> der_vec;
	vector<double> a_vec;
	vector<double> b_vec;
	vector<double> fraction_e_vec;

	int der_cntr = 0;
	// First point in the queue will be full fraction 100%
	queue<double> gradient_q;
	gradient_q.push( 1. );

	while ( a > 0. && b > 0.) {

		// scorers
		n_inside = 0;
		e_inside = 0;
		for (int j = 1 ; j <= ybins ; j++) {
			for (int i = 1 ; i <= xbins ; i++) {

				// Only deal with hit pixels
				if( h->GetBinContent(i,j) == 0 ) continue;
				// Coordinates in the new coordinate system
				x = i-1;
				y = j-1;
				//Log << MSG::INFO << x << "," << y << " -- > ";
				// Apply translation which brings weighted center to the origin
				x -= wc_tr_x;
				y -= wc_tr_y;
				// Calculate how many of these point lie inside a given ellipse with
				// axis a and b as calculated by the clustering
				// Get theta
				if (x != 0.) theta = TMath::ATan( y / x );
				else theta = 0.;
				// Get r
				point_r = TMath::Sqrt( x*x + y*y );
				// Calculate the limit radius at a given Theta
				r2 = a*a*b*b;
				r2 /= (   ( b*b * TMath::Cos(theta) * TMath::Cos(theta) )
						+ ( a*a * TMath::Sin(theta) * TMath::Sin(theta) )  );
				// inside the ellipse condition
				if( point_r*point_r < r2 ) {
					n_inside++;
					e_inside += h->GetBinContent(i,j);
					//xv.push_back( x );
					//yv.push_back( y );
				}
				//Log << MSG::INFO << x << "," << y << " : " << h->GetBinContent(i,j) << endreq;
			}
		}
		// Calculate fractions
		fraction_pixels = (double)n_inside / (double)c.bP.clusterSize;
		fraction_e = (double)e_inside / c.bP.clusterEnergy;

		Log << MSG::DEBUG << "Ellipse : " << TString::Format("%.1f", a) << "," << TString::Format("%.1f", b)
		<< " --> Fraction pixels = " << fraction_pixels
		<< " | fraction E = " << fraction_e << endreq;

		// Calculate the derivative.  Pop the first item in the queue
		// Save the derivative to a vector
		if( gradient_q.size() == 5 ) {

			der = MAFTools::DerivativeFivePointsStencil(gradient_q);
			Log << MSG::DEBUG << "der = " << der << endreq;
			gradient_q.pop(); // first in first out

			// save the derivative to make a plot
			x_vec.push_back( der_cntr );
			der_vec.push_back( der );
			a_vec.push_back( a );
			b_vec.push_back( b );
			fraction_e_vec.push_back( fraction_e );
		}

		// Push to the queue
		gradient_q.push( fraction_e );

		// Increment various counters and values
		der_cntr++;
		a -= 0.5;
		b -= 0.5;

	}

	// Determine the limits here
	int min_der_size = 4;
	if(der_vec.size() < min_der_size) {
		Log << MSG::WARNING << "Can't determine the trend with less than 3 items in the derivative." << endreq;
		Log << MSG::WARNING << "This means the clusters are too small for this type of study." << endreq;
		return 0.;
	}

	int N = (int) der_vec.size();
	bool foundSkirtLimit = false;
	double threshold_fraction = 0.;
	int shift_back = min_der_size - 1;

	for(int der_i = 0 ; der_i < N ; der_i++) {

		Highlighter * hl = new Highlighter(
				c.bP.weightedCenter_x,
				c.bP.weightedCenter_y,
				a_vec[der_i], b_vec[der_i],
				this);
		hl->SetLineWidth(1);
		hl->SetLineColor(kWhite);
		//PullToStoreGateAccess(hl, MPXDefs::DO_NOT_SERIALIZE_ME);

		// Define skirt limit
		if( der_vec[der_i] >= m_sel.bleedConcentricDiffLimit) {

			foundSkirtLimit = true;
			Highlighter * hl = new Highlighter(
					c.bP.weightedCenter_x,
					c.bP.weightedCenter_y,
					a_vec[der_i - shift_back],
					b_vec[der_i - shift_back],
					this);
			hl->SetLineWidth(1);
			hl->SetLineColor(kWhite);
			//PullToStoreGateAccess(hl, MPXDefs::DO_NOT_SERIALIZE_ME);

			// Pick up the value two points before
			threshold_fraction = fraction_e_vec[der_i - shift_back];

		}

		if(foundSkirtLimit) break;

	}

	double bleed_frac = 0., bleed = 0.;
	if (foundSkirtLimit) {
		TGraph * g = new TGraph( x_vec.size(), &x_vec[0], &der_vec[0] );
		m_graphVector.push_back( g );
		bleed_frac = 1. - threshold_fraction;
		bleed = c.bP.clusterEnergy*bleed_frac;
		m_ebleedOverCluster.push_back( bleed_frac ); // The skirt fraction := bleed
		m_ebleed.push_back( bleed );
	}

	return bleed;
}

double HimacFeb2012::ConcentricBleeding(cluster c, TH2 * h, double sigmas) {

	// weighted center
	int wc_tr_x = TMath::Floor( c.bP.width_x/2. + (c.bP.weightedCenter_x - c.bP.geoCenter_x) );
	int wc_tr_y = TMath::Floor( c.bP.width_y/2. + (c.bP.weightedCenter_y - c.bP.geoCenter_y) );

	Log << MSG::INFO << "Translated weighted center : (" << wc_tr_x << "," << wc_tr_y << ")" << endreq;

	// Draw concentric circles pick up
	// Init the spiral
	steerSpiral spiral;
	MAFTools::RewindSpiral(spiral);

	// Start adding up making a spiral, keep a fifo of what's coming through and decide
	// where to stop.
	queue<double> gradient_q;
	// first point
	pair<int, int> pos = make_pair(wc_tr_x + 1, wc_tr_y + 1); // +1 ! i will be working in bins coordinates
	double content = h->GetBinContent(pos.first, pos.second);
	gradient_q.push( content );
	Log << MSG::INFO << "Point : " << pos.first << ", " << pos.second << " : " << content << endreq;

	bool lowgradient = false;
	double der = 0.;
	int cntr = 0, cntr_smoothness = 0;
	vector<double> x_vec;
	vector<double> der_vec;

	set<pair<int, int> > core_bins;

	while ( ! lowgradient ) {

		Log << MSG::DEBUG << "-----------------------------------------" << endreq;
		Log << MSG::DEBUG << "Queue size = " << gradient_q.size() << endreq;

		if( gradient_q.size() == 5 ) {
			der = MAFTools::DerivativeFivePointsStencil(gradient_q);
			Log << MSG::DEBUG << "der = " << der << endreq;
			gradient_q.pop(); // first in first out

			// save the derivative to make a plot
			x_vec.push_back( cntr );
			der_vec.push_back( der );

			// Stop if the relax condition have been met or
			// if more than 70% of the cluster has been studied already
			if( TMath::Abs(der) < m_sel.bleeding_derivative_thl ) {
				cntr_smoothness++;
			} else {
				cntr_smoothness = 0;
			}
			if(cntr_smoothness > m_sel.bleeding_derivative_counts) lowgradient = true;

			if( cntr > TMath::Floor(c.bP.clusterSize*0.7) ) lowgradient = true;

			// If still running in this loop this means that these pixels belong to the core
			// of the cluster.  Store the coordinates
			core_bins.insert( make_pair(pos.first, pos.second) );

		}

		pos = MAFTools::GetNextPosition(pos, spiral);
		content = h->GetBinContent(pos.first, pos.second);
		Log << MSG::DEBUG << "Point : " << pos.first << ", " << pos.second << " : " << content << endreq;
		gradient_q.push( content );

		cntr++;
	}

	TGraph * g = new TGraph( x_vec.size(), &x_vec[0], &der_vec[0] );
	m_graphVector.push_back( g );

	// Bleeding calculation.
	// Consider the whole cluster.  Take out the pixels belonging to the core of the cluster.
	int Nx = h->GetNbinsX();
	int Ny = h->GetNbinsY();
	double bleeding = 0.;
	for(int j = 1 ; j <= Ny ; j++) {
		for(int i = 1 ; i <= Nx ; i++) {
			if( core_bins.find( make_pair(i,j) ) == core_bins.end() ) { // all those pixels not in the core
				bleeding += h->GetBinContent(i,j);
			}
		}
	}

	Log << MSG::INFO << "Energy bleed = " << bleeding << " keV" << endreq;

	m_ebleed.push_back( bleeding );
	if ( c.bP.clusterEnergy > 0. ) {
		m_ebleedOverCluster.push_back( bleeding/c.bP.clusterEnergy );
	} else {
		m_ebleedOverCluster.push_back( 0. );
	}

	return bleeding;
}

vector< pair<double, double> > HimacFeb2012::GetCritialPoints (TF2 * f, TH2 * h, double sigmas) {

	vector< pair<double, double> > crit;
	crit.push_back( make_pair(0., 0.) ); // first pair is the maximum

	double Nx = h->GetNbinsX();
	double Ny = h->GetNbinsY();
	double step = 0.2;
	double val = 0., max = 0.;
	for(double iy = step ; iy < Ny ; iy += step) {
		for(double ix = step ; ix < Nx ; ix += step) {

			val = f->Eval(ix, iy);
			if(val > max) {
				max = val;
				crit[0].first = ix;
				crit[0].second = iy;
			}

		}
	}

	// Search for radius at a certain number of sigmas.
	// This will define an ellipse
	// Amplitude, mean, sigma, mean, sigma
	double sigmar1 = f->GetParameter(2);
	double sigmar2 = f->GetParameter(4);

	// Decide on the size of the ellipse using sigma
	sigmar1 *= sigmas;
	sigmar2 *= sigmas;

	Log << MSG::INFO << "Max in fitted function = " << max << endreq;
	Log << MSG::INFO << "             center at : " << crit[0].first << " , " << crit[0].second << endreq;
	Log << MSG::INFO << "           sigmar1 * N = " << sigmar1 << endreq;
	Log << MSG::INFO << "           sigmar2 * N = " << sigmar2 << " where N = " << sigmas << endreq;

	// The ellipse is defined by
	// ( ( x - crit[0].first )/( sigma*r1 ) )^2 + ( ( y - crit[0].second )/( sigma*r2 ) )^2  =  1
	// I will calculate the amount of charge out of that region
	// I will pick up two components.  Bleeding + delta rays
	// Do this operation on the 2D histogram here instead.  Watch the discrete space
	double xcoor, ycoor, ycalc;
	double bledding = 0.;
	for(int j = 1 ; j <= Ny ; j++) {

		for(int i = 1 ; i <= Nx ; i++) {

			// Use as coordinate the center of the bin to sample the continous function
			xcoor = (double)(i-1) + 0.5;
			ycoor = (double)(j-1) + 0.5;
			// Inside the ellipse condition
			ycalc = 1 - ( (xcoor - crit[0].first)*(xcoor - crit[0].first) / (sigmar1*sigmar1) );
			// ycalc = sigmar2*TMath::Sqrt(ycalc) + crit[0].second;
			// Out the ellipse
			if( ycalc < 0. ) {
				bledding += h->GetBinContent(i,j);
			}

		}

	}

	Log << MSG::INFO << "              Bleeding = " << bledding << " keV" << endreq;
	m_ebleed.push_back( bledding );


	return crit;
}


bool HimacFeb2012::NContourQueueIsFull() {

	if( m_nContoursQueue.size() == __max_in_nContours_queue ) return true;

	return false;
}

// This function keeps the queue up to date and calculated it's sum
//  which is always kept available in the member var m_nContoursQueueSum
bool HimacFeb2012::PushToNContourQueue ( int nConts ) {

	// Add to the summatory of the queue
	m_nContoursQueueSum += nConts;

	bool ret = false; // false if not full

	// If the queue is full then remove the oldest element and calculate the Sum of all elements inside
	if( m_nContoursQueue.size() == __max_in_nContours_queue ) {
		int frontval = m_nContoursQueue.front();       // Access the oldest value
		m_nContoursQueue.pop();                        // Get rid of it
		m_nContoursQueueSum -= frontval;			   // Remove it from the sum
		ret = true;									   // Queue is full
	}

	m_nContoursQueue.push( nConts );                   // Push the new val

	return ret;  // True if the queue it is full
}

double HimacFeb2012::IsoEnergyContoursAnalysisBleeding(cluster b, int clusterIndx){

	// I may modify the TOT values in the cluster.  Let's make a copy with such pattern
	cluster b_s;
	double weightslice = 0.;

	// Get the current cluster description
	list< pair < pair<Int_t, Int_t>, Int_t > > cl = b.GetClusterDescriptionCalibrated();

	// Square root (or something else) this description
	list< pair < pair<Int_t, Int_t>, Int_t > >::iterator cl_itr = cl.begin();
	int newtot = 0;
	for ( ; cl_itr != cl.end() ; cl_itr++ ) {
		// newtot = (int) TMath::Ceil( TMath::Log( (double) ((*cl_itr).second) ) );
		// newtot = (int) TMath::Ceil( TMath::Sqrt( (double) ((*cl_itr).second) ) );
		newtot = (*cl_itr).second;
		// if (newtot <= 0) newtot = 0;
		(*cl_itr).second = newtot;
	}

	// Load this description into the cloned cluster
	b_s.SetClusterDescription( cl );
	b_s.CalculateProperties(m_sizex, m_sizey, 1.0);

	// Get TGraph2D from this cluster
	TGraph2D * g = MAFTools::ConvertClusterToGraph2D(b_s, m_sel.ndivitionsTGraph2D, clusterIndx, "TOT", "test");
	m_graph2DVector.push_back( g );

	// FIXME !
	//m_graph2DVector.push_back( g );

	// This forces TGraphDelaunay to calculate the triangles ! IMPORTANT !
	g->GetXaxis();
	// Get the hotest pixel value.  A pixel can be noisy
	cl_itr = MAFTools::FindHotestPixel( &cl );
	double highest_val = (*cl_itr).second;
	cl_itr = MAFTools::FindColdestPixel( &cl );
	double coldest_val = (*cl_itr).second;
	if(coldest_val == 0) coldest_val = 1;
	// Define first a very fine step dedicated to the search of the first iso energy line
	double step = ( highest_val - coldest_val ) / highest_val;
	//double step = ( highest_val - coldest_val ) / 500.;
	Log << MSG::INFO << "Cluster [" << clusterIndx << "] Range : " << coldest_val << " --> " << highest_val << " | step = " << step << " keV"<< endreq;

	// Vector to hold all the contours
	vector<TGraph *> contours;
	// Get the convex hull.  Notice that the TGraph2D has zeros.  Let{s get rid of the zeroes first.
	TGraph * g_zs = MAFTools::GetTGraphZeroSuppression(g);
	vector<pair<int, int> > convh = MAFTools::ConvexHull(g_zs);
	// Convert the convex hull into a TGraph to calculate the area
	int itN = (int)convh.size();
	TGraph * gch = new TGraph( itN );
	for(int itc = 0 ; itc < itN ; itc++ ) { gch->SetPoint(itc, convh[itc].first, convh[itc].second); }
	double charea = b.bP.clusterSize;
	double clusterVolume = b.bP.clusterEnergy;
	// The area of the whole cluster in pixels squared is nothing else than the cluster Size
	// which is very similar here to the area of the convex hull := gch->Integral();
	// I use the convex hull here just to have a full Axis range.  I will schedule for drawing if requested.
	gch->SetLineColor(kWhite);
	gch->SetMarkerSize(0);
	contours.push_back( gch );

	// Now I can ask for the contours
	int Ncount = 0;
	int level_cntr = 0, valid_level_cntr = 0;
	bool firstContourFound = false;
	vector< vector<TGraph *> > allContoursOneCluster; // store all contours at all selected levels
	vector<double> allContoursEnergyLevels;

	//vector<double> internalVolumeAtLevel;

	for (double cont = coldest_val ; cont < highest_val ; cont += step) {

		TList * listContour = g->GetContourList(cont);
		if(listContour) {
			Ncount = listContour->GetEntries();
		} else {
			Ncount = 0;
		}
		Log << MSG::INFO << "N countours = " << Ncount << endreq;
		if ( Ncount == 0 ) break; // If no more contours, end the loop.

		// If there is more than 1 contour available many things can be happening
		// 1) It is due to volcano effect.  I need to check if the secondary contour is inside the first contour.
		//    To do so I calculate the Hull of the outter contour, and then the hull of the joint set of points from
		//    the first and second contour.  If the hull is the same the second contour is interior to the first one
		//    and this is the finger print of volcano.  In this case I use only the first contour for the running
		//    derivative and keep the second one to study the volcano.

		// Save all polygons at this level
		vector<TGraph * > pol_one_level;
		int index_first_contour = -1;
		double remainsLevel0Saved = false;

		for ( int i = 0 ; i < Ncount ; i++ ) {

			TString gname = "cluster_";
			gname += clusterIndx;
			gname += "_cont_";
			gname += valid_level_cntr;   // Contour level
			gname += "_";                  // If more than one contour at that level
			gname += i;
			TGraph * g_1d = (TGraph *) listContour->At( i );
			g_1d->SetName(gname);
			double polygonarea = g_1d->Integral(); // Calcuate area of this polygon
			pol_one_level.push_back( g_1d );

			// See if this can be the first contour
			// Compare the area of all of the contours here with the are of the convex hull of this cluster

			if ( ! firstContourFound ) {

				Log << MSG::DEBUG << "polygonarea / ch area = " << polygonarea / charea << " | level = " << cont << endreq;

				// If this level is reached, here's the location of the first enveloping iso energy line
				// Use a different condition here.  See if the polygon is closed.
				if( PolygonIsClosed(g_1d, polygonarea/charea) ) {
					//if( polygonarea / charea > m_sel.areaFractionPolyFull ) {

					Log << MSG::INFO << "Found first enveloping contour. Polygon area / full area = " << polygonarea / charea << " | level = " << cont << endreq;
					firstContourFound = true;
					index_first_contour = i;
					// Change the step to a coarser value
					double newstep = TMath::Ceil( TMath::Log(highest_val - coldest_val) );
					newstep = ( highest_val - coldest_val ) / (newstep*10.);
					Log << MSG::INFO << "Change step : " << step << " --> " << newstep << endreq;
					step = newstep;

					// Here's the first polygon.  Calculate the volume out of that polygon
					//Log << MSG::INFO << "Internal volume fraction = " << internal_volume/clusterVolume << endreq;
					//Log << MSG::INFO << "External volume fraction = " << external_volume/clusterVolume << endreq;


				}

			}

		}
		// Done with all contours at this level

		// Check if the first contour was found and save all contours at this level
		if ( firstContourFound ) {

			// Save all contours at this level
			allContoursOneCluster.push_back( pol_one_level );
			// and the energy level
			allContoursEnergyLevels.push_back( cont );

		}

		// clear polygons in this level
		pol_one_level.clear();

		//Log << MSG::INFO << "Level area = " << allpolygons_inthislevel_area << endreq;
		level_cntr++; 								// counter on all contours
		if( firstContourFound ) valid_level_cntr++;   // counter on good usable contours

	}

	// And I will use this one for calculations : allContoursOneCluster
	// This on is filled once per level and refers only to one cluster

	int nLevels = (int) allContoursOneCluster.size();
	Log << MSG::INFO << "Number of levels = " << allContoursOneCluster.size() << endreq;

	vector<double> volumePerLevel;
	vector<double> outVolumePerLevel;
	vector<double> areaPerLevel;
	vector<double> energyLevelPerLevel;

	double prevVolume = b.bP.clusterEnergy; // start with the full volume
	queue<int> nContoursDrag;
	int selectedLevelsCounter = 0;
	int selectedLevelsCrownIndex = -1;
	double selectedLevelsCrownEnergy = 0.;

	for ( int lev_i = 0 ; lev_i < nLevels ; lev_i++ ) {

		vector<TGraph *> levconts = allContoursOneCluster[lev_i];
		double energyAtLevel = allContoursEnergyLevels[lev_i];

		int nConts = (int) levconts.size();
		Log << MSG::INFO << "   [" << lev_i << "] : contours : " << nConts << endreq;

		// The levels such follow the following trend in terms of contours per level
		// 5 4 3 2 1 1 1 1 1 1 1 1 1 (1 for a long time) and then break again into multiple contours 1 1 1 1 2 3 3 3 3 etc.
		// I want to determine at which level this happens
		// The following code will keep a queue with 5 entries and checks if it is full as well.
		// The summatory is available at m_nContoursQueueSum
		if ( PushToNContourQueue( nConts ) ) { // It is full

			if ( m_nContoursQueueSum == __steady_sum_in_nContours_queue && m_nContoursQueueState == __nContour_initState ) { // Got to the stable plateau

				m_nContoursQueueState = __nContour_reached_plateau;

			} else if ( m_nContoursQueueSum > __excess_sum_in_nContours_queue && m_nContoursQueueState == __nContour_reached_plateau ) { // Left plateau and reached crown

				m_nContoursQueueState = __nContour_reached_crown;
				m_nContoursCrownLevel = lev_i - 1;  // This is the crown level
				// One time !
				TGraph2D * g_stripped = ExtractParOfTGraph2D( g, allContoursOneCluster[m_nContoursCrownLevel][0] );
				m_graph2DVector.push_back( g_stripped );

			}
		}

		// If the Queue reached the crown here then the previous level is the location of the tip of the mountain which is cut !

		// I need to identify if in the case of multiple contours per level
		// Do the combinatorics here.  For every polygon check with its partners if
		// there are polygons inside.  The insider polygon won't add to the calculated volume
		vector<int> badIndxs;
		for ( int cont_i = 0 ; cont_i < nConts - 1 ; cont_i++ ) {

			for ( int cont_j = cont_i + 1 ; cont_j < nConts; cont_j++ ) {

				// If the polygon is inside then   levconts[cont_j]   won't be taken into account for
				// the volume calculation and it'll be associated to a special list.
				bool inside = PolygonIsInside(levconts[cont_i], levconts[cont_j]);
				if ( inside ) badIndxs.push_back( cont_j );

				//Log << MSG::DEBUG << "            " << cont_i << " , " << cont_j << " | inside = " << inside << endreq;

			}

		}

		// Calculate the volume at this level
		double levelVolume = 0.;
		double levelOutVolume = 0.;
		double levelArea = 0.;

		for ( int cont_i = 0 ; cont_i < nConts; cont_i++ ) {

			// See if the polygon is to be taken into account for volume calculation, i.e. not in the badIndxs list
			if ( std::find(badIndxs.begin() , badIndxs.end(), cont_i) == badIndxs.end() ) {

				// Volumes
				double internal_volume = 0.;
				double external_volume = 0.;
				GetInternalExternalVolume(levconts[cont_i], g, internal_volume, external_volume);
				levelVolume += internal_volume;

				// Area
				levelArea += levconts[cont_i]->Integral();

			}		

		}
		// Outter volume for this level
		levelOutVolume += b.bP.clusterEnergy - levelVolume;

		// Schedule all of the contours at this level for drawing
		// if the volume changed
		if(levelVolume != prevVolume) {

			// Change the name of the level
			//TString gname = "cluster_";
			//gname += clusterIndx;
			//gname += "_cont_";
			//gname += valid_level_cntr;   // Contour level
			//gname += "_";                // If more than one contour at that level
			//gname += i;
			//g_1d->SetName(gname);

			for ( int cont_i = 0 ; cont_i < nConts; cont_i++ ) contours.push_back( levconts[cont_i] );

			Log << MSG::INFO << "Volume at level " << lev_i << " = " << levelVolume << " |  energy level = " << energyAtLevel << endreq;
			volumePerLevel.push_back( levelVolume );
			outVolumePerLevel.push_back( levelOutVolume );
			areaPerLevel.push_back( levelArea );
			energyLevelPerLevel.push_back( energyAtLevel );

			//////////////////////////////////////////////////////////////////////////////////////////////////////
			// hot !
			if ( m_nContoursQueueState == __nContour_reached_crown ) {
				selectedLevelsCrownIndex = selectedLevelsCounter;
				selectedLevelsCrownEnergy = energyAtLevel;
				// m_ehot.push_back( prevVolume );
			}

			// Actual levels selected
			selectedLevelsCounter++;
		}
		prevVolume = levelVolume;

		//for( int cont_i = 0 ; cont_i < levconts.size() )
	}

	if ( m_nContoursQueueState == __nContour_initState ) { // Deal with the cases where we couldn't even reach a plateau.  Probably a cut cluster.
		return 0;
	} else if ( m_nContoursQueueState != __nContour_reached_crown ) { // If crown not found because the cluster is regular.  Pick the last one.
		selectedLevelsCrownIndex = selectedLevelsCounter - 1;
		selectedLevelsCrownEnergy = energyLevelPerLevel[selectedLevelsCrownIndex];
		Log << MSG::INFO << "Crown NOT found.  The crown level will be the last contour.  Level index = " << selectedLevelsCounter << endreq;
	} else {
		Log << MSG::INFO << "Crown level at " << m_nContoursCrownLevel << " --> " << selectedLevelsCrownIndex << " in selection." << endreq;
	}
	// This one is used for drawing
	// Filled once per cluster with an array of contours.
	m_clusterCountours.push_back( contours );

	double skirtVolume = 0.;
	if ( (int)volumePerLevel.size() > 2 ) {

		// Calculate the differences now
		vector<double> volumeDifferences;
		int Nlevels = (int) volumePerLevel.size();
		int Ndiffs = Nlevels + 1;
		// First level.  Outer first level

		skirtVolume = b.bP.clusterEnergy - volumePerLevel[0];

		m_ebleed.push_back( b.bP.clusterEnergy - volumePerLevel[0] );
		m_ebleed2.push_back( b.bP.clusterEnergy - volumePerLevel[1] );
		m_ehot.push_back( volumePerLevel[selectedLevelsCrownIndex] );
		//m_ehot2.push_back( volumePerLevel[m_nContoursCrownLevel] );

		// Calculate differences.  This vector is one entry bigger than the vector containing the volumes per level.
		volumeDifferences.push_back( b.bP.clusterEnergy - volumePerLevel[0]);
		for ( int idiff = 1 ; idiff < Ndiffs ; idiff++ ) {
			volumeDifferences.push_back( volumePerLevel[idiff-1] - volumePerLevel[idiff] );
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Apply correction

		// Calculate the first few levels as a function of the hot energy compared to the level at which the crown is found !
		//int first_few_levels = (int) TMath::Floor( clusterVolume / volumePerLevel[selectedLevelsCrownIndex] );
		//int first_few_levels = volumePerLevel[selectedLevelsCrownIndex]/selectedLevelsCrownEnergy;
		int first_few_levels = (int) TMath::Ceil ( b.bP.clusterSize / (b.bP.clusterSize - areaPerLevel[0]) );
		m_areaSkirtTotalRatioFactors.push_back( b.bP.clusterSize / (b.bP.clusterSize - areaPerLevel[0]) ); // output

		Log << MSG::INFO << "- Crown volume = " <<  volumePerLevel[selectedLevelsCrownIndex]
                         << " | Crown level = " << selectedLevelsCrownEnergy
                         << " | Number of levels = " << first_few_levels << endreq;
		Log << MSG::INFO << "- Skirt volume = " <<  skirtVolume
						 << " | Skirt level = " << areaPerLevel[0] << endreq;
		Log << MSG::INFO << "- Cluster size = " << b.bP.clusterSize
					     << " | hot size = " << areaPerLevel[selectedLevelsCrownIndex]
					     << " | cold size =  " << areaPerLevel[0] << endreq;
		double acc = 0.;
		//double fact = ( 1. - (skirtVolume/highest_val) );
		double fact = ( 1. - (skirtVolume / energyLevelPerLevel[selectedLevelsCrownIndex]) );
		Log << MSG::INFO << "Factor = " << fact << endreq;
		if( m_nContoursQueueState == __nContour_reached_crown ) { // Volcano.  Explicit of mild.
			for ( int ilev = 0 ; ilev < Nlevels ; ilev++ ) {
				//if (ilev < first_few_levels) {
				if( ilev < 3 ) {
					acc += volumePerLevel[ilev] * fact;
					Log << MSG::INFO << "Level volume = " << volumePerLevel[ilev] << " | acc = " << acc << endreq;
				}
			}
		} else { // Saturation effect
			acc += outVolumePerLevel[selectedLevelsCrownIndex];
			//acc += volumePerLevel[selectedLevelsCrownIndex] * fact;
		}


		// Dump differences
		//for ( int idiff = 0 ; idiff < Ndiffs ; idiff++ ) {
		//	Log << MSG::INFO << "Diff = " << volumeDifferences[idiff] << endreq;
		//}

		double reg_let = 0.;
		// Good physics volume
		// reg_let = b.bP.clusterEnergy - volumePerLevel[selectedLevelsCrownIndex]; // Extract the hot crown"cylinder"
		// reg_let = reg_let * ( b.bP.clusterEnergy / (skirtVolume*20) );
		reg_let = acc;
		// volumePerLevel[selectedLevelsCrownIndex] +
		// (b.bP.clusterEnergy/skirtVolume);

		Log << MSG::INFO << "Skirt = " << skirtVolume << endreq;
		Log << MSG::INFO << "Regularized L.E.T. = " << reg_let << endreq;
		m_clusterEnergyReg.push_back( reg_let );
	}


	// Cleaning up
	m_nContoursQueueSum = 0;
	// Rewind queue
	while ( ! m_nContoursQueue.empty() ) m_nContoursQueue.pop();
	m_nContoursQueueState = __nContour_initState;
	m_nContoursCrownLevel = -1;

	return 0.;
}

double HimacFeb2012::HeightMissingFactor(){

	return __missing_height;
}

TGraph2D * HimacFeb2012::ExtractParOfTGraph2D ( TGraph2D * g, TGraph * g_extract) {

	TGraph2D * g2D = (TGraph2D *)g->Clone("temp_stripp_temp");

	int N2D = g2D->GetN();
	double * x2D = g2D->GetX();
	double * y2D = g2D->GetY();
	double * z2D = g2D->GetZ();

	int N = g_extract->GetN();
	double * x = g_extract->GetX();
	double * y = g_extract->GetY();

	for ( int i = 0 ; i < N2D ; i++ ) {

		if( g_extract->IsInside( x2D[i], y2D[i] ) % 2 == 0 ) { // Pair --> Outside
			z2D[i] = 0;
		}
	}

	g2D->Delete();

	return new TGraph2D(N2D, x2D, y2D, z2D);
}

bool HimacFeb2012::PolygonIsClosed(TGraph * poly, double areafraction) {

	// DumpPolygon(poly);

	// Look at the distance between the first and last point.  If it is way too different from the average
	//  distance in the entire set of points then this polygon is probably not closed.
	// The points come in order in "poly".

	// Calc average distance
	int N = poly->GetN();
	double * xx = poly->GetX();
	double * yy = poly->GetY();

	double averageDist = 0.;
	for ( int i = 0 ; i < N - 1 ; i++ ) {
		if( xx[i] != xx[i] ) return false; // Catch nan.  Can come from the contour algo.
		averageDist += MAFTools::CalcDistance(xx[i], yy[i], xx[i+1], yy[i+1]);
	}
	if( xx[N-1] != xx[N-1] ) return false; // Catch nan.  Can come from the contour algo.

	averageDist /= N-1;
	double distanceExtremes = MAFTools::CalcDistance(xx[0], yy[0], xx[N-1], yy[N-1]);

	Log << MSG::INFO << "Average distance = " << averageDist << " | extremes = " << distanceExtremes << endreq;

	// If the distance between the extremes is bigger than twice the average distance. Those two points are too far
	//  and we say the polygon is not closed.  But the contour can be very small.  For that reason I OR a second
	//  condition where i request the are of the polygon to be at least 20% of the total area of the cluster.
	if ( distanceExtremes > m_sel.distanceMultiplierClosingPoly*averageDist || areafraction < 0.2 ) return false;

	return true;
}

void HimacFeb2012::DumpPolygon(TGraph * poly){

	int N = poly->GetN();
	double * xx = poly->GetX();
	double * yy = poly->GetY();

	Log << MSG::DEBUG << "---- Dump polygon ----------------------------" << endreq;
	for ( int i = 0 ; i < N ; i++ ) {
		Log << MSG::DEBUG << "[" << i << "](" << xx[i] << "," << yy[i] << ")  ";
	}
	Log << endreq << MSG::DEBUG << "----------------------------------------------" << endreq;

}

// Will check if polyA is inside polyB
bool HimacFeb2012::PolygonIsInside ( TGraph * polyA_is_inside, TGraph * polyB ) {

	//int NA = polyA_is_inside->GetN();
	//double * xxA = polyA_is_inside->GetX();
	//double * yyA = polyA_is_inside->GetX();

	int NB = polyB->GetN();
	double * xxB = polyB->GetX();
	double * yyB = polyB->GetX();

	//for ( int iA = 0 ; iA < NA ; iA++ ) {

	int nPointsInside = 0;
	for ( int iB = 0 ; iB < NB ; iB++ ) {

		if( polyA_is_inside->IsInside( xxB[iB], yyB[iB] ) % 2 != 0 ) { // Odd
			nPointsInside++;
		}

	}

	// If all points in B are inside A
	if(nPointsInside == NB) return true;

	//}
	return false;
}

void HimacFeb2012::GetInternalExternalVolume(TGraph * contour, TGraph2D * cluster2D, double & int_volume, double & ext_volume) {

	// "cluster" is the TGraph2D representation of the cluster which shares coordinates
	// system with "contour".
	int N = cluster2D->GetN();
	double * xx = cluster2D->GetX();
	double * yy = cluster2D->GetY();
	double * zz = cluster2D->GetZ();

	int x, y, e;
	for (int i = 0 ; i < N ; i++) {

		x = xx[i];
		y = yy[i];
		e = zz[i];

		// Taken from
		//______________________________________________________________________________
		// Int_t TGraph::IsInside(Double_t x, Double_t y) const
		// Algorithm:
		// The loop is executed with the end-point coordinates of a line segment
		// (X1,Y1)-(X2,Y2) and the Y-coordinate of a horizontal line.
		// The counter inter is incremented if the line (X1,Y1)-(X2,Y2) intersects
		// the horizontal line. In this case XINT is set to the X-coordinate of the
		// intersection point. If inter is an odd number, then the point x,y is within
		// the polygon.
		if( contour->IsInside( (double)x, (double)y ) % 2 != 0 ) { // Odd
			int_volume += e;
		} else {
			ext_volume += e;
		}

	}

}

double HimacFeb2012::ClusterBleeding (cluster cl, int clusterIndex, TF2 * f, TH2 * h) {

	// Identify the maximum of the fit
	// This doesn't work correctly with delta rays
	//vector< pair<double, double> > crit = GetCritialPoints(f, h, 3.0);

	// Second technique to identify the bleeding.
	//return ConcentricBleeding(cl, h, 0.0);

	// Third technique based on concentric circles or iso-TOT curves
	//return IsoVolumeAnalysisBleeding(cl, h, 0.0);

	return IsoEnergyContoursAnalysisBleeding(cl, clusterIndex);

	/*
	list< pair < pair<Int_t, Int_t>, Int_t > > cdes = cl.GetClusterDescriptionCalibrated();
	list< pair < pair<int, int>, int > >::iterator itr = MAFTools::FindHotestPixel(&cdes);
	pair<int, int> hotpix = (*itr).first;
	int ehot = (*itr).second;
	list< pair < pair<int, int>, int > >::iterator itrcold = MAFTools::FindColdestPixel(&cdes, ehot);
	int ecold = (*itrcold).second;
	// estimator of bleeding level
	int bleeding_e = ehot/1000;
	// In case there is nothing belong the estimation
	// if( bleeding_e <= ecold ) bleeding_e = ecold*2;

	Log << MSG::INFO << "Hottest pixel (" << hotpix.first << " , " << hotpix.second << ") : "
			<< ehot << " | " << bleeding_e << endreq;

	int e_bleed = 0;
	for ( itr = cdes.begin() ; itr != cdes.end() ; itr++ ) {
		//Log << MSG::INFO << (*itr).first.first << " , " << (*itr).first.second << " : " << (*itr).second << endreq;
		if ( (*itr).second < bleeding_e ) {
			e_bleed += (*itr).second;
		}
	}

	Log << MSG::INFO << "Energy bleed = " << e_bleed << " keV" << endreq;
	//m_ebleed.push_back( e_bleed );

	if(e_bleed < 10.0) {
		Highlighter * hl = new Highlighter(cl.bP.geoCenter_x,
				cl.bP.geoCenter_y,
				"circle", this);
		hl->SetLineWidth(1);
		FillValuesForDisplay(hl, cl);

		PullToStoreGateAccess(hl, MPXDefs::DO_NOT_SERIALIZE_ME);
	}
	 */

}


vector<TF1 * > HimacFeb2012::GetKernelDensityFunctions(int indx, vector<TGraph *> gs) {

	vector<TF1 * > kernels;
	vector<TGraph *>::iterator i = gs.begin();

	vector<double> histo;
	for( ; i != gs.end() ; i++) {

		int N = (*i)->GetN();
		double * y = (*i)->GetY();
		for(int j = 0 ; j < N ; j++) {

			histo.push_back( y[j] );
			//Log << MSG::INFO << " --> " << y[j] << endreq;
		}

		TF1 * f = MAFTools::CreateKernelDensityFunction( indx, histo, 0.4 );
		kernels.push_back( f );

		// clear for next histo
		histo.clear();

	}

	return kernels;

}

void HimacFeb2012::FindAllCritialPoints(vector<TGraph * > g, vector<TF1* > f ) {

	vector<double> min;
	vector<double> max;
	vector<double>::iterator itr;

	vector<TGraph * >::iterator ig = g.begin();
	vector<TF1 * >::iterator ifunc = f.begin();

	int critpoints;
	for( ; ig != g.end() ; ig++ ) {

		critpoints = MAFTools::GetCriticalPoints((*ifunc), (*ig)->GetN(), min, max);

		Log << MSG::INFO << "Number of critial points [" << (*ig)->GetName() << "] = " << critpoints << endreq;
		Log << MSG::INFO << "Minimums at : ";
		for(itr = min.begin() ; itr != min.end() ; itr++) {
			Log << MSG::INFO << (*itr) << ", ";
		}
		Log << MSG::INFO << endreq;

		Log << "Maximums at : ";
		for(itr = max.begin() ; itr != max.end() ; itr++) {
			Log << MSG::INFO << (*itr) << ", ";
		}
		Log << MSG::INFO << endreq;

		// clear the vectors
		min.clear();
		max.clear();

		// second iterator
		ifunc++;
	}

}

vector<TGraph *> HimacFeb2012::SelectSlicesAboveMean(vector<TGraph *> slices, double mean) {

	vector<TGraph * > selection;

	vector<TGraph *>::iterator i = slices.begin();
	bool store = false;
	for( ; i != slices.end() ; i++) {

		store = false;
		int N = (*i)->GetN();
		double * y = (*i)->GetY();

		for(int j = 0 ; j < N ; j++) {
			if(y[j] > mean) {
				store = true;
			}
		}

		if(store) selection.push_back( *i );

	}

	return selection;
}

double HimacFeb2012::GetBiggestPerpDistanceSpecialCaseSlopeZero_Y(double ycut, vector<pair<int,int> > h){

	double dist = 0., dte = 0.;
	vector<pair<int,int> >::iterator itr = h.begin();

	for( ; itr != h.end() ; itr++){
		dte = TMath::Abs( ycut - (*itr).second ) ; // Just use the Y coordinate in this case
		if(dte > dist) dist = dte;
	}

	return dist;
}

double HimacFeb2012::GetBiggestPerpDistanceSpecialCaseSlopeZero_X(double xcut, vector<pair<int,int> > h){

	double dist = 0., dte = 0.;
	vector<pair<int,int> >::iterator itr = h.begin();

	for( ; itr != h.end() ; itr++){
		dte = TMath::Abs( xcut - (*itr).first ) ; // Just use the Y coordinate in this case
		if(dte > dist) dist = dte;
	}

	return dist;


}

double HimacFeb2012::GetBiggestPerpDistance(double slope, double cut, vector<pair<int,int> > h){

	double dist = 0., dte = 0.;
	vector<pair<int,int> >::iterator itr = h.begin();

	for( ; itr != h.end() ; itr++){
		dte = MAFTools::CalcPerpDistanceToLine(slope, cut, *itr);
		if(dte > dist) dist = dte;
	}

	return dist;
}

double HimacFeb2012::GetSmallestPerpDistance(double slope, double cut, vector<pair<int,int> > h){

	double dist = 1024., dte = 0.;
	vector<pair<int,int> >::iterator itr = h.begin();

	for( ; itr != h.end() ; itr++){
		dte = MAFTools::CalcPerpDistanceToLine(slope, cut, *itr);
		if(dte < dist) dist = dte;
	}

	return dist;
}

void HimacFeb2012::Finalize() {

	if( CalibIsOK() ) {

		MAFTools::TimePixCalibrationHandler * calib = GetCalibrationHandler();

		// store a surrogate function for monitoring
		calib->GetSurrogateTF1( make_pair(GetMaxWidthCalib()/2, GetMaxHeightCalib()/2) )->Write();

		calib->GetSurrogateTF1(make_pair(221, 174))->Write();
		calib->GetSurrogateTF1(make_pair(221, 175))->Write();
		calib->GetSurrogateTF1(make_pair(234, 211))->Write();

		calib->GetSurrogateTF1(make_pair(223 , 200))->Write();
		calib->GetSurrogateTF1(make_pair(224 , 199))->Write();
		calib->GetSurrogateTF1(make_pair(224 , 200))->Write();
		calib->GetSurrogateTF1(make_pair(224 , 201))->Write();
		calib->GetSurrogateTF1(make_pair(225 , 200))->Write();
		calib->GetSurrogateTF1(make_pair(225 , 201))->Write();

		calib->GetSurrogateTF1(make_pair(224 , 67))->Write();
		calib->GetSurrogateTF1(make_pair(111 , 70))->Write();

		calib->GetSurrogateTF1(make_pair(183 , 135))->Write();
		calib->GetSurrogateTF1(make_pair(198 , 137))->Write();

		// get a calib map
		Log << MSG::INFO << "Generating a calibration map" << endreq;
		int xDim = GetMatrixXdim();
		int yDim = GetMatrixYdim();
		TH2F * h_calib = new TH2F("calibMap","calib",xDim, 0, xDim-1,
				yDim, 0, yDim-1);
		pair<int, int> pix;
		for (int colItr = 0; colItr < xDim-1 ; colItr++) {
			for (int rowItr = 0; rowItr < yDim-1 ; rowItr++) {
				pix.first = colItr;
				pix.second = rowItr;
				// get a calibration for a TOT of 100
				h_calib->Fill(pix.first, pix.second, calib->GetE(pix, 100) );
			}
		}
		h_calib->Write();

	}

	// write the cuts
	//TString cutS = BuildCutGeo1();
	//TString cutS = BuildCutGeo_Hull();
	TString cutS = "";
	TCut cut(cutS);
	cut.Write( GetDataSetNumber()+"_cut" );

	// Normalize ClusterSize map
	int sx = m_clusterSizeMap->GetNbinsX();
	int sy = m_clusterSizeMap->GetNbinsY();
	int val, entries;
	double vald = 0.;

	for(int i = 1 ; i < sx ; i++) {
		for(int j = 1 ; j < sy ; j++) {

			// cluster size
			val = m_clusterSizeMap->GetBinContent(i, j);
			entries = m_clusterSizeEntries->GetBinContent(i, j);
			if(entries > 0) m_clusterSizeMap->SetBinContent(i, j, val/entries);

			// cluster TOT
			val = m_clusterTOTMap->GetBinContent(i, j);
			if(entries > 0) m_clusterTOTMap->SetBinContent(i, j, val/entries);

			// cluster E
			vald = m_clusterEMap->GetBinContent(i, j);
			if(entries > 0) m_clusterEMap->SetBinContent(i, j, vald/entries);

		}
	}

	m_clusterSizeMap->Write();
	m_clusterSizeEntries->Write();
	m_clusterTOTMap->Write();
	m_clusterEMap->Write();

}

TString HimacFeb2012::BuildCutGeo1(){
	TString cutS = "( innerPixels >= ";
	cutS += m_sel.minInner;
	cutS += " && innerPixels <= ";
	cutS += m_sel.maxInner;
	cutS += " )";

	cutS += " && "; // second part

	cutS += "( clusterSize >= ";
	cutS += m_sel.minClusterSize;
	cutS += " && clusterSize <= ";
	cutS += m_sel.maxClusterSize;
	cutS += " )";

	cutS += " && "; // third part

	cutS += "( clusterTOT >= ";
	cutS += m_sel.minClusterTOT;
	cutS += " && clusterTOT <= ";
	cutS += m_sel.maxClusterTOT;
	cutS += " )";
	return cutS;
}

TString HimacFeb2012::BuildCutGeo_Hull(){

	TString cutS = "( innerPixels >= ";
	cutS += m_sel.minInner;
	cutS += " && innerPixels <= ";
	cutS += m_sel.maxInner;
	cutS += " )";

	cutS += " && "; // second part

	cutS += "( clusterSize >= ";
	cutS += m_sel.minClusterSize;
	cutS += " && clusterSize <= ";
	cutS += m_sel.maxClusterSize;
	cutS += " )";

	cutS += " && "; // third part

	cutS += "( hullCordsFraction >=";
	cutS += TString::Format("%.3f", m_sel.minHullCord);
	cutS += " && hullCordsFraction <=";
	cutS +=TString::Format("%.3f", m_sel.maxHullCord);
	cutS += ")";

	return cutS;
}

int HimacFeb2012::CheckDiscontinuityWithModifiedMatrix(blob bl, map<int,int> totmap, int discontinuityTolerance){

	// New object containing all the blobs
	AllBlobsContainer * reclusteringContainer = new AllBlobsContainer((MediPixAlgo *) this);
	Log << MSG::DEBUG << "Re-clustering container --> " << reclusteringContainer << endreq;

	MAFTools::DigitalBlobReclustering(discontinuityTolerance, reclusteringContainer, bl, GetMatrixWidth(), GetMatrixHeight(), totmap);

	Log << MSG::DEBUG << "Re-clustering ... found " << reclusteringContainer->GetBlobsVector().size() << " clusters in the blob" << endreq;

	int nSubClusters = reclusteringContainer->GetBlobsVector().size();
	delete reclusteringContainer;

	return nSubClusters;
}

bool HimacFeb2012::ParticleTagging (blob b) {

	// Find the pixels below the average threshold, i.e. cold pixels
	double meanTOT = (double) b.bP.clusterTOT / (double) b.bP.clusterSize;
	double hottest = MAFTools::GetTOTOfHotestPixel( b );
	// cut-off half way between the mean TOT and the max
	double cutoff_cold = ( meanTOT + hottest ) / m_sel.multipeak_divider;
	set<pair<int, int> > coldpixels = MAFTools::FindPixelsBelow(b, cutoff_cold);
	//set<pair<int, int> >::iterator hitr = coldpixels.begin();

	// Make a copy of the cluster but masking the cold pixels
	cluster b2(b, coldpixels);
	//b2.AllDump();

	// Re-cluster and see if the hot spot is fragmented --> If that is the case this is probably an overlap
	map<int, int> totmap = GetTOTMap();
	MAFTools::MaskInMap(totmap, coldpixels, GetWidth());
	int subclusters = CheckDiscontinuityWithModifiedMatrix(b2, totmap, 0);
	Log << MSG::DEBUG << "number of sub-clusters : " << subclusters << endreq;

	// Check if the cluster touches the border
	bool touchesBorder = MAFTools::ClusterTouchesBorder(b, GetWidth(), GetHeight(), m_sel.borderExclusion);

	double circleOverEllipse = b.bP.circleArea / b.bP.ellipseArea;

	if (
			b.bP.nInnerPixels >= m_sel.minInner
			&&
			subclusters == 1   // Only one hot peak found in the cluster
			&&
			circleOverEllipse <= m_sel.maxCircleOverElipse // Roundness
			&&
			! touchesBorder    // Border selection
			&&
			b.bP.clusterSize <= m_sel.maxClusterSize
	) {

		// Extra parameters to display
		//int nextras = 1;
		//vector<TString> names(1, "hullr");
		//vector<double> values(1, hullr);

		Highlighter * hl = new Highlighter(b.bP.weightedCenter_x, b.bP.weightedCenter_y,
				b.bP.ellipseA, b.bP.ellipseB,
				this);
		hl->SetLineWidth(1);
		FillValuesForDisplay(hl, b);
		PullToStoreGateAccess(hl, MPXDefs::DO_NOT_SERIALIZE_ME);


		return true;
	}

	return false;
}

bool HimacFeb2012::ClusterInBorder(blob b){

	// Pairs [(X,Y), counts]
	list< pair < pair<int, int>, int > >::iterator listItr;

	pair<int, int> pix;
	for ( listItr = b.blobContent.begin() ; listItr != b.blobContent.end() ; listItr++ ) {
		pix = (*listItr).first;

		if ( pix.first <= 0 || pix.first >= GetMatrixWidth()-1 ) {
			return true;
		}
		if ( pix.second <= 0 || pix.second >= GetMatrixHeight()-1 ) {
			return true;
		}
	}

	return false;
}

void HimacFeb2012::GetMaxMin(vector<pair<int,int> > h, int * min_x, int * max_x, int * min_y, int * max_y) {

	vector<pair<int,int> >::iterator i = h.begin();
	*min_x = 1024;
	*max_x = 0;
	*min_y = 1024;
	*max_y = 0;
	int x, y;
	for( ; i != h.end() ; i++) {
		x = (*i).first;
		y = (*i).second;
		if( x < *min_x ) *min_x = x;
		if( x > *max_x ) *max_x = x;
		if( y < *min_y ) *min_y = y;
		if( y > *max_y ) *max_y = y;
	}

}

double HimacFeb2012::GetSmallestDistanceTo(pair<int, int> pix, vector<pair<int,int> > h){
	vector<pair<int,int> >::iterator i = h.begin();
	double dist = 1024.0, d = 0.;
	for( ; i != h.end() ; i++) {
		d = MAFTools::CalcDistance(*i, pix);
		if(d < dist) dist = d;
	}
	return dist;
}

double HimacFeb2012::GetBiggestDistanceTo(pair<int, int> pix, vector<pair<int,int> > h){

	vector<pair<int,int> >::iterator i = h.begin();
	double dist = 0.0, d = 0.;
	for( ; i != h.end() ; i++) {
		d = MAFTools::CalcDistance(*i, pix);
		if(d > dist) dist = d;
	}
	return dist;
}


void HimacFeb2012::ProcessHIMAC_MC(){



}

double paraboloid_f(double * xx, double * par){

	double a = par[0];
	double b = par[1];
	double c = par[2];
	double x = xx[0];
	double y = xx[1];

	return c*((x*x/a) + (y*y/b));
}

#endif


// Approach to calculate the area of a polygon based on Delaunay triangles.  TGraph::Integral turns out to do
// the trick much faster and realiably.
/*
		//if ( contour_cntr == 8 && i == 0 ) Log << MSG::INFO << "Integral = " << g_1d->Integral() << endreq;

		g_1d->SetName(gname);
		g_1d->SetLineWidth(3);
		g_1d->SetMarkerStyle(20);
		// push it back
		if ( contour_cntr == 8 && i == 0 ) contours.push_back( g_1d );

		int N = g_1d->GetN();
		double * x = g_1d->GetX();
		double * y = g_1d->GetY();
		double * z = new double[N];
		for(int jj = 0 ; jj < N ; jj++) z[jj] = 0.;
		TGraph2D * gdel = new TGraph2D(N, x, y, z);
		TGraphDelaunay * del = new TGraphDelaunay( gdel );
		del->SetMaxIter(10000);
		del->ComputeZ(x[0], y[0]); // Force the Delaunay data structure to get prepared
		del->FindAllTriangles();
		int deltriangs = del->GetNdt();
		if(contour_cntr == 8 && i == 0) Log << MSG::INFO << "Delaunay triangles = " << deltriangs << endreq;
		int * p = del->GetPTried();
		int * n = del->GetNTried();
		int * m = del->GetMTried();

		double triangle_area = 0., polygon_area = 0.;
		for ( int kk = 0 ; kk < deltriangs ; kk++ ) {

			// p, n, m are the vertices
			//Log << MSG::INFO
				//	<< "(" << x[ p[kk] ] << ", " << y[ p[kk] ] << "), "
				//	<< "(" << x[ n[kk] ] << ", " << y[ n[kk] ] << "), "
				//	<< "(" << x[ m[kk] ] << ", " << y[ m[kk] ] << ")"
				//	<< endreq;

			double xx[4] =  { x[ p[kk] ], x[ n[kk] ], x[ m[kk] ] , x[ p[kk] ] }; // repeat the first point so the triangle closes
			double yy[4] =  { y[ p[kk] ], y[ n[kk] ], y[ m[kk] ] , y[ p[kk] ] };

			//double xx[4] =  { 1, 7, 2, 1 }; // repeat the first point so the triangle closes
			//double yy[4] =  { 1, 1, 4, 1 };

			TGraph * g_del = new TGraph(4, xx, yy);
			if(contour_cntr == 8 && i == 0) contours.push_back( g_del );

			// Calculate the area of this triangle
			triangle_area = MAFTools::TriangleArea( xx, yy );
			polygon_area += triangle_area;
			if(contour_cntr == 8 && i == 0) Log << MSG::INFO << "Triangle area = " << triangle_area << endreq;

		}
		if(contour_cntr == 8 && i == 0) Log << MSG::INFO << "Polygon area = " << polygon_area << endreq;
 */
