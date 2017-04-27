/**
 *  Author: John Idarraga <idarraga@cern.ch>
 *  TOT calibration for a Timepix device
 */

#define TOTCalib_cxx
#include "TOTCalib.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TLegend.h>

#include <time.h>

TOTCalib::TOTCalib(){
	m_TOTMultiplierFactor = 1.0;
}

void TOTCalib::Loop()
{
	//   In a ROOT session, you can do:
	//      Root > .L TOTCalib.C
	//      Root > TOTCalib t
	//      Root > t.GetEntry(12); // Fill t data members with entry number 12
	//      Root > t.Show();       // Show values of entry 12
	//      Root > t.Show(16);     // Read and show values of entry 16
	//      Root > t.Loop();       // Loop on all entries
	//

	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	if( m_nFrames > 0 && m_nFrames < nentries ) nentries = m_nFrames;

	Long64_t nbytes = 0, nb = 0;

	// these lines open blank canvas
	//TCanvas * c2 = new TCanvas("histo" + TString( m_calhandler->GetSourcename()) );
	//c2->cd();

	// prepare output ROOT
	// m_output_root = new TFile("output_calib.root", "RECREATE");

	double percentage = 0;
	cout << "Retrieving single hit information from the ROOT file ..." << endl;
	for (Long64_t jentry = 0; jentry < nentries ; jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		if(jentry%1000 == 0) { // refresh progress bar every 1000 frames
			percentage = ( (double)jentry / (double)nentries ) * 100;
			printProgBar( (int) percentage );
		}

		vector<int>::iterator coorItr = m_SingleHitCoor->begin();
		vector<int>::iterator totItr = m_SingleHitTOT->begin();
		pair<int, int> pix;
		int tot;

		for ( ; coorItr != m_SingleHitCoor->end() ; ) {

			//pix = XtoXY( *coorItr, __matrix_width );
			//m_calibhistos[ XYtoX( pix, __matrix_width ) ][ tot ]++;

			tot = *totItr;
			tot = tot * m_TOTMultiplierFactor; // should be 1.0 otherwise changed by the user

			if(tot < m_nbins) {
				m_calibhistos[ *coorItr ][ tot ]++;
			} else {
				if(m_verbose == __VER_DEBUG) cout << " [WARNING] Pixel contains tot off the requested number of bins : " << *coorItr << " | tot = " << tot << endl;
			}

			coorItr++; totItr++;
		}


		// if (Cut(ientry) < 0) continue;
	}
	printProgBar( (int) 100 );
	cout << endl; // extra endl to finish correctly the progress bar

	/////////////////////////////////////////////////////////////////////////////
	// Create kernel density functions
	// Plus peak finding

	cout << "Creating kernel functions ..." << endl;
	int barstep = ( m_maxpix - m_minpix ) / 100;
	if ( barstep == 0 ) barstep = 100;

	for (int i = m_minpix ; i <= m_maxpix ; i++) {

		// Skip the bad pixels
		if ( PixelInBadPixelList(i) ) continue;

		// See first of all if this particular pixel has valid data.
		// It can be a masked pixel or something noisy
		if( VectorSum( m_calibhistos[i] ) == 0. ) continue;

		if ( i % barstep == 0 ) {
			percentage = ( (double)(i - m_minpix) / (double)(m_maxpix - m_minpix) ) * 100;
			printProgBar( (int) percentage );
		}
		//cout << "On function : " << i << " | size = " << m_calibhistos[i].size() << endl;

		// Create kernel density function and sample it (using a stencil)
		//  in order to find the critical points
		vector<double> min, max;
		GetCriticalPoints(i, min, max);

		// Fill the maps with the critical points.  The
		//  maximums should correspond to the matching
		//  values in energy.  We use the handler to
		//  get those energies.  This establishes the
		//  relation : TOT(E)
		m_critPointsMax[i] = max;
		m_critPointsMin[i] = min;

	}
	printProgBar( (int) 100 );
	cout << endl; // extra endl to finish correctly the progress bar

}

double TOTCalib::VectorSum(vector<double> v){

	vector<double>::iterator i = v.begin();
	double sum = 0.;
	for( ; i != v.end() ; i++ ) {
		sum += *i;
	}

	return sum;
}

enum {
	__s_pos = 0,
	__s_neg,
};

void TOTCalib::FillHisto(int x, int tot) {

	m_calibhistos[x][tot]++;

}
/**
 * Simple linear fit using the points in p
 */
void TOTCalib::GetLinearFit(double & a, double & b, vector<pair<double,double> > p){

	int N = (int) p.size();
	
	if (N != 0){
        
        TGraph * g = new TGraph(N);
        if(m_verbose == __VER_DEBUG) cout << "Linear fit using : ";
        for( int i = 0 ; i < N ; i++ ) {
            g->SetPoint( i, p[i].first, p[i].second );
            if(m_verbose == __VER_DEBUG) cout << "(" << p[i].first << ", " <<  p[i].second << ") ";
        }
        if(m_verbose == __VER_DEBUG) cout << endl;
    
        TF1 * f = new TF1("lin_temp_surr", "[0]*x + [1]", 0, 100000);
        f->SetParameters(1., 1.);
        // "N"  Do not store the graphics function, do not draw
        // "R"  Use the Range specified in the function range
        // "M"  More. Improve fit results.
        // "S"  The result of the fit is returned in the TFitResultPtr
        // "Q"  Quiet
        TString fitconfig = "RNQS";
        if(m_verbose == __VER_DEBUG) fitconfig = "RNS";
        g->Fit(f, fitconfig.Data(), "");
    
        // Return values
        a = f->GetParameter(0);
        b = f->GetParameter(1);
    
        // delete the objects
        delete f;
        delete g;
	}

}

void TOTCalib::RandomFitParametersSurrogate (TF1 * f, double a, double b) {

	// p0 is a (a from Am241 +/- 2)
	f->SetParameter(0, ( (2. - a) * m_rand1->Rndm() / 4. ) );
	// p1 is b (b from Am241 +/- 50)
	f->SetParameter(1, ( (50. - b) * m_rand1->Rndm() / 100. ));
	// p2 is c ]10 500]
	f->SetParameter(2, ( m_rand1->Rndm()*490 ) + 10 );
	// p3 is t ]1 5]
	f->SetParameter(3, ( m_rand1->Rndm()*4 ) + 1 );

}

/**
 * Fill conviniently a point in the fit function parameter space
 * to start the fit. Some of the vars will be random.
 */

void TOTCalib::RandomFitParameters (TF1 * f, TH1 * h, int tot, TOTCalib* s) {

	//double * pars = new double (__npars_lowe_fitfunc);
	double loc_bandwidth = s->GetKernelBandWidth();
	if(m_verbose == __VER_DEBUG_LOOP) { cout << "[RAND] Producing a set of random fit parameters." << endl; }    
    
	if ( TString(f->GetName()).Contains("gf_lowe") ) { // FIXME: this case is Ad hoc

		//pars = new double (__npars_lowe_fitfunc);
		// return random numbers in ]0,1]
		//m_rand1->RndmArray(__npars_lowe_fitfunc, pars);

        // p0 is the amplitude and can always be taken from the histogram
        // For now --> (20, 50)
        f->SetParameter(0, m_rand1->Rndm()*30. + 20. );
        // p1 is the mean --> fixed in PeakFit()
        // p2 is the sigma 
        f->SetParameter(2, 2.);
        // p3 is a. Random --> (2, 4)
        f->SetParameter(3, m_rand1->Rndm()*2. + 2.);
        // p4 is b. Random --> (20, 100)
        f->SetParameter(4, m_rand1->Rndm()*80. + 20.);
        // p5 is c. Random --> (-50, 350)
        f->SetParameter(5, m_rand1->Rndm()*100. + 150.);//175.
        // p6 is t. Random --> (1, 6)
        f->SetParameter(6,0.);
        
        //cout<<f->GetParameter(0)<<" "<<f->GetParameter(1)<<" "<<f->GetParameter(2)<<" "<<f->GetParameter(3)<<" "<<f->GetParameter(4)<<" "<<f->GetParameter(5)<<" "<<f->GetParameter(6)<<endl;
        
        
	} else if ( TString(f->GetName()).Contains("gf_linear") ) {

		// p0 is the amplitude and can always be taken from the histogram
		f->SetParameter(0, h->GetBinContent( h->FindBin(tot) ) );
		// p1 is the mean tot and can be just the tot, or somewhere close towards to the right
		f->SetParameter(1, m_rand1->Rndm()*loc_bandwidth + tot);
		// p2 is the sigma of the gaussian part.  Use the bandwidth
		f->SetParameter(2, loc_bandwidth);

	}

	//delete pars;

}

vector<double> TOTCalib::LowStatsPeakSelection( vector<double> peaks,unsigned int size, TOTCalib * source, double bandwidth){ //EDIT:JS019
	vector<double> global_p = source->GetGlobalMaximumPoints();
	vector<double> selected_p;

	unsigned int i = 0;

	if (m_verbose != __VER_QUIET){
		cout << "[INFO] ----- Low Stats Peak Selection -----" << endl;
		cout << "Expected " << size << " peaks. | Peaks found : " ;
		for (vector<double>::iterator k=peaks.begin(); k!=peaks.end() ; k++){ cout << (*k) << "  "; }
		cout << endl << "Global kernel function's maximums : ";
		for (vector<double>::iterator l=global_p.begin(); l !=global_p.end() ; l++){ cout << (*l) << "  ";}
		cout << endl;
	}


	if (global_p.size() != size){
		cout << "[ERROR] The expected number of peaks is not equal to the number of maximums in the global kernel function!" <<endl;
		cout << "Try to work with more pixels in your next run." << endl;
		return selected_p; // return empty
	}

	for ( ; i<global_p.size() ; i++){
		unsigned int j = 0;
		vector<double> diff;

		for ( ; j<peaks.size(); j++){
			diff.push_back(TMath::Abs( peaks[j] - global_p[i] ));
		}

		int min_ind = distance (diff.begin(), min_element(diff.begin(),diff.end()) );
		
        if (diff[min_ind] > bandwidth*2.){
            cout << "[ERROR] Selected peak is too far from the global kernel function's maximum (bandwidth*2 = "<< bandwidth*2 <<"). Selection cancelled." << endl;
			return selected_p; // return empty
		}

		selected_p.push_back(peaks[min_ind]);
	}

	if (m_verbose != __VER_QUIET){
		cout << "Selected peaks : ";
		for (vector<double>::iterator m = selected_p.begin() ; m!= selected_p.end() ; m++){cout << (*m) << "  ";}
		cout << endl;
	}

	return selected_p;
}

int TOTCalib::PeakFit(TOTCalib * src, int /*pix*/, int tot, TF1 * f, TH1 * h, store * sto, double energy) {

	// This was a suggestion but didn't quite work correctly.
	// Use a,b from previous fits in the linear region.
	// double a = 0., b = 0.;
	// GetLinearFit(a, b, sto->linearpairs);
	// if(m_verbose == __VER_DEBUG) cout << " -- NOT USED -- a = " << a << ", b = " << b << endl;

	double loc_bandwidth = src->GetKernelBandWidth();
	double minf;
	double maxf;

	if (src->GetCalibMethod() == __lowStats){ // optimized fit interval
		minf = tot - 2*loc_bandwidth;
		maxf = tot + 2*loc_bandwidth;

	} else { // standard fit interval

		minf = tot - loc_bandwidth;
		maxf = tot + loc_bandwidth;

		int centerBin = h->FindBin(tot);
		double heightAtKernelHint = h->GetBinContent(centerBin);

		unsigned int nBins = h->GetXaxis()->GetNbins();
		unsigned int rightBin = h->FindBin( maxf );
		unsigned int leftBin = h->FindBin( minf );

		int underThresholdCntr = 0;
		for(unsigned int i=centerBin; i<=nBins; ++i ) {
			if( h->GetBinContent(i) < heightAtKernelHint * __fraction_of_height_range_id ) {
				underThresholdCntr++;
				rightBin = i;
				if(underThresholdCntr > 1) break;
			}
		}
		underThresholdCntr = 0;
		for(unsigned int i=centerBin; i>1; --i ) {
			if( h->GetBinContent(i) < heightAtKernelHint * __fraction_of_height_range_id ) {
				underThresholdCntr++;
				leftBin = i;
				if(underThresholdCntr > 1) break;
			}
		}
		// take one more
		if(rightBin < nBins) rightBin++;
		if(leftBin > 1) leftBin--;

		// define maxf and minf
		maxf = h->GetBinCenter(rightBin);
		minf = h->GetBinCenter(leftBin);

		//int width = rightBin - centerBin;
	}


	if( TString(f->GetName()).Contains("gf_lowe") ) {
		minf = 1;
		maxf = tot + loc_bandwidth*3;
	}
	if(minf < 1) minf = 1; // correct for negative or zero minf value

	TString fitconfig = "NQS";
	if(m_verbose == __VER_DEBUG_LOOP) fitconfig = "NS";

	// "N"  Do not store the graphics function, do not draw
	// "R"  Use the Range specified in the function range
	// "M"  More. Improve fit results.
	// "S"  The result of the fit is returned in the TFitResultPtr
	// "Q"  Quiet

	if( (m_verbose >= __VER_INFO) && (m_verbose != __VER_QUIET) ) {
		cout << "[FIT] Fit in the interval : " << minf << ", " << maxf << " with options : " << fitconfig << endl;
	}

	// Save this histogram FIXME
	//TString extrasavename = h->GetName();
	//extrasavename += "_tosave";
	//m_histosToSave.push_back( static_cast <TH1 *>( h->Clone( extrasavename ) ) );

	int status = -1, fittries = 0;
	double sprob = 0.;
	int fit_max_rand_tries = 0;

	// Save all the tries and if the limit __fit_pars_randomization_max is reached just pick up the best
	map<int, vector<double> > calibTriesMap;
	map<int, vector<double> > calibTriesErrMap;
	vector<double> calibTriesProb;
	vector<int> calibTriesStatus;
	int indexmax = 0; // Index of best fit
	bool continueFitting = true;

    Double_t a = 0.;
    Double_t b = 0.; 
    
    if (TString(f->GetName()).Contains("gf_lowe")) {
        if (sto->linearpairs.size() >= 2){
            GetLinearFit(a,b,sto->linearpairs);        
            if( m_verbose != __VER_QUIET ) {
                cout<< "[FIT] Using a and b coeff from linear region: a="<<a<<" and b="<<b<<endl;
            }
        }else{
            if( m_verbose != __VER_QUIET ) {
                cout<< "[FIT] Warning: trying to fit lowen func without a and b parameters fixed"<<endl;
            }
        }

    }

	//while ( sprob < __min_tmathprobtest_val && fit_max_rand_tries < __fit_pars_randomization_max ) {
	while ( continueFitting ) {

		RandomFitParameters( f, h, tot, src );

        // keep track of the parameters to choose the better set in case the fit is not good enough

        if ( TString(f->GetName()).Contains("gf_lowe") ){
           if (a!=0. && b!=0.){
               f->FixParameter(1,energy);
               f->FixParameter(3,a);
               f->FixParameter(4,b);
           }else {
               f->FixParameter(1,energy); 
           }
        }

		fittries = 0; // rewind
		status = -1;
		while ( status != 0 && fittries < __max_fit_tries ) {
		
			if (src->GetCalibMethod() == __lowStats){
				f->FixParameter(1,tot); // fix mean of the gaussian (otherwise incorrect fit in most cases)
				f->SetParLimits(2,0, 1.5*loc_bandwidth); // set a max value for sigma (otherwise incorrect in most cases)
			}

			TFitResultPtr fitr = h->Fit(f, fitconfig.Data(), "" , minf, maxf);
			status = int ( fitr );
			// special case where something very bad happens like trying to fit with empty data
			if( status == -1 ) break;

			// Get the stat prob of the fit
			sprob = TMath::Prob( fitr.Get()->Chi2() / fitr.Get()->Ndf() , fitr.Get()->Ndf() );

			fittries++;
		}

		// special case where something very bad happens like trying to fit with empty data
		if( status == -1 ) break;

		if(m_verbose == __VER_DEBUG_LOOP) cout << "Try " << fit_max_rand_tries << " { status : " << status << " } [PROB] = " << sprob << " " << endl;

		// Save the results of this fit
		calibTriesMap[fit_max_rand_tries] = vector<double>(f->GetNpar(), 0.); // Save as many parameters as the fit has.  It does change.
		calibTriesErrMap[fit_max_rand_tries] = vector<double>(f->GetNpar(), 0.);
		calibTriesProb.push_back( sprob );
		calibTriesStatus.push_back( status );
		// Recover the resulting fit parameters
		for(int i = 0 ; i < f->GetNpar() ; i++) {
			calibTriesMap[fit_max_rand_tries][i] = f->GetParameter(i);
			calibTriesErrMap[fit_max_rand_tries][i] = f->GetParError(i);
		}

		fit_max_rand_tries++;

		// Decide wether we try another fit of not
		if ( sprob < __min_tmathprobtest_val && fit_max_rand_tries < __fit_pars_randomization_max ) { // keep trying if the minimum has not been reached up to a maximum number of tries
			continueFitting = true;
		} else if ( sprob > __min_tmathprobtest_val && fit_max_rand_tries < __fit_pars_randomization_min) { // even if the minimum was found try a minimum number of times
			continueFitting = true;
		} else { // otherwise stop
			continueFitting = false;
		}
	}

	//if ( fit_max_rand_tries > __fit_pars_randomization_max ) { // If this happens the minimum C.L. for the fit was not reached.  Find the best fit.

	// Find the best fit

	// Search the max
	vector<double>::iterator iB = calibTriesProb.begin();
	vector<double>::iterator iE = calibTriesProb.end();
	vector<double>::iterator imax = max_element( iB, iE );

	// See what the index is
	vector<double>::iterator i = iB;
	for( ; i != iE ; i++) {
		if( i == imax ) break;
		indexmax++;
	}

	// Print out the best fit values
	if(m_verbose == __VER_DEBUG) {
		cout << "[MAX] Fit = " << f->GetName() << " | best fit index = " << indexmax << " | prob = "
				<< calibTriesProb[indexmax] << " | status : " << calibTriesStatus[indexmax] << " | pars : ";
		for(int j = 0 ; j < f->GetNpar() ; j++) {
			cout << f->GetParName(j) << " = " << calibTriesMap[indexmax][j] << "+/-" << calibTriesErrMap[indexmax][j];
			if( j < f->GetNpar() - 1 ) cout << ", ";
		}
		cout << endl;
	}

	// Change the fit parameters in the fitting function to the best found.
	// The parameters are return to the calling function using the TF1 pointer.
	for ( int j = 0 ; j < f->GetNpar() ; j++ ) {
		f->SetParameter( j, calibTriesMap[indexmax][j] );
	}

	// save the status
	sto->peakFitStatus.push_back( calibTriesStatus[indexmax] );

	//} else {
	// Otherwise the fit reached the expected C.L. and it will be taken just as it is
	// Save the status
	//	sto->peakFitStatus.push_back( status );

	//}

	// See if after all the goal was reached.
	// If the statistical test wasn't reached
	// then I will fake the status of the fit.
	// This will be taken as a bad pixel.
	//if(sprob < __min_tmathprobtest_val) {
	//status = -1;
	//}


	return status;
}




void TOTCalib::ProcessOneSource(TOTCalib * s, store * sto, TGraphErrors * g, int pix, int & cntr) {

	// Before proceeding with the sets of points identified
	//  I need a gaussian fit in each peak and extract the sigma.
	//  This sigma will be the error in the TGraphError to fit.

	// Sets of points
	vector< pair<double, double> > points = Extract_E_TOT_Points ( pix, s ) ; // energies in this vector are in absolute value (for peak order)
	vector<pair<double, double> >::iterator i;

	// region type, linear, low energy, undefined
	map<int, int> region = s->GetCalibHandler()->GetCalibPointsRegion();
	map<int, int>::iterator regionItr = region.begin();
    
	map<int, double> Epoint = s->GetCalibHandler()->GetCalibPoints(); // energies in this vector may be < 0 (to skip)
	map<int, double>::iterator EpointItr = Epoint.begin();

	// The selected fitting function.  Do not delete in this scope !
	TF1 * gf;

	// Obtain the histogram for this pixel
	int totval = 0;
	double totmeanfit = 0., sigmafit = 0., constantfit = 0.;

	// The data histogram
	TH1I * hf = s->GetHisto(pix, "blender");

	if(m_verbose == __VER_DEBUG) {
		cout << "After Fit, points for : " <<  s->GetCalibHandler()->GetSourcename() <<  "  |  ";
	}

	int calibPointIterator = 0;
	for ( i = points.begin() ; i != points.end(); i++ ) {

		// Make the fit around the peak
		totval = (*i).second;

		// Select the appropiated function
		gf = FittingFunctionSelector( (*i).first, s , calibPointIterator );
		if(m_verbose == __VER_DEBUG) cout << " [ fit func --> " << gf->GetName() << "] ";

		// Fit in the peak
        int status = 0;
        if ( TString(gf->GetName()).Contains("gf_lowe") ){
            status = PeakFit(s, pix, totval, gf, hf, sto, (*i).first);
        }else{
            status = PeakFit(s, pix, totval, gf, hf, sto);
        }
        if(m_verbose == __VER_DEBUG) cout << " { status : " << status << " } ";
        Double_t func_TOTatMax = gf->GetMaximumX();                    

		// These are the points for the surrogate function fit
		constantfit = gf->GetParameter(0);
		totmeanfit = gf->GetParameter(1);
		sigmafit = TMath::Abs ( gf->GetParameter(2) );
        
        if ( (*EpointItr).second > 0. ){ // only fill graph if positive energy (the graph is used for surrogate fit in Blender)
            
            if ( TString(gf->GetName()).Contains("gf_lowe") ) {
                g->SetPoint(cntr, (*i).first, func_TOTatMax ); // for gf_lowe the tot point is not a gaussian mean but the tot (x coordinate) at the function maximum
                g->SetPointError(cntr, 0., sigmafit );
            }else{
                g->SetPoint(cntr, (*i).first, totmeanfit ); // E, TOT (from the fit in this context)
                g->SetPointError(cntr, 0., sigmafit );
            }
            cntr++; 
        }


		sto->pointsSaveSigmas.push_back( sigmafit );                       // The sigma of the fit
		sto->pointsSaveConstants.push_back( constantfit );                 // The constant of the fit
		sto->calibTOTPeaks.push_back( totval );                            // The original TOT val where the fit starts
		sto->peakFitStatus.push_back( status );

		if( TString(gf->GetName()).Contains("gf_lowe") ) {  // in this case store the extra params
            sto->pointsSave.push_back( make_pair( (*EpointItr).second, func_TOTatMax ) );            
			sto->pointsSave_ia.push_back( gf->GetParameter(3) );
			sto->pointsSave_ib.push_back( gf->GetParameter(4) );
			sto->pointsSave_ic.push_back( gf->GetParameter(5) );
			sto->pointsSave_it.push_back( gf->GetParameter(6) );
		} else {
            sto->pointsSave.push_back( make_pair( (*EpointItr).second, totmeanfit ) );  // The mean of the fit             
			sto->pointsSave_ia.push_back( 0. );
			sto->pointsSave_ib.push_back( 0. );
			sto->pointsSave_ic.push_back( 0. );
			sto->pointsSave_it.push_back( 0. );
		}

		//////////////////////////////////////////////////////////////////
		// Verify if this fit has been requested for the linear region
		// cout << "[REGION] Calib region : " << (*regionItr).first << " --> " << (*regionItr).second << endl;
		if ( (*regionItr).second == CalibHandler::__linear_reg && (*EpointItr).second > 0. ) {
			sto->linearpairs.push_back( make_pair( (*i).first , totmeanfit ) );  // (E, TOT)
		}
		regionItr++;
        EpointItr++;

		//////////////////////////////////////////////////////////////////

		if(m_verbose == __VER_DEBUG) {
			cout << " (" << (*i).first << " , " << totmeanfit << ") ";
		}

		calibPointIterator++;
	}
	if(m_verbose == __VER_DEBUG) cout << endl;
	delete hf;

}

/**
 *  Reorder the sources used for calibration in such a way
 *  that the first source processed gives at least two points
 *  in the linear region which aid the distribution fits in
 *  the non-linear region.
 */
void TOTCalib::ReorderSources() {

	vector<TOTCalib *>::iterator i = m_allSources.begin();
	map<int, int>::iterator reg_i;
	int nlinear = 0;
    int nlowenergy = 0;
    int nSources = m_allSources.size();

	// Indexes of sources that must be processed first and last
	int firstIndex = 0;
	int lastIndex = 0;

	int index = 0;
	for( ; i != m_allSources.end() ; i++ ) {

		nlinear = 0; // rewind
        nlowenergy = 0; // rewind        

		map<int, int> reg = (*i)->GetCalibHandler()->GetCalibPointsRegion();
		reg_i = reg.begin();
		for ( ; reg_i != reg.end() ; reg_i++ ) {
			if( (*reg_i).second == CalibHandler::__linear_reg ) nlinear++;
            if( (*reg_i).second == CalibHandler::__lowenergy_reg ) nlowenergy++;
		}

		if ( nlinear >= 2 ) {    // this can be the first source
			firstIndex = index;
		}
        
        if ( nlowenergy >= 1 ) {    // this must be at the end
			lastIndex = index;
		}       

		index++;
	}

	// Switch positions with the pointer that should be processed first
	TOTCalib * temp_Ptr = m_allSources[0];
	m_allSources[0] = m_allSources[firstIndex];
	m_allSources[firstIndex] = temp_Ptr;
    
    // Switch positions with the pointer that should be processed last
	TOTCalib * temp_Ptr2 = m_allSources[nSources-1];
	m_allSources[nSources-1] = m_allSources[lastIndex];
	m_allSources[lastIndex] = temp_Ptr2;

	// Report order
	cout << "Processing sources in the following order : ";
	for(i = m_allSources.begin() ; i != m_allSources.end() ; i++ ) {
		cout << (*i)->GetCalibHandler()->GetSourcename() << " ";
		if( i+1 !=  m_allSources.end() ) cout << ", ";
	}
	cout << endl;
    cout<<"-----------------------------------------------------------"<<endl;    
}

void TOTCalib::CreateGlobalKernelAndGetCriticalPoints(){ 


	//Create histogram (a vector, not a TH1) for the global spectrum of each source (whole map, not single pixels)
	vector<TOTCalib *>::iterator itr = m_allSources.begin();
	for ( ; itr!=m_allSources.end() ; itr++){
		
		int method = (*itr)->GetCalibMethod();
		if (method != __lowStats) continue;

		if (m_verbose !=__VER_QUIET) {
            cout << "[INFO]	Gathering critical points for the whole pixel selection " << endl;
        }
		int nbins = (*itr)->GetNBins();

        vector<vector<double> > m_histo = (*itr)->Get_m_histo();
		vector<double> global_ker(nbins,0.);

		for (int pix = m_minpix ; pix <= m_maxpix ; pix++) {
			vector<double> temp = m_histo[pix];
			vector<double>::iterator i = temp.begin();
			int cntr = 0;

			for ( ; i!=temp.end() ; i++){
				global_ker[cntr] += *i;
				cntr++;
			}
		}

		(*itr)->SetGlobalHisto(global_ker);

		//from CreateKernelDensityFunction()
		//create the global kernel density function for each source

		double * par = new double[nbins * 3 + 1];
		par[0] = nbins;
		double bandwidth = (*itr)->GetKernelBandWidth();

		int j = 1;
		for ( ; j <= nbins ; j++ ) {
			par[j] = global_ker[j-1] / bandwidth; // constant
			par[j + nbins] = j-1; // mean
			par[j + 2*nbins] = bandwidth; // sigma
		}

		TString kernelname = "global_kernel_";
		kernelname += (*itr)->GetCalibHandler()->GetSourcename();
		TF1 * fker = new TF1(kernelname, GausFuncAdd, 0, nbins, 3 * nbins + 1);
		fker->SetParameters( par );

		//from GetCriticalPoints()
		//find the maximum and minimum of the global kernel

		short sign = __s_pos;
	
		double der = 0.;
		double step = 1.;
		vector<double> temp_max;
		vector<double> temp_min;
		for (double x = 1. ; x <= (double)nbins ; x+= step) {
	
			//der = f->Derivative(x, 0x0, 0.1);
			der = DerivativeFivePointsStencil(fker, x, 0.1);
	
			if(der < 0 && sign == __s_pos) { // Flip from pos to neg --> maximum
				temp_max.push_back( x ); // Critical point here
			}
			if(der > 0 && sign == __s_neg) { // Flip from neg to pos --> minimum
				temp_min.push_back( x );  // Critical point here
			}
	
			if(der > 0) sign = __s_pos;
			if(der < 0) sign = __s_neg;

		}

		(*itr)->SetGlobalCriticalPoints(temp_max,temp_min);

		if (m_verbose !=__VER_QUIET) {
            cout << "Source : " << (*itr)->GetCalibHandler()->GetSourcename() << endl;
            cout << "Maximum (" <<temp_max.size()<< ") : ";
            for (vector<double>::iterator k = temp_max.begin() ; k!=temp_max.end() ; k++){ cout << *k << " " ;}

            cout << endl << "Minimum (" <<temp_min.size()<< ") : ";
        	for (vector<double>::iterator k = temp_min.begin() ; k!=temp_min.end() ; k++){ cout << *k << " " ;}
        	cout << endl << "-------------------------" << endl;
        }
	
		delete fker; 
		temp_min.clear(); temp_max.clear();

	}
}

void TOTCalib::Blender (TOTCalib * s2, TOTCalib * s3, TOTCalib * s4, TOTCalib * s5, TString outputName, int calibMethod) {

    cout<<endl<<"-----------------------------------------------------------"<<endl;    
	cout << "Blender ... making all the fits" << endl;

	// The first source to be processed needs to be the source where two,
	// or more points defining the linear region are to be found
	m_allSources.push_back( this );
	m_allSources.push_back( s2 );
	m_allSources.push_back( s3 );
	m_allSources.push_back( s4 );
    m_allSources.push_back( s5 );

	Blender( outputName, calibMethod );
}

void TOTCalib::Blender (TOTCalib * s2, TOTCalib * s3, TOTCalib * s4, TString outputName, int calibMethod) {

    cout<<endl<<"-----------------------------------------------------------"<<endl;    
	cout << "Blender ... making all the fits" << endl;

	// The first source to be processed needs to be the source where two,
	// or more points defining the linear region are to be found
	m_allSources.push_back( this );
	m_allSources.push_back( s2 );
	m_allSources.push_back( s3 );
	m_allSources.push_back( s4 );

	Blender( outputName, calibMethod );
}

void TOTCalib::Blender (TOTCalib * s2, TOTCalib * s3, TString outputName, int calibMethod) {

    cout<<endl<<"-----------------------------------------------------------"<<endl;    
	cout << "Blender ... making all the fits" << endl;

	// The first source to be processed needs to be the source where two,
	// or more points defining the linear region are to be found
	m_allSources.push_back( this );
	m_allSources.push_back( s2 );
	m_allSources.push_back( s3 );

	Blender( outputName, calibMethod );
}

void TOTCalib::Blender (TOTCalib * s2, TString outputName, int calibMethod) {
    
    cout<<endl<<"-----------------------------------------------------------"<<endl;    
	cout << "Blender ... making all the fits" << endl;

	// The first source to be processed needs to be the source where two,
	// or more points defining the linear region are to be found
	m_allSources.push_back( this );
	m_allSources.push_back( s2 );

	Blender( outputName, calibMethod );
}


void TOTCalib::Blender (TString outputName, int calibMethod) {
    
    if (calibMethod==__jakubek){
        cout<<"You have chosen the calibration method described in NIMPR A 633 (2011) S262–S266"<<endl;
        cout<<"You should have given one point in the low energy region and at least two points in the linear region."<<endl;
        cout<<"The a and b coeff will be determined by a fit in the linear region."<<endl;
        cout<<"The c and t coeff will be determined by a fit with the low energy function to the last spectrum in the order list."<<endl;
        cout<<"The threshold energy will be ignored, and no fit with the surrogate function will be done."<<endl; 
    }
    
	ReorderSources();
	CreateGlobalKernelAndGetCriticalPoints(); // will only be used for sources with m_method == __lowStats
    
	// Files to write
	TString fn_a = outputName + "_a.txt";
	TString fn_b = outputName + "_b.txt";
	TString fn_c = outputName + "_c.txt";
	TString fn_t = outputName + "_t.txt";
	TString fn_prob = outputName + "_prob.txt";

	ofstream f_a(fn_a, ostream::out);
	ofstream f_b(fn_b, ostream::out);
	ofstream f_c(fn_c, ostream::out);
	ofstream f_t(fn_t, ostream::out);
	ofstream f_prob(fn_prob, ostream::out);

	// Fit function for distributions (not the surrogate)
	// "gf" is just a gaussian and is good in all the linear region of the calibration
	// "gf_lowe" is a gaussian + (surr)^-1.  See fitfunc_lowen.

	// first find who has the biggest range
	vector<TOTCalib *>::iterator iir = m_allSources.begin();
	double maxrange = 0.;
	for ( ; iir != m_allSources.end() ; iir++ ) {
		if( (*iir)->GetNBins() > maxrange ) maxrange = (double) ( (*iir)->GetNBins() );
	}

	m_gf_linear = new TF1("gf_linear", "gaus(0)", 0., maxrange);
	m_gf_linear->SetParameters(1, 1, m_bandwidth);

	m_gf_lowe = new TF1("gf_lowe", fitfunc_lowen, 0., maxrange, __fitfunc_lowen_npars);
	m_gf_lowe->SetParameters(1, 1, m_bandwidth, 1,1,1,1);

	double percentage = 0.;
	// iterate over pixels
	for (int pix = m_minpix ; pix <= m_maxpix ; pix++) {

		// Skip the bad pixels
		if ( PixelInBadPixelList(pix) ) continue;

		// Vectors to save the fit constants and properties
		vector<double> calibConst;
		vector<double> calibProperties;

		// Name for the surrogate function
		TString fn = "surr_pix_";
		fn += pix;

		int totalNPoints = 0;
		vector<TOTCalib *>::iterator i;
		for(i = m_allSources.begin() ; i != m_allSources.end() ; i++ ) {
			totalNPoints += GetNumberOf_E_TOT_Points_Positive( *i );
		}

		// Set of vectors used to store info
		store * st = new store;
		TGraphErrors * g = 0x0;
		// Values for linear fit
		double a = 1., b = 1.;
		// Save all the tries and if the limit __fit_pars_randomization_max is reached just pick up the best
		map<int, vector<double> > calibTriesMap;
		vector<double> calibTriesProb;
		vector<int> calibTriesStatus;
		int indexmax = 0; // Index of best fit

		// Proceed with the local fits
		if ( totalNPoints > 0 ) {

			// Create an object with this points and perform a fit.  A TGraph.
			g = new TGraphErrors(totalNPoints+1); //1 point for every peak expected + detector treshold

			// Counter for points in TGraphErrors
			int cntr = 0;

			// Process source in the defined order
			vector<TOTCalib *>::iterator i;
            
            if(m_verbose != __VER_QUIET) {            
                cout<<endl<< "**************** Processing pixel : "<<pix<<" **************** "<<endl;
            }/*else{
                if (pix % 1000 == 0) cout<<endl<<"Processing pixel : "<<pix<<endl;
            }    */      
            
			for (i = m_allSources.begin() ; i != m_allSources.end() ; i++ ) {

				ProcessOneSource(*i, st, g, pix, cntr);
				if(m_verbose == __VER_DEBUG_LOOP) {
					cout << "Got " << st->linearpairs.size() << " pairs in the linear part." << endl;
				}

			}

			// Get the linear part to help the fit later
			GetLinearFit( a, b, st->linearpairs );

            if (calibMethod != __jakubek){
                
                // Global threshold from THL calibration
                g->SetPoint( cntr, m_thresholdEnergy, 0.0 );        // E = threshold_energy --> 0 TOT counts
                g->SetPointError(cntr, m_thresholdEnergy_Err, 0.0 );
                cntr++;
            }            
		}

		// Check if any of the sources had a fit status -1 := no data
		vector<int> allstatus = st->peakFitStatus;
		bool missingInfo = false;
		if ( find(allstatus.begin(), allstatus.end(), -1) !=  allstatus.end() ) { // found -1
			missingInfo = true;
		}


		// If all fits were good
		if( totalNPoints > 0 && !missingInfo ) {

			// store points
			m_calibTOTPeaks[pix] = st->calibTOTPeaks;
			m_calibPoints[pix] = st->pointsSave;
			m_calibPointsSigmas[pix] = st->pointsSaveSigmas;
			m_calibPointsConstants[pix] = st->pointsSaveConstants;
			m_calibPoints_ia[pix] = st->pointsSave_ia;
			m_calibPoints_ib[pix] = st->pointsSave_ib;
			m_calibPoints_ic[pix] = st->pointsSave_ic;
			m_calibPoints_it[pix] = st->pointsSave_it;

            if (calibMethod == __standard){
                      
                // [0] --> a
                // [1] --> b
                // [2] --> c
                // [3] --> t
                TF1 * surr = new TF1(fn, surrogatefunc_calib, m_thresholdEnergy, 60.0, 4); // range in keV, 4 parameters
                //TF1 * surr = new TF1(fn, surrogatefunc_calib, 0.0, 60.0, 4); // range in keV, 4 parameters
    cout<<m_thresholdEnergy<<" "<<endl;
    g->Print();
                //surr->SetParameters( 1, 1, 100, 1 );
                //surr->SetParLimits(2, 10, 250);
                //surr->SetParLimits(3, 0.1, 5.0);
                
                // "N"  Do not store the graphics function, do not draw
                // "R"  Use the Range specified in the function range
                // "M"  More. Improve fit results.
                // "S"  The result of the fit is returned in the TFitResultPtr
                // "E"  Perform better Errors estimation using Minos technique
                // "Q"  Quiet
                int status = -1; int initFit = 0;
                TString fitconfig = "NRSQ";
                if(m_verbose == __VER_DEBUG_LOOP) fitconfig = "NRS";
                double sprob = 0.0;
                int fit_max_rand_tries = 0;
    
    
                while ( sprob < __min_tmathprobtest_val && fit_max_rand_tries < __fit_pars_randomization_max ) {
    
                    RandomFitParametersSurrogate( surr, a, b );
                    initFit = 0;
                    status = -1;
    
                    while ( status != 0 && initFit < __max_fit_tries ) {
    
                        TFitResultPtr fitr = g->Fit(fn, fitconfig.Data(), "");
                        TFitResult * fitp = fitr.Get();
                        status = fitp->Status();
                        sprob = TMath::Prob( surr->GetChisquare()/surr->GetNDF(), surr->GetNDF() );
                        initFit++;
    
                    }
    
                    if(m_verbose == __VER_DEBUG_LOOP) cout << "Try " << fit_max_rand_tries << " { status : " << status << " } [PROB] = " << sprob << " " << endl;
    
                    // Save the results of this fit
                    calibTriesMap[fit_max_rand_tries] = vector<double>(__npars_surrogate, 0.);
                    calibTriesProb.push_back( sprob );
                    calibTriesStatus.push_back( status );
                    for(int i = 0 ; i < __npars_surrogate ; i++) {
                        calibTriesMap[fit_max_rand_tries][i] = surr->GetParameter(i);
                    }
    
                    fit_max_rand_tries++;
                }
    
                // If the search reached the maximum numbers of tries
                vector<double> surr_pars(__npars_surrogate, 0.);
    
                if(fit_max_rand_tries == __fit_pars_randomization_max) {
    
                    // Search the max
                    vector<double>::iterator iB = calibTriesProb.begin();
                    vector<double>::iterator iE = calibTriesProb.end();
                    vector<double>::iterator imax = max_element( iB, iE );
                    // See what the index is
                    vector<double>::iterator i = iB;
                    for( ; i != iE ; i++) {
                        if( i == imax ) break;
                        indexmax++;
                    }
    
                    // This is the best set of parameters found
                    surr_pars[0] = calibTriesMap[indexmax][0];
                    surr_pars[1] = calibTriesMap[indexmax][1];
                    surr_pars[2] = calibTriesMap[indexmax][2];
                    surr_pars[3] = calibTriesMap[indexmax][3];
    
                    if (m_verbose !=__VER_QUIET){
                        cout<<"------- Calibration results -------"<<endl;
                        cout << "[MAX] Pixel = " << pix << " | best fit index = " << indexmax << " | prob = ";
                        cout << calibTriesProb[indexmax] << " | status : " << calibTriesStatus[indexmax];
                        cout << " | a, b, c, t : ";
                        cout << calibTriesMap[indexmax][0] << ", ";
                        cout << calibTriesMap[indexmax][1] << ", ";
                        cout << calibTriesMap[indexmax][2] << ", ";
                        cout << calibTriesMap[indexmax][3] << endl;

                    }
    
                } else {
    
                    surr_pars[0] = surr->GetParameter(0);
                    surr_pars[1] = surr->GetParameter(1);
                    surr_pars[2] = surr->GetParameter(2);
                    surr_pars[3] = surr->GetParameter(3);
    
                }
    
    
                //cout << "Status : " << status << endl;
                //cout << "p value of the fit : " << fitp->Prob() << endl;
    
                // Save the fit constants
                if( ( fit_max_rand_tries == __fit_pars_randomization_max && calibTriesStatus[indexmax] == 0 ) 	 // Fit which never reached the requested C.L.
                        || status == 0 )                                                                         // Fit which converged with high C.L.
                {
                    calibConst.push_back( surr_pars[0] );
                    calibConst.push_back( surr_pars[1] );
                    calibConst.push_back( surr_pars[2] );
                    calibConst.push_back( surr_pars[3] );
    
                    calibProperties.push_back( calibTriesProb[indexmax] );
    
                } else {
    
                    calibConst.push_back( 0.0 );
                    calibConst.push_back( 0.0 );
                    calibConst.push_back( 0.0 );
                    calibConst.push_back( 0.0 );
    
                    calibProperties.push_back( 0.0 );
    
                }
                
                // delete the function and TGraph
                delete surr;
                                
            } else if (calibMethod == __jakubek){  // In this case no fit with surrogate is needed
              
                int size = st->pointsSave_ic.size();
                
                if ( size != 0 ){
                    calibConst.push_back( a );
                    calibConst.push_back( b );
                    
                    // For c and t, look for the results obtained from low energy fits                
                    double c = 0.;
                    double t = 0.;
                    vector< double >::iterator st_itr;
                    int counter = 0;
                    int position = 0;
                    for(st_itr = st->pointsSave_ic.begin() ; st_itr != st->pointsSave_ic.end() ; st_itr++ ) {
                        
                        if ((*st_itr)!=0.){c = *st_itr; position = counter;}
                        counter ++;
                    }
                    t = st->pointsSave_it.at(position);
                    calibConst.push_back( c );
                    calibConst.push_back( t );
    
                    calibProperties.push_back( 0.0 );
                    calibTriesProb.push_back( 0. );
                    
                    if (m_verbose !=__VER_QUIET) {
                        cout<<"------- Calibration results -------"<<endl;
                        cout << "a, b, c, t : "<< a << ", "<< b << ", "<< c << ", "<< t << ", "<<endl;
                    }
                    
                } else{
                	calibConst.push_back( 0.0 );
                    calibConst.push_back( 0.0 );
                    calibConst.push_back( 0.0 );
                    calibConst.push_back( 0.0 );
                	calibProperties.push_back( 0.0 );
                	calibTriesProb.push_back( 0. );

                }
            }           

		} else {

			// This will get print and interrupt the progress bar
			// ad a couple of extra endl
			cout << "Surrogate function [ " << fn << " ]" << " can not be built.  This can be a masked pixel, noisy, or too low stats." << endl;
			// Save the fit constants
			calibConst.push_back( 0.0 );
			calibConst.push_back( 0.0 );
			calibConst.push_back( 0.0 );
			calibConst.push_back( 0.0 );

			calibProperties.push_back( 0.0 );

		}

		// Delete TGraph and store object
		if(g) delete g;
        if(st) delete st;
        
		// Fill the calib constants and properties
		m_calibSurrogateConstants[pix]  = calibConst;
		// Properties
		m_calibSurrogateProperties[pix] = calibProperties;

		if(pix%100 == 0) { // refresh progress bar every 1000 frames
			percentage = ( (double)(pix - m_minpix) / (double)(m_maxpix - m_minpix) ) * 100;
			printProgBar( (int) percentage );
		}

		/*
		cout << "Surrogate function for pixel " << pix << " has the following set of parameters (a,b,c,t) : "
				<< m_calibConstants[pix][0] << ", "
				<< m_calibConstants[pix][1] << ", "
				<< m_calibConstants[pix][2] << ", "
				<< m_calibConstants[pix][3] << endl;
		 */

		// Write to files
		/*
		f_a << m_calibSurrogateConstants[pix][0] << " ";
		f_b << m_calibSurrogateConstants[pix][1] << " ";
		f_c << m_calibSurrogateConstants[pix][2] << " ";
		f_t << m_calibSurrogateConstants[pix][3] << " ";
		f_prob << m_calibSurrogateProperties[pix][0] << " ";

		if ( pix > 0 && pix % __matrix_width == 0 ) { // introduce a \n every 256 pixels
			f_a << endl;
			f_b << endl;
			f_c << endl;
			f_t << endl;
			f_prob << endl;
		}
		 */

		// Store info for final histograms
		m_par_a_v.push_back( m_calibSurrogateConstants[pix][0] );
		m_par_b_v.push_back( m_calibSurrogateConstants[pix][1] );
		m_par_c_v.push_back( m_calibSurrogateConstants[pix][2] );
		m_par_t_v.push_back( m_calibSurrogateConstants[pix][3] );
		m_surr_prob_v.push_back( calibTriesProb[indexmax] );

	}

	// Send the whole matrix to the output files
	for( int pixout = 0 ; pixout < __matrix_size ; pixout++ ) {

		if ( pixout > 0 && pixout % __matrix_width == 0 ) { // introduce a \n every 256 pixels
			f_a << endl;
			f_b << endl;
			f_c << endl;
			f_t << endl;
			f_prob << endl;
		}

		if( pixout >= m_minpix && pixout <= m_maxpix && !PixelInBadPixelList(pixout) ) {
			f_a << m_calibSurrogateConstants[pixout][0] << " ";
			f_b << m_calibSurrogateConstants[pixout][1] << " ";
			f_c << m_calibSurrogateConstants[pixout][2] << " ";
			f_t << m_calibSurrogateConstants[pixout][3] << " ";
			f_prob << m_calibSurrogateProperties[pixout][0] << " ";
		} else {
			f_a << '0' << " ";
			f_b << '0' << " ";
			f_c << '0' << " ";
			f_t << '0' << " ";
			f_prob << '0' << " ";
		}

	}


	// finish the progress bar
	printProgBar( (int) 100 );
	cout << endl;

	// close all files
	f_a.close();
	f_b.close();
	f_c.close();
	f_t.close();
	f_prob.close();

}

TF1 * TOTCalib::FittingFunctionSelector(double /*E*/, TOTCalib * s, int pointIndex){

	// Assume linear first
	TF1 * g = m_gf_linear;

	//// E in keV
	//if(E < 10.0) { g = m_gf_lowe; }

	// Pull out the calibration handler and check if a different funcion is needed
	CalibHandler * ch = s->GetCalibHandler();
	map<int, int> regions = ch->GetCalibPointsRegion();
	if ( regions[pointIndex] == CalibHandler::__lowenergy_reg ) g = m_gf_lowe;

	//cout << "[" << pointIndex << "] " << regions[pointIndex] << endl;

	return g;
}

void TOTCalib::Finalize(){

	TH1F * ha = CreateParameterHistogram(m_par_a_v, "a");
	TH1F * hb = CreateParameterHistogram(m_par_b_v, "b");
	TH1F * hc = CreateParameterHistogram(m_par_c_v, "c");
	TH1F * ht = CreateParameterHistogram(m_par_t_v, "t");
	TH1F * hprob = CreateParameterHistogram(m_surr_prob_v, "Prob");

	//TH1I * h_5000_Fe = GetHisto(50000, "Fe55_50000");

	// prepare output ROOT
	m_output_root = new TFile("output_calib.root", "RECREATE");
	m_output_root->cd();

	vector<TH1*>::iterator i = m_histosToSave.begin();
	for( ; i != m_histosToSave.end() ; i++ ) {
		(*i)->Write();
	}

	ha->Write();
	hb->Write();
	hc->Write();
	ht->Write();
	hprob->Write();
	//h_5000_Fe->Write();

	//h_5000_Fe->Write();

	m_output_root->Close();

}


TH1F * TOTCalib::CreateParameterHistogram(vector<double> v, TString name){

	double vmin = *min_element( v.begin(), v.end() );
	double vmax = *max_element( v.begin(), v.end() );

	TString hname = "h_";
	hname += name;
	int nbins = (int) ( TMath::Ceil( vmax - vmin ) * 50. );

	TH1F * h = new TH1F(hname, hname, nbins, vmin, vmax);

	vector<double>::iterator i = v.begin();
	for( ; i != v.end() ; i++ ) {
		h->Fill( *i );
	}

	return h;
}

// Uses the information from the previously fitted surrogate func "[0]*x + [1] - ([2]/(x - [3]))"
TF1 * TOTCalib::GetSurrogateFunction(int pix) {

	vector<double> vc = m_calibSurrogateConstants[pix];

	if( vc.empty() ) return 0x0;

	cout << "Request of surrogate function for pixel " << pix << " with parameters : "
			<< vc[0] << ", " << vc[1] << ", " << vc[2] << ", " << vc[3] << endl;

	TString fn = "surr_pix_";
	fn += pix;
	TF1 * surr = new TF1(fn, "[0]*x + [1] - ([2]/(x - [3]))", 0., 100.); // range in keV
	surr->SetParameters(vc[0], vc[1], vc[2], vc[3]);

	return surr;
}

// The number of expected points
int TOTCalib::GetNumberOf_E_TOT_Points (TOTCalib * s) { // returns every point defined for a source, even ignored points

	map<int, double> calibPoints = s->GetCalibHandler()->GetCalibPoints();
	return (int)calibPoints.size();

}

int TOTCalib::GetNumberOf_E_TOT_Points_Positive (TOTCalib * s) { // returns only points that are used for surrogate function's fit

	map<int, double> calibPoints = s->GetCalibHandler()->GetCalibPoints();
	int size=0;
	map<int, double>::iterator k = calibPoints.begin();
	for ( ; k != calibPoints.end(); k++ ){
		if( (*k).second>0){
			size+=1;
		}
	} 
	return (int)size;
}

void TOTCalib::PushToLowActivityList(int pix) {
	m_lowActivityList.push_back( pix );
}

void TOTCalib::PushToBadPixelList(int pix) {
	m_badPixelList.push_back( pix );
}

void TOTCalib::PushToBadPixelListRangeInclusive(int pixi, int pixf) {
	for(int i = pixi ; i <= pixf ; i++) {
		m_badPixelList.push_back( i );
	}
}


void TOTCalib::PushToBadPixelList(int x, int y) {
	PushToBadPixelList( XYtoX( make_pair(x,y), __matrix_width) );
}

bool TOTCalib::PixelInBadPixelList(int pix){

	vector<int>::iterator i = m_badPixelList.begin();

	for ( ; i != m_badPixelList.end() ; i++ ) {
		if ( *i == pix ) return true;
	}

	return false;
}

bool TOTCalib::PixelInLowActivityList(int pix){

	vector<int>::iterator i = m_lowActivityList.begin();

	for ( ; i != m_lowActivityList.end() ; i++ ) {
		if ( *i == pix ) return true;
	}

	return false;
}

// Get a point (E, TOT)
vector<pair<double, double> > TOTCalib::Extract_E_TOT_Points (int pix, TOTCalib * s ) {

	// These are the identified peaks.  There could be one more than expected which is usually artificial.
	map<int, vector<double> > s_tot = s->GetMaxPeaksIdentified();
	vector<double> peaks = s_tot[pix];

    string source_name = s->GetCalibHandler()->GetSourcename();
    
	double loc_bandwidth = s->GetKernelBandWidth(); 

	// This is coming from the kernel funciton
	// If this pixels presents low energy activity, get rid of a certain number of artificial peaks
	//if ( PixelInLowActivityList(pix) ) {
	//peaks.pop_back(); // remove the first peak
	//}

	// These are the calib points expected per source
	map<int, double> calibPoints = s->GetCalibHandler()->GetCalibPoints();

    if ( m_verbose != __VER_QUIET ){ 
        cout<<"------- Source: "<<source_name<<" -------"<<endl;
        cout<<"N peaks found = " << peaks.size() << " | Expected peaks = " << calibPoints.size() << endl;
    }

	// If this situation is present determine which peaks to remove
	// remove as many as necesary
	if (s->GetCalibMethod() == __lowStats && peaks.size()> calibPoints.size() ) peaks = LowStatsPeakSelection(peaks, calibPoints.size(), s, loc_bandwidth); // better peak selection
	while( peaks.size() > calibPoints.size() ) { // regular peak selection, should probably be modified
		TH1I * th = s->GetHisto(pix, "test");
		vector<double> integ;
		int npeaks = (int)peaks.size();
		for ( int i = 0 ; i < (int)peaks.size() ; i++ ) {
			integ.push_back( th->Integral( th->FindBin( peaks[i] - loc_bandwidth ), th->FindBin( peaks[i]+loc_bandwidth ) ) );
			if (m_verbose != __VER_QUIET){ 
                cout << th->GetName() << " thl=" << peaks[i] << " [" << i << "]=" << integ[i] << " | ";
            }
		}
		// If the last peak has low weight (data close to it, integral) then remove the last one
		if ( integ[npeaks-1] < 0.1*integ[0] ) { // it the last is less than 10% the first integral
            if (m_verbose != __VER_QUIET){ 
                cout << " | last item removed ";
			}
			peaks.pop_back();
		} else {
			// Otherwise remove the first one
            if (m_verbose != __VER_QUIET){ 
                cout << " | first item removed ";
			}
			peaks.erase( peaks.begin() );
		}

		delete th;
	}

	// Resulting points ( E , TOT )
	vector<pair<double, double> > points;

	// If there is no peaks information coming from the distribution
	//  this particular pixel can be noisy or simply masked.
	// Also is there is less peaks identified than calibration points
	//  this pixel can not be processed.
	if( peaks.empty() || peaks.size() < calibPoints.size() ) {
		return points; // return empty
	}

	if(m_verbose == __VER_DEBUG) cout << "First guess, points for : " <<  s->GetCalibHandler()->GetSourcename() <<  "  |  ";

	for (int p = 0 ; p < (int)calibPoints.size() ; p++) {

		if(m_verbose == __VER_DEBUG) cout << " (" << calibPoints[p] << " , " << peaks[p] << ") ";

		// Energy, TOT
		points.push_back(
				make_pair (
						// Energy information for that point
						TMath::Abs(calibPoints[p]), 	// must take absolute value, otherwise order may not be respected;
						// TOT 							// ignored points will be associated to wrong peaks
						peaks[p]
				)
		);
	}

	if(m_verbose == __VER_DEBUG) cout << endl;

	/*
	for (int p = 0 ; p < (int)peaks.size() ; p++) {
		// Energy, TOT
		points.push_back(
				make_pair (
						// Energy information for that point
						s->GetCalibHandler()->GetOneEnergyMatch( p ) ,
						// TOT
						peaks[p]
				)
		);
	}
	 */

	return points;
}

int TOTCalib::GetCriticalPoints(int i, vector<double> & min, vector<double> & max){

	TF1 * f = CreateKernelDensityFunction( i, m_calibhistos[i], m_bandwidth );

	int ncrit = 0;
	short sign = __s_pos;

	double der = 0.;
	double step = 1.;
	for (double x = 1. ; x <= (double)m_nbins ; x+= step) {

		//der = f->Derivative(x, 0x0, 0.1);
		der = DerivativeFivePointsStencil(f, x, 0.1);

		if(der < 0 && sign == __s_pos) { // Flip from pos to neg --> maximum
			//max.push_back( x - step );  // Critical point in the previous step
			max.push_back( x ); // Critical point here
			ncrit++;
		}
		if(der > 0 && sign == __s_neg) { // Flip from neg to pos --> minimum
			//min.push_back( x - step );  // Critical point in the previous step
			min.push_back( x );  // Critical point here
			ncrit++;
		}

		if(der > 0) sign = __s_pos;
		if(der < 0) sign = __s_neg;

		//cout << " x = " << x << " --> der = " << der << endl;
	}

	delete f;

	//cout << "ncrit = " << ncrit;
	//vector<double>::iterator i;
	//cout << " | minimums at : ";
	//for(i = min.begin() ; i != min.end() ; i++) cout << *i << "  ";
	//cout << " | maximums at : ";
	//for(i = max.begin() ; i != max.end() ; i++) cout << *i << "  ";
	//cout << endl;
	//cout << "------------------------------------------------------- " << endl;

	return ncrit;
}

double TOTCalib::DerivativeFivePointsStencil(TF1 * f, double x, double h) {

	double der = 0.;

	der -=      f->Eval( x + 2*h );
	der += 8. * f->Eval( x +   h );
	der -= 8. * f->Eval( x -   h );
	der +=      f->Eval( x - 2*h );
	der /= 12.*h;

	return der;
}

/*
double GausFuncAdd(double * x, double * par) {

	double xx = x[0];

	// very first parameter is the nbins which is the
	// same range of the kernel density function
	int N = par[0];

	double * a = new double[N]; // constant
	double * b = new double[N]; // mean
	double * c = new double[N]; // sigma

	int i = 1;
	int j = 0;
	for ( ; i <= N ; i++ ) {
 *( a + j++ ) = par[i];
	}
	j = 0;
	for ( ; i <= 2*N ; i++ ) {
 *( b + j++ ) = par[i];
	}
	j = 0;
	for ( ; i <= 3*N ; i++ ) {
 *( c + j++ ) = par[i];
	}
	j = 0;

	// sum func
	double func = 0;
	for ( i = 0 ; i < N ; i++ ) {
		func += a[i] * exp( - ( ( xx - b[i] )*( xx - b[i] ) ) / ( 2*c[i]*c[i] ) );
	}

	delete [] a;
	delete [] b;
	delete [] c;

	return func;
}
 */

double GausFuncAdd(double * x, double * par) {

	double xx = x[0];

	// very first parameter is the nbins which is the
	// same range of the kernel density function
	int N = par[0];

	// C-array "par" contains 3N + 1 parameters
	// N = par[               0 ] : Number of entries per parameter
	// a = par[       1 -->   N ] : N constants
	// b = par[   N + 1 --> 2*N ] : N mean vaues
	// c = par[ 2*N + 1 --> 3*N ] : N sigma values

	// sum func
	double func = 0;
	int Nd = 2*N;
	for ( int i = 0 ; i < N ; i++ ) {
		func += par[ i + 1 ] * exp( - ( ( xx - par[ N + i + 1 ] )*( xx - par[ N + i + 1 ] ) ) / ( 2.*par[ Nd + i + 1 ]*par[ Nd + i + 1 ] ) );
	}

	return func;
}


TF1 * TOTCalib::CreateKernelDensityFunction(int pix, vector<double> hist, double bandwidth) {

	// Set of parameters for the full kernel density function
	// Npars = m_nbins*( constants + sigmas + means ) + size
	// ex: for 100 bins we need 301 parameters
	double * par = new double[m_nbins * 3 + 1];
	par[0] = m_nbins;
	//double distmax = *max_element(hist.begin(), hist.end());

	int i = 1;
	for ( ; i <= m_nbins ; i++ ) {

		//par[i] = hist[i-1]; // constant
		//par[i] = hist[i-1] / distmax; // constant
		par[i] = hist[i-1] / bandwidth; // constant

		par[i + m_nbins] = i-1; // mean

		par[i + 2*m_nbins] = bandwidth; // sigma

	}

	// Final kernel density
	TString kernelname = "kernel_";
	kernelname += m_calhandler->GetSourcename();
	kernelname += "_";
	kernelname += pix;
	//cout << "Creating kernel function : " << kernelname << endl;
	TF1 * fker = new TF1(kernelname, GausFuncAdd, 0, m_nbins, 3 * m_nbins + 1);
	fker->SetParameters( par );

	return fker;
}

/*
//TF1 * TOTCalib::CreateKernelDensityFunction(TH1I * hd, double bandwidth) {
TF1 * TOTCalib::CreateKernelDensityFunction(int pix, vector<double> hist, double bandwidth) {

	//int nbins = hd->GetNbinsX();

	// name of the kernel sum
	TString kernelsum = "";
	vector<TString> functionsToErase;

	for(int i = 1 ; i <= m_nbins ; i++){

		TString funcname = "kernelbit_";
		//funcname += hd->GetName();
		funcname += pix;
		funcname += "_";
		funcname += i;

		// this will be erased from memory later
		functionsToErase.push_back( funcname );

		// Prepare kernel sum
		kernelsum += funcname;
		if(i < m_nbins) kernelsum += " + ";

		// Construct the gaussians to build the kernel density estimation
		TF1 * f = new TF1( funcname, "gaus", 0, m_nbins );
		//f->SetParameter(0, hd->GetBinContent(i) ); // constant
		//f->SetParameter(1, hd->GetBinLowEdge(i) ); // mean
		//f->SetParameter(2, bandwidth ); // sigma

		f->SetParameter(0, hist[i] ); // constant
		f->SetParameter(1, i ); // mean
		f->SetParameter(2, bandwidth ); // sigma

	}

	// Final kernel density
	TString kernelname = "kernel_";
	//kernelname += hd->GetName();
	kernelname += m_calhandler->GetSourcename();
	kernelname += "_";
	kernelname += pix;
	//cout << "Creating kernel function : " << kernelsum << endl;
	cout << "Creating kernel function : " << kernelname << endl;
	TF1 * fker = new TF1(kernelname, kernelsum, 0, m_nbins);

	// Erase the kernelbit functions.  Not needed anymore
	vector<TString>::iterator itr = functionsToErase.begin();
	for( ; itr != functionsToErase.end() ; itr++){
		TF1 * f = static_cast<TF1 *> ( gROOT->FindObject( *itr ) );
		delete f;
	}

	return fker;
}
 */

int TOTCalib::XYtoX(pair<int, int> pix, int dimX){
	return pix.second * dimX + pix.first;
}

pair<int, int> TOTCalib::XtoXY(int X, int dimX){
	return make_pair(X % dimX, X/dimX);
}

TH2I * TOTCalib::EntriesPlots(int i){

	// Histos to plot
	TFile * fc = m_inputfile[0];

	TH2I * she = (TH2I *) fc->Get("entriesSingle");
	TH2I * dhe = (TH2I *) fc->Get("entriesDouble");
	TH2I * the = (TH2I *) fc->Get("entriesTriple");
	TH2I * qhe = (TH2I *) fc->Get("entriesQuad");

	TH2I * avs = (TH2I *) fc->Get("averageTOTSingle");


	TH2I * histo = 0x0;

	if (i==0) { histo = she; }
	if (i==1) { histo = dhe; }
	if (i==2) { histo = the; }
	if (i==3) { histo = qhe; }
	if (i==4) { histo = avs; }

	if(!histo) {
		cout << "Couldn't retreive histo" << endl;
		return histo;
	}

	cout << "Draw " << histo->GetName() << " : " << histo << endl;

	// plus some info on top
	// average entries per pixel

	double average = 0.;
	int ntotalpixels = histo->GetNbinsX() * histo->GetNbinsY();
	for(int j = 1 ; j <= histo->GetNbinsX() ; j++) {
		for(int k = 1 ; k <= histo->GetNbinsY() ; k++) {
			average += histo->GetBinContent(j,k);
		}
	}
	average /= (double)ntotalpixels;
	cout << "Entries : " << TString::Format( "%ld", (Long_t)histo->GetEntries() ) << " | " << "average per pixel : " << average << endl;

	//s = TString::Format("Total entries = %d, Average entries per pixel = %.6f", (int)histo->GetEntries(), average );
	//TLatex * l1 = new TLatex();
	//l1->DrawLatex(0, 256, s);

	return histo;
}

//double TOTCalib::GetMean(int i){
//return m_calibhistos[i]->GetMean();
//}

TH1I * TOTCalib::GetHisto(int pix, TString extraName){

	// Create an histo from the vector<double> only if requested here
	TString name = m_calhandler->GetSourcename();
	name += "_";
	if(extraName.Length() != 0) {
		name += extraName;
		name += "_";
	}
	name += pix;
	TH1I * h = new TH1I(name, name, m_histoRebinning, 0, m_nbins);

	vector<double> hist = m_calibhistos[pix];
	vector<double>::iterator i = hist.begin();

	int cntr = 0;
	for ( ; i != hist.end() ; i++) {
		h->Fill(cntr, *i);
		cntr++;
	}

	return h;
}

TF1 * TOTCalib::GetKernelDensityFunction(int i){

	return TOTCalib::CreateKernelDensityFunction(i, m_calibhistos[i], m_bandwidth);

	//return m_kerneldensityfunctions[i];
}



int TOTCalib::SmallMin(int a, int b){
	if(a < b) return a;
	if(b < a) return b;
	return a;
}
int TOTCalib::GetSign(double slope){
	if(slope > 0) return 1;
	if(slope < -1) return -1;
	if(slope == 0) return 0;
	return 0;
}

double TOTCalib::FivePointsStencil(queue<double> q) {

	double sten = 0.;

	sten  =      q.front();     q.pop(); // get oldest element and remove
	sten -= 8. * q.front();     q.pop();

	q.pop();                             // middle element not needed in the stencil

	sten += 8. * q.front();     q.pop();
	sten -=      q.front();     q.pop();

	return sten/12.;
}


///
/// This is only good for the peaks identification.
/// The cut information here is meaningless
///
void TOTCalib::SimpleLinearReggresion(queue<int> q, double * slope){

	int N = (int) q.size();
	double x = 0;
	vector<pair<double, double> > q2;
	// build the x,y pairs
	while(!q.empty()) {
		// retrieve oldest element
		q2.push_back( make_pair( x , (double) q.front() ) );
		// and remove it
		q.pop();
		x += 1; // next x value ... arbitrary
	}

	vector<pair<double, double> >::iterator i = q2.begin();
	vector<double> xy, x2;
	double sumx =  0.;
	double sumy =  0.;
	double sumxy = 0.;
	double sumx2 =  0.;
	for( ; i != q2.end() ; i++) {
		xy.push_back( (*i).first * (*i).second );
		x2.push_back( (*i).first * (*i).first );
		sumx  += (*i).first;
		sumy  += (*i).second;
		sumxy += (*i).first * (*i).second;
		sumx2 += (*i).first * (*i).first;
	}

	// slope
	*slope = ( N*sumxy - sumx*sumy ) / ( N*sumx2 - sumx2*sumx2);

}

//TOTCalib::TOTCalib(TString fn, TString source) : fChain(0) {

//}

TOTCalib::TOTCalib(TString fn, TString source, int minpix, int maxpix, int maxtot, Long64_t nFrames, TOTCalib * oc, int method) : fChain(0) {

	// Copy certain properties from another calib and continue
	m_badPixelList = oc->GetBadPixelList();

	// Call the actual function doing the job
	SetupJob(fn, source, minpix, maxpix, maxtot, nFrames, method);

}

TOTCalib::TOTCalib(TString fn, TString source, int minpix, int maxpix, int maxtot, Long64_t nFrames, int method) : fChain(0) {

	// Call the actual function doing the job
	SetupJob(fn, source, minpix, maxpix, maxtot, nFrames, method);

}

void TOTCalib::SetupJob(TString fn, TString source, int minpix, int maxpix, int maxtot, Long64_t nFrames, int method) {

	m_TOTMultiplierFactor = 1.0;

	// Detector global THL
	m_thresholdEnergy = 3.6;      // 3.60 keV --> 1000e- : This can be set by the user
	m_thresholdEnergy_Err = 0.36; // 0.36 keV -->  100e- : This can be set by the user
	// Kernel density
	m_bandwidth = 5.0;
	// default verbose level
	m_verbose = __VER_INFO;

	cout << endl;
	cout << "----------------------------------------------------------------" << endl;

	// The number of bins per histogram will be set equal to the max tot in the distributions
	m_nbins = maxtot;
	m_histoRebinning = m_nbins; //m_nbins/8;
	m_nFrames = nFrames;
	m_method = method;

	// File
	TFile * f = new TFile(fn);
	m_inputfile.push_back( f );

    // Look at the MetaData
    TVectorF * md = static_cast<TVectorF *> ( f->Get("MetaData") );

    if ( md != 0x0 ) {
        TVectorF md1 = *md;
        __matrix_height = md1[0];
        __matrix_width = md1[1];
        __matrix_size = __matrix_height * __matrix_height;
        cout << "MetaData. Size of pad. Width = " << md1[0] << ", Height = " << md1[1] << endl;

        if ( minpix < 0 || maxpix > __matrix_size ) {
            cout << "[ERROR] The matrix has " << __matrix_size << "pixels. You have requested " << endl;
            cout << "[ERROR]   a range between " << minpix << " and " << maxpix << ". Giving up ... exit(1)" << endl;
            exit(1);
        }

        if ( maxpix - minpix < 0 ) {
            cout << "[ERROR] Request more pixels.  Zero have been scheduled for processing." << endl;
            exit(1);
        }
    } else {
        cout << "[WARNING] Can't verify metadata !!!  Running at your own risk !!!" << endl;
        __matrix_height = 256;
        __matrix_width = 256;
        __matrix_size = __matrix_height * __matrix_height;
        //cout << "Size of pad. Width = " << md1[0] << ", Height = " << md1[1] << endl;
    }


	// min and max pixels
	m_minpix = minpix;
	m_maxpix = maxpix;

	// TChain --> Tree
	TTree * tree = static_cast<TTree *> ( f->Get("TOTCalibrationPreparation") );

	// Calib handler
	m_calhandler = new CalibHandler( source.Data() );

	// Preparing density function vector
	m_kerneldensityfunctions = new TF1 * [__matrix_size];
	// initialize all pointers at 0
	for (int i = 0 ; i < __matrix_size ; i++) {
		m_kerneldensityfunctions[i] = 0x0;
    }

	// Creating histos in a std::vector
	/*
	TString hname;
	for ( int i = 0 ; i < __matrix_size ; i++ ) {
		hname =  "h_";
		hname += m_calhandler->GetSourcename().c_str();
		hname += "_";
		hname += i;
		//m_calibhistos[i] = new TH1F(hname, hname, 100, 0, 100);
		m_calibhistos.push_back( new TH1I(hname, hname, 100, 0, 100) );
	}
	 */
	// Initialize the representation of an histogram in vectors.
	// This is faster than having TH1 objects.  Only when needed
	// the conversion to TH1 will per performed.
	for ( int i = 0 ; i < __matrix_size ; i++ ) {
		m_calibhistos.push_back( vector<double> (m_nbins, 0.) );
	}

 	// Look at the DACs
    if ( md != 0x0 ) {
        TVectorF * v = static_cast<TVectorF *> ( f->Get("DACs") );
        TVectorF v1 = *v;
        cout << "DACs : ";
        for (int i = 0 ; i < v1.GetNoElements() - 1 ; i++) {
            cout << v1(i) << " " ;
        }
        cout << " | clock = " << v1(v1.GetNoElements() - 1) << " MHz \n";
    }
    

	Init(tree);

	// Some randon numbers
	// Seed the random number generator using the time
	time_t rawtime;
	time(&rawtime);
	m_ranseed_time = unsigned ( rawtime );
	cout << "[RAND] The random seed (localtime): " << m_ranseed_time << endl;
	m_rand1 = new TRandom1(m_ranseed_time);

}

TOTCalib::~TOTCalib() {

	// Close input files
	vector<TFile* >::iterator iff = m_inputfile.begin();
	for( ; iff != m_inputfile.end() ; iff++){
		(*iff)->Close();
		delete *iff;
	}

	cout << "Erasing objects ... " << m_calhandler->GetSourcename() << endl;

	// delete calibration handler
	delete m_calhandler;

	// Calib histograms is a std::vector<double>

	// and kernel functions
	for (int i = 0 ; i < __matrix_size ; i++) {
		if( m_kerneldensityfunctions[i] ) delete m_kerneldensityfunctions[i];
	}
	delete [] m_kerneldensityfunctions;

	//if (!fChain) return;
	//delete fChain->GetCurrentFile();

	if(m_gf_linear) delete m_gf_linear;
	if(m_gf_lowe) delete m_gf_lowe;

	vector<TF1 *>::iterator i = m_extra_tf1_to_erase.begin();
	for( ; i != m_extra_tf1_to_erase.end() ; i++) delete *i;

	// delete histos to store
	vector<TH1*>::iterator i2 = m_histosToSave.begin();
	for( ; i2 != m_histosToSave.end() ; i2++ ) {
		delete *i2;
	}

}

Int_t TOTCalib::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t TOTCalib::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void TOTCalib::Init(TTree * tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set object pointer
	m_SingleHitCoor = 0;
	m_SingleHitTOT = 0;
	m_DoubleHitCoor = 0;
	m_DoubleHitTOT = 0;
	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("timePerFrame", &timePerFrame, &b_timePerFrame);
	fChain->SetBranchAddress("m_SingleHitCoor", &m_SingleHitCoor, &b_m_SingleHitCoor);
	fChain->SetBranchAddress("m_SingleHitTOT", &m_SingleHitTOT, &b_m_SingleHitTOT);
	fChain->SetBranchAddress("m_DoubleHitCoor", &m_DoubleHitCoor, &b_m_DoubleHitCoor);
	fChain->SetBranchAddress("m_DoubleHitTOT", &m_DoubleHitTOT, &b_m_DoubleHitTOT);
	Notify();

}

Bool_t TOTCalib::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void TOTCalib::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}

Int_t TOTCalib::Cut(Long64_t /*entry*/)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}

void TOTCalib::DrawFullPixelCalib(int x, int y) {

	DrawFullPixelCalib( XYtoX( make_pair(x,y), __matrix_width) );

}

void TOTCalib::DrawFullPixelCalib(int pix) {

	TString cname = "pix_";
	cname += pix;

	pair<int, int> pix_xy = XtoXY(pix, __matrix_width);

	cout << endl;
	cout << "---------------------------------------------------------------" << endl;
	cout << "-- Requesting calibration summary for pixel " << pix << "(" << pix_xy.first << ","
			<< pix_xy.second << ")" <<  endl;

	// fits on data
	vector<double> fit_const = m_calibPointsConstants[pix];
	vector< pair<double, double> > fit_mean = m_calibPoints[pix];
	vector<double> fit_sigmas = m_calibPointsSigmas[pix];
	vector<double> fit_ia = m_calibPoints_ia[pix];
	vector<double> fit_ib = m_calibPoints_ib[pix];
	vector<double> fit_ic = m_calibPoints_ic[pix];
	vector<double> fit_it = m_calibPoints_it[pix];

	if( fit_const.empty() ) {
		cout << "[WARNING] No mix available for this pixel. Check if you are running on a subset of the matrix." << endl;
		return;
	}

	TString ctitle = "Data pixel ";
	ctitle += pix;
	ctitle += " (";
	ctitle += pix_xy.first;
	ctitle += ",";
	ctitle += pix_xy.second;
	ctitle += ")";

	// number of sources
	if(m_allSources.empty()) {
		cout << "This is not the Blender object. Use the Blender object to execute this function" << endl;
		return;
	}

	int nSources = (int) m_allSources.size();
	// number of division in data canvas
	int xdiv = TMath::Ceil( TMath::Sqrt(nSources) );
	int ydiv = TMath::Ceil( (double)nSources / (double)xdiv );

	// Check if everything fits. I need at least one extra box
	if( (xdiv*ydiv) == nSources) ydiv++;

	TCanvas * c1 = new TCanvas(cname, ctitle);
	c1->Divide(xdiv, ydiv);

	TH1 * h; TF1 * kf;

	bool legenddone = false;
	TLegend * leg1 = new TLegend(0.6, 0.6, 0.9, 0.9);
	leg1->SetFillColor(kWhite);
	leg1->SetBorderSize(1);

	vector<int> original_totval = m_calibTOTPeaks[pix];

	TF1 * gf; TF1 * gf_clone;

	int sour = 0;
	int orderCntr = 0;
	for ( ; sour < nSources ; sour++) {

		// Get in the canvas
		c1->cd(sour + 1);

		// Data
		h = m_allSources[sour]->GetHisto(pix, "summary");
		h->Draw("HIST");
		h->GetXaxis()->SetTitle("TOT");
		h->GetYaxis()->SetTitle("entries");

		// Kernel density
		kf = m_allSources[sour]->GetKernelDensityFunction(pix);
		kf->Draw("same");
		kf->SetLineColor(kBlack);
		kf->SetLineStyle(2);
		kf->SetLineWidth(1);

		kf->GetXaxis()->SetTitle("TOT");
		kf->GetYaxis()->SetTitle("entries");

		// Info about the source
		TString sourceName = m_allSources[sour]->GetCalibHandler()->GetSourcename();

		if ( ! legenddone ) {
			leg1->AddEntry(kf, "k.d.f. (envelope)", "L");
			leg1->AddEntry(h, "Calibration data", "L");
		}

		TLatex * l1 = new TLatex();
		l1->DrawLatex(m_nbins/2, h->GetMaximum() * 0.5 , sourceName);

		// Get list of peaks identified
		map<int, double> calibPoints = m_allSources[sour]->GetCalibHandler()->GetCalibPoints();
		vector<double> peaks = (m_allSources[sour]->GetMaxPeaksIdentified())[pix];

		// Points used in the fit
		vector< pair<double, double> > calibFitPoints = GetCalibPoints(pix);

		int nCalibPoints = (int)calibPoints.size();
		double pos_offset = 1.;

		// Check if it's possible to draw
		if ( (int) peaks.size() < nCalibPoints ) {
			cout << "Not enough peaks were identified for this pixel !" << endl;
			return;
		}

		// loop over calip points per source
		for(int p = 0 ; p < nCalibPoints ; p++) {

			if(calibPoints[p]<0){ // skips energies < 0 (ignored points)
				orderCntr++;
				continue;
			}

			// Points used in the fit
			//cout << "---> " << calibFitPoints[orderCntr].first << " , " << calibFitPoints[orderCntr].second << endl;

			//c1->cd(sour + 1);

			// Fit on data
			TString funcName = "f_histofit_";
			funcName += pix;
			funcName += "_";
			funcName += p;
			gf = FittingFunctionSelector( calibPoints[p], m_allSources[sour], p );
            
            gf->FixParameter( 0, fit_const[orderCntr] );
            gf->FixParameter( 2, fit_sigmas[orderCntr] );
            //gf->SetNpx(1000);

            if( TString(gf->GetName()).Contains("gf_lowe") ) {
                gf->FixParameter( 1, fit_mean[orderCntr].first ); // fixed energy
                gf->FixParameter( 3, fit_ia[orderCntr] );
                gf->FixParameter( 4, fit_ib[orderCntr] );
                gf->FixParameter( 5, fit_ic[orderCntr] );
                gf->FixParameter( 6, fit_it[orderCntr] );

            }else{
                gf->SetParameter( 1, fit_mean[orderCntr].second ); // mean tot from fit
            }
            
            
			cout<<"--------- Source: " << sourceName<< " ---------" <<endl << "The fitting function for this source is : " << gf->GetName() << endl;
			//gf->Dump();
			gf_clone = static_cast<TF1 * > ( gf->Clone( funcName ) );
			//gf_clone->Dump();

			cout << "Point : " << calibPoints[p] << " | order : " << orderCntr << " [ ";
            
            if( TString(gf->GetName()).Contains("gf_lowe") ) {

                cout << "Constant: "<< fit_const[orderCntr] << ", " << "Mean: " <<fit_mean[orderCntr].first << ", " <<"Sigma: "<< fit_sigmas[orderCntr] << ", ";
                cout <<"a: "<< fit_ia[orderCntr] << ", " <<"b: "<< fit_ib[orderCntr] << ", " <<"c: "<< fit_ic[orderCntr] << ", " <<"t: "<< fit_it[orderCntr];

            }else{

                cout << "Constant: "<< fit_const[orderCntr] << ", " << "Mean: " <<fit_mean[orderCntr].second << ", " <<"Sigma: "<< fit_sigmas[orderCntr] << ", ";

            }
			cout << " ]" << endl;
            
			gf_clone->SetLineColor(kRed);
			gf_clone->Draw("same");                     // The drawing is not really happening in this scope.  That's why we need clones
			m_extra_tf1_to_erase.push_back( gf_clone ); // schedule to erase at the end

			if ( ! legenddone ) {
				leg1->AddEntry(gf_clone, "Fit on spectrum", "L");
			}

			//TString peak = TString::Format( "%.1f  -> %.3f keV", calibFitPoints[orderCntr].second, calibPoints[p] );
			TString peak = TString::Format( "%.1f  -> %.3f keV", fit_mean[orderCntr].second, calibPoints[p] );
			l1->DrawLatex( peaks[p], h->GetMaximum() * pos_offset , peak );

			pos_offset -= 0.2;

			orderCntr++;
		}

		if ( ! legenddone ) {
			leg1->Draw();
			legenddone = true;
		}

	}

	// If there is still room in the same pad
	if(sour <= nSources) {

		// Get in the canvas
		c1->cd(sour + 1);
		c1->GetPad(sour + 1)->SetGridx();
		c1->GetPad(sour + 1)->SetGridy();

		// Graph
		TGraphErrors * g = GetCalibGraph(pix);
		g->Draw("A*");

		TF1 * s = GetSurrogateFunction(pix);
		s->Draw("same");
		s->GetXaxis()->SetTitle("Energy (keV)");
		s->GetYaxis()->SetTitle("TOT");
		s->SetNpx(1000);

		// Set the maximum in Y
		double * ym = g->GetY();
		vector<double> ym_v;
		for(int i = 0 ; i < g->GetN() ; i++) ym_v.push_back(ym[i]);
		double maxel_y = *max_element( ym_v.begin(), ym_v.end() );
		s->GetYaxis()->SetRangeUser( 0., maxel_y*1.1 );
		g->GetYaxis()->SetRangeUser( 0., maxel_y*1.1 );


		// Set the maximum in X
		double * xm = g->GetX();
		vector<double> xm_v;
		for(int i = 0 ; i < g->GetN() ; i++) xm_v.push_back(xm[i]);
		double maxel_x = *max_element( xm_v.begin(), xm_v.end() );
		//s->GetXaxis()->SetRangeUser( 0., maxel_x*1.1 );

		TLatex * l2 = new TLatex();
		TString parS;
		int npar = s->GetNpar();
		for(int i = 0 ; i < npar ; i++) {
			parS = "";
			if ( i == 0 ) parS = "a = ";
			if ( i == 1 ) parS = "b = ";
			if ( i == 2 ) parS = "c = ";
			if ( i == 3 ) parS = "t = ";
			parS += TString::Format("%.2f"/* +/- %.2f*/, s->GetParameter(i)/*, s->GetParError(i) */); 	// error on fit parameters is not computed
			l2->DrawLatex(maxel_x/2, maxel_y * (1 - (i/10.)), parS); 									// maybe add soon?
		}

		//
		sour++;
	}

	c1->Update();

}

TGraphErrors * TOTCalib::GetCalibGraph(int pix){

	vector<pair<double, double> > points = m_calibPoints[pix];
	vector<double> err = m_calibPointsSigmas[pix];

	TGraphErrors * g = new TGraphErrors(); // we don't know how many points are not ignored

	vector<pair<double, double> >::iterator i = points.begin();
	vector<double>::iterator ie = err.begin();

	int cntr = 0;
	for( ; i != points.end() ; i++,ie++) {
        
        //cout<<(*i).first<<" "<<(*i).second<<" "<<*ie<<endl;
		if( (*i).first > 0){
			g->SetPoint(cntr, (*i).first, (*i).second );
			g->SetPointError(cntr, 0.,  *ie );
			cntr++;

		}
	}

	g->SetPoint( cntr, m_thresholdEnergy, 0.0 );
	g->SetPointError(cntr, m_thresholdEnergy_Err, 0.0 );

	return g;
}

CalibHandler::CalibHandler (string source) {

	m_sourceName = source;
	cout << "Source is " << m_sourceName;

	if( ! TString(source).CompareTo( "Am241", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Am-241", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Am_241", TString::kIgnoreCase )
	) {

		// Am-241
		// Energy, Intensity
		// XR l	    13.9	    37 % 3
		//          59.5409 1   35.9 % 4
		//m_calibPoints[0] = 13.9;    // Am-241 13.9    keV
		//m_calibPoints[1] = 26.3446; // Am-241 26.3446    keV
		//m_calibPoints[2] = 59.5409; // Am-241 59.5409 keV

		//m_calibPoints[0] = 13.9;      // Am-241 13.9    keV
		m_calibPoints[0] = 9.0;      // Closed source in the lab !!

		//m_calibPoints[1] = 26.3446;   // Am-241 26.3446    keV
		m_calibPoints[1] = 59.5409;   // Am-241 59.5409 keV

		m_calibPointsRegion[0] = __linear_reg; // this point can be used to define the linear region
		m_calibPointsRegion[1] = __linear_reg; // this point can be used to define de linear region
		//m_calibPointsRegion[2] = __linear_reg; // this point can be used to define de linear region

	} else if (
			! TString(source).CompareTo( "Fe55", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Fe-55", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Fe_55", TString::kIgnoreCase )
	) {

		// Fe-55
		// Energy, Intensity
		// XR kα2	     5.888	      8.2 % 4 	  4.85E-4 21
		// XR kα1	     5.899	     16.2 % 7 	  9.6E-4 4
		m_calibPoints[0] = 5.899;  // Fe-55 5.899 KeV
		m_calibPointsRegion[0] = __lowenergy_reg;
		//m_calibPointsRegion[0] = __linear_reg;

	} else if ( ! TString(source).CompareTo( "In_fluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "In-fluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Influo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "In", TString::kIgnoreCase )
	) {

		// Indium
		// Fluorescence  24.2 keV
		m_calibPoints[0] = 24.2; //
		m_calibPointsRegion[0] = __linear_reg;

	} else if ( ! TString(source).CompareTo( "Cd_fluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Cd-fluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Cdfluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Cd", TString::kIgnoreCase )
	) {

		// Indium
		// Fluorescence  24.2 keV
		m_calibPoints[0] = 23.1; //
		m_calibPointsRegion[0] = __linear_reg;

	} else if ( ! TString(source).CompareTo( "Cd_fluo_TPXGaAs", TString::kIgnoreCase )
    ) {
        
        m_calibPoints[0] = 10.5; // K-alpha peak of As
        m_calibPointsRegion[0] = __linear_reg;
        m_calibPoints[1] = 23.1; // K-alpha peak of Cd 
        m_calibPointsRegion[1] = __linear_reg;

    } else if ( ! TString(source).CompareTo( "Cu_fluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Cu-fluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Cufluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Cu", TString::kIgnoreCase )
	) {

		// Fluorescence  8.05 keV
		m_calibPoints[0] = 8.05; //
		m_calibPointsRegion[0] = __lowenergy_reg;

	} else if ( ! TString(source).CompareTo( "Am241CutLow_Plus_SnFluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Am241CutLow-Plus-SnFluo", TString::kIgnoreCase )
	) {

		// Fluorescence  Sn
		m_calibPoints[0] = 25.27;
		m_calibPointsRegion[0] = __linear_reg;

		m_calibPoints[1] = 59.5409;   // Am-241 59.5409 keV
		m_calibPointsRegion[1] = __linear_reg; // this point can be used to define de linear region

	} else if ( ! TString(source).CompareTo( "Cd109", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Cd-109", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Cd_109", TString::kIgnoreCase )
	) {

		m_calibPoints[0] = 23.0; //
		m_calibPointsRegion[0] = __linear_reg;

	} else if ( ! TString(source).CompareTo( "Am241In", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Am241_In", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Am241-In", TString::kIgnoreCase 	)
	) {
		// X-ray Am241
		m_calibPoints[0] = 13.9;	// kev
		m_calibPointsRegion[0] = __linear_reg;

		// Fluorescence In
		m_calibPoints[1] = 24.2;	// kev
		m_calibPointsRegion[1] = __linear_reg;

		// X-ray Am241
		m_calibPoints[2] = 59.5409;   // kev
		m_calibPointsRegion[2] = __linear_reg; // this point can be used to define de linear region

	} else if ( ! TString(source).CompareTo( "ZrFluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Zr-Fluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Zr_Fluo", TString::kIgnoreCase )
	) {
		// Fluorescence Zr
		m_calibPoints[0] = 15.78; // kev
		m_calibPointsRegion[0] = __linear_reg;

    } else if ( ! TString(source).CompareTo( "ZnFluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Zn-Fluo", TString::kIgnoreCase )
			||
			! TString(source).CompareTo( "Zn_Fluo", TString::kIgnoreCase )
	) {
    	// Fluorescence Zn
		m_calibPoints[0] = 8.64; // kev
		m_calibPointsRegion[0] = __linear_reg; // could be __lownenergy_reg in some cases
    }

	cout << " and contains " << m_calibPoints.size() << " calibration point(s)" << endl;

}

void printProgBar( int percent ){

	string bar;

	for(int i = 0; i < 50; i++){
		if( i < (percent/2)){
			bar.replace(i,1,"=");
		}else if( i == (percent/2)){
			bar.replace(i,1,">");
		}else{
			bar.replace(i,1," ");
		}
	}

	cout<< "\r" "[" << bar << "] ";
	cout.width( 3 );
	cout<< percent << "%     " << std::flush;

}

double surrogatefunc_calib(double * x, double * par) {

	// independent var
	double xx = x[0];

	// pars
	double a = par[0];
	double b = par[1];
	double c = par[2];
	double t = par[3];

	double func = a * xx + b;
	func -= ( c / ( xx - t) );

	return func;
}


double fitfunc_lowen(double * x, double * par) {

    // Function proposed in J. Jakubek / Nuclear Instruments and Methods in Physics Research A 633 (2011) S262–S266

    // independent var
    Double_t xx = x[0]; // TOT

    // parameters for gaussian
    Double_t gconst = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];

    // surrogate ^ -1
    Double_t a = par[3];
    Double_t b = par[4];
    Double_t c = par[5];
    Double_t t = par[6];

    // Inverse of the surrogate ( = energy in keV)
    Double_t surrogate_inverse, delta, num1, num2, denum;
    delta = TMath::Power((b-a*t-xx),2) - 4 * a * (xx*t-c-b*t);
    num1 = a*t + xx - b;
    num2 = TMath::Sqrt(delta) ;
    denum = 2 * a;
    surrogate_inverse = (num1 +  num2) / denum;

    // Gaussian of the inversed surrogate
    Double_t arg = (surrogate_inverse - mean)/sigma;
    Double_t func = gconst * TMath::Exp(-0.5*arg*arg);

    return func;
}

void TOTCalib::GetInputStats() {

		/////////////////////////////////////////////////////////////////////////
		// stats
		TString source = this->GetCalibHandler()->GetSourcename();
		TCanvas * c1 = new TCanvas("Stats_"+source, "Stats_"+source);
		TH2I * entrieshisto = this->EntriesPlots(0);
		c1->cd();
		entrieshisto->Draw("colz");

}
