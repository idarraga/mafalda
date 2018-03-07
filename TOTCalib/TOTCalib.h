/**
 *  Author: John Idarraga <idarraga@cern.ch>
 *  TOT calibration for a Timepix device
 */

#ifndef TOTCalib_h
#define TOTCalib_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TLatex.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2I.h>
#include <TF1.h>
#include <TVectorF.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TRandom1.h>
#include "Math/MultiRootFinder.h"
#include "Math/WrappedMultiTF1.h"

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <map>
#include <algorithm>
#include <unordered_map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;


#define __npars_surrogate			    4
#define __npars_lowe_fitfunc            7
#define __min_tmathprobtest_val      0.90
#define __max_fit_tries               100
#define __fit_pars_randomization_max 1000
#define __fit_pars_randomization_min   20

#define __fraction_of_height_range_id 0.5

// Parameter hints for low energy function fit
// FIXME: this should be done without user input
// Could fit the global spectrum with high fit tries and choose the best found parameters
#define __lowen_sigma_hint 2.0
#define __lowen_para_hint 3.0 // useless if at least 2 fits in the linear region succeeds
#define __lowen_parb_hint 80.0 // useless if at least 2 fits in the linear region succeeds
#define __lowen_parc_hint 200.0
#define __lowen_part_hint 2.0
#define __lowen_par_fraction_random 0.4 // range for randomization around the given hint (i.e. put 0.0 to set the hint without randomization) 

// Prototypes
double GausFuncAdd(double * x, double * par);
void printProgBar( int );
#define __fitfunc_lowen_npars  7
double fitfunc_lowen(double * x, double * par);
double fitfunc_lowen2(double * x, double * par);
double fitfunc_lowen_ZERO(double * x, double * par);
double surrogatefunc_calib(double * x, double * par);
double surrogatefunc_calib_ZERO(double * x, double * par);
pair<double,double> Calculate_ab_From_ct_e1s1_e2s2(double,double,double,double,double,double);


// Calibration Handler for each source
class CalibHandler {

public:

	CalibHandler(string);
	~CalibHandler(){};
	string GetSourcename () { return m_sourceName; };
	map<int, double> GetCalibPoints(){ return m_calibPoints; };
	map<int, int> GetCalibPointsRegion(){ return m_calibPointsRegion; };
	double GetOneEnergyMatch(int i){return m_calibPoints[i]; };

	void CalibIgnorePoint(int i){ m_calibPoints[i]*= -1; };
	void CalibSetPointRegion(int i, int reg){ m_calibPointsRegion[i] = reg; }

	void Set_m_calibPoints( map< int, double> p){ m_calibPoints = p;};
	void Set_m_calibPointsRegion( map< int, int> r){ m_calibPointsRegion = r;};

	enum {
		__linear_reg = 0,
		__lowenergy_reg,
		__other_reg
	};

private:


	// TOT, Energy for the particular source or fluorescence data
	map<int, double> m_calibPoints; // key: peak number in the source's histogram, value: energy of the peak (as defined in CalibHandler constructor) 
	map<int, int> m_calibPointsRegion;
	string m_sourceName;
};


class TOTCalib {

public :

	vector<TFile *> m_inputfile; // input files

	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain

	// Declaration of leaf types
	Double_t        timePerFrame;
	vector<int>     *m_SingleHitCoor;
	vector<int>     *m_SingleHitTOT;
	vector<int>     *m_DoubleHitCoor;
	vector<int>     *m_DoubleHitTOT;

	// List of branches
	TBranch        *b_timePerFrame;   //!
	TBranch        *b_m_SingleHitCoor;   //!
	TBranch        *b_m_SingleHitTOT;   //!
	TBranch        *b_m_DoubleHitCoor;   //!
	TBranch        *b_m_DoubleHitTOT;   //!

	//TTree * m_tree;

    // Methods for peak selection/fitting
    enum {
        __peakStandard = 0,
        __peakLowStats 	// Fit range and gaussian mean are fixed, better peak selection
    };

    // Methods for calibration/surrogate fitting
    enum {
    	__calibStandard = 0,
    	__calibJakubek, // Method proposed in J. Jakubek / Nuclear Instruments and Methods in Physics Research A 633 (2011) S262–S266
    	__calibJakubekAlt // Alternate method, use new parametrization of the surrogate
	};

	TOTCalib();   
	TOTCalib(TString, TString, int minpix, int maxpix, int maxtot, Long64_t nFrames, int peakMethod = __peakStandard);
	TOTCalib(TString, TString, int minpix, int maxpix, int maxtot, Long64_t nFrames, TOTCalib *, int peakMethod = __peakStandard);

	void SetupJob(TString, TString, int minpix, int maxpix, int maxtot, Long64_t nFrames, int peakMethod = __peakStandard);

	virtual ~TOTCalib();
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual void     Loop();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
	double FivePointsStencil(queue<double>);
	void SimpleLinearReggresion(queue<int> q, double * slope);
	int SmallMin(int a, int b);
	void SetTOTMultiplierFactor(double f) { m_TOTMultiplierFactor = f;};
	void SetHistoRebinning(int v){ m_histoRebinning = v; };
	void GetInputStats ();

	// points to store
	struct store {

		vector< pair<double, double> > pointsSave_E_TOTfit;
		vector< double > pointsSaveSigmas;
		vector< double > pointsSaveConstants;
		vector< double > pointsSave_ia;
		vector< double > pointsSave_ib;
		vector< double > pointsSave_ic;
		vector< double > pointsSave_it;
        vector< double > pointsSave_chi2ndf;
		vector< int > calibTOTPeaks;
		vector<pair<double,double> > linearpairs; //(E, TOT)
		vector<int> peakFitStatus;
	};

	TH1I * GetHisto(int, TString extra = "");
	TF1 * GetKernelDensityFunction(int);
	int GetCriticalPoints(int i , vector<double> &, vector<double> &);
    int GetCriticalPoints2(int pixID , vector<double> &, vector<double> &);
    int GetSign(double slope);
	unordered_map<int, vector<double> > GetMaxPeaksIdentified(){ return m_critPointsMax;}
    unordered_map<int, vector<double> > GetMaxPeaksIdentified_amplitude(){ return m_critPointsMax_amplitude; }	
    vector <vector<double> > GetMaxPeaksIdentified_vec(){ return m_critPointsMax_vec;}
    vector < vector<double> > GetMaxPeaksIdentified_amplitude_vec(){ return m_critPointsMax_amplitude_vec; }	
    CalibHandler * GetCalibHandler(){ return m_calhandler; }
	vector< pair<double, double> > GetCalibPoints(int pix){ return m_calibPoints_E_TOTfit[pix]; }
	TF1 * GetSurrogateFunction(int);
	TGraphErrors * GetCalibGraph(int pix);
	void RandomFitParametersSurrogate (TF1 * f, double, double);
	bool PixelInLowActivityList(int pix);
	void PushToLowActivityList(int pix);
	void PushToBadPixelList(int pix);
	void PushToBadPixelListRangeInclusive(int pixi, int pixf);
	void PushToBadPixelList(int x, int y);
	bool PixelInBadPixelList(int pix);
	vector<int> GetBadPixelList() { return m_badPixelList; };
    
    void Blender(TOTCalib * , TOTCalib *, TOTCalib *, TOTCalib *, TString, int = 0);
	void Blender(TOTCalib * , TOTCalib *, TOTCalib *, TString, int = 0);
	void Blender(TOTCalib * , TOTCalib *, TString, int = 0);
	void Blender(TOTCalib * s2, TString outputName, int = 0);
	void Blender(TString, int = 0);
    void Blender2(TOTCalib * , TOTCalib *, TString, int = 0);

    void SavePixelResolution(TString = "", TString = "", TString = "", TString = "");
    void GetCoeffFromFiles(double *, double *, double *, double *, const char *, const char *, const char *, const char *, int, int);
    void GetCuadraticSolutions(pair<int,int>, int, double &, double &, double &, double *, double *, double *, double *);
    double GetE(pair<int,int>, int tot, double *, double *, double *, double *);
    TH1I *GetHistoCalibrated(int, TString, double, double *, double *, double *, double *);
    
	void ProcessOneSource(TOTCalib * s, store * sto, TGraphErrors * g, int pix, int & cntr);
    void ProcessOneSource2_gaussian(TOTCalib *, store *, int);
    void ProcessOneSource2_lowen(TOTCalib *, store *, int, double &, double &, double &, double &);	
    void ReorderSources();

	//int PeakFit(TOTCalib *, int, int, TF1 *, TH1 *, store *);
    int PeakFit(TOTCalib *, int, int, TF1 *, TH1 *, store *, double = 0.0); 
    int PeakFit2_gaussian(TOTCalib *, int, int, TF1 *, TH1 *, store *);    
    int PeakFit2_lowen(TOTCalib *, int, int, TF1 *, TH1 *, store *, double = 0.0);
	TF1 * FittingFunctionSelector(double, TOTCalib *, int);
	void GetLinearFit(double & a, double & b, vector< pair<double,double> > p);
	void RandomFitParameters(TF1 * f, TH1 * h, int tot, TOTCalib *);

	vector<pair<double, double> > Extract_E_TOT_Points(int, TOTCalib * );
    vector<pair<double, double> > Extract_E_TOT_Points2(int, TOTCalib * );	
    int GetNumberOf_E_TOT_Points (TOTCalib * s);
	int GetNumberOf_E_TOT_Points_Positive (TOTCalib * s);

	double DerivativeFivePointsStencil(TF1 *, double, double);
	double VectorSum(vector<double>);

	int XYtoX(pair<int, int> pix, int dimX);
	pair<int, int> XtoXY(int X, int dimX);
	TH2I * EntriesPlots(int);
	//TF1 * CreateKernelDensityFunction(TH1I *, double);
	TF1 * CreateKernelDensityFunction(int, vector<double>, double);

	void SetGlobalThresholdEnergyAndErr(double th, double th_err) { m_thresholdEnergy = th; m_thresholdEnergy_Err = th_err; };
	pair<double, double> GetGlobalThresholdEnergyAndErr() {return make_pair(m_thresholdEnergy, m_thresholdEnergy_Err);};

	void FillHisto(int, int);

	void Finalize();
	TH1F * CreateParameterHistogram(vector<double>, TString);

	void DrawFullPixelCalib(int);
	void DrawFullPixelCalib(int, int);
    void DrawFullPixelCalib2(int);

	void SetVerboseLevel(int v) {m_verbose = v;};
	enum {
		__VER_DEBUG_LOOP = 0,
		__VER_DEBUG,
		__VER_INFO,
		__VER_QUIET
	};

	double GetKernelBandWidth(){ return m_bandwidth; };
	double SetKernelBandWidth(double b) {return m_bandwidth = b; };

	int GetNBins() { return m_nbins; };

	void IgnorePoint(int i){ m_calhandler->CalibIgnorePoint(i); };
	void SetPointRegion(int i, int reg){ m_calhandler->CalibSetPointRegion(i,reg); };

	int GetPeakMethod() {return m_peakMethod; };
	int GetCalibMethod() { return m_calMethod; };
	void SetPeakMethod(int m) {m_peakMethod = m; };
	void SetCalibMethod(int m) {m_calMethod = m; };

	vector<vector<double> > Get_m_histo(){return m_calibhistos;}; // each pixel has its own histogram
	void SetGlobalHisto(vector<double> v){m_globalhisto=v;}; // every histogram are merged in this vector
	vector<double> GetGlobalHisto(){return m_globalhisto;};

	void SetGlobalCriticalPoints(vector<double> max, vector<double> min){ 
		m_global_max=max;
		m_global_min=min;	};
	vector<double> GetGlobalMaximumPoints(){return m_global_max;};
	vector<double> LowStatsPeakSelection( vector<double> peaks, unsigned int s, TOTCalib * source, double bandwidth);

	void CreateGlobalKernelAndGetCriticalPoints(); 

	enum surr_status{ 	//surrogate function status
		__no_data = -1, // no source has sufficient data, surrogate cannot be built
		__partial_data, // one or many sources has insufficient data, but surrogate was built anyway (possible bad fit)
		__good_data,  	// every peak of every source was identified and the results are probably reliable
	};

	void DumpCalibParametersFromSavedFile(	map<int, vector< pair<double, double> > > points,  
											map<int, vector<double> > sig, 
											map<int, vector<double> > consts,
											map<int, vector<double> > ia,
											map<int, vector<double> > ib,
											map<int, vector<double> > ic,
											map<int, vector<double> > it,
											map<int, vector<double> > param,
											map<int, int> status,
											int calibMethod){
		m_calibPoints_E_TOTfit = points;
		m_calibPointsSigmas = sig;
		m_calibPointsConstants = consts;
		m_calibPoints_ia = ia;
		m_calibPoints_ib = ib;
		m_calibPoints_ic = ic;
		m_calibPoints_it = it;
		m_calibSurrogateConstants = param;
		m_surrogateStatus = status;
		m_maxpix = status.rbegin()->first;
		m_minpix = status.begin()->first;
		m_calMethod = calibMethod;

	};

	void DumpSpectrumVectorFromSavedFile(	vector< vector<double> > spectrum){
		m_calibhistos = spectrum;
		m_nbins = (spectrum[0]).size();
		m_histoRebinning = m_nbins;
	};

	void CreateCalibHandlerFromSavedFile( TString s){
		m_calhandler = new CalibHandler(s.Data());
	};

	void DumpSourceInfoFromSavedFile(map<int, double> p, map<int, int> r, unordered_map<int, vector<double> > max, int peakMethod){
		m_calhandler->Set_m_calibPoints(p);
		m_calhandler->Set_m_calibPointsRegion(r);
		m_critPointsMax = max;
		m_peakMethod = peakMethod;
	};

	void AddSingleSourceFromSavedFile( TOTCalib * s){
		m_allSources.push_back(s);
	};

	void CreateGausAndLowEnergFitFunction(){
		vector<TOTCalib *>::iterator iir = m_allSources.begin();
		double maxrange = 0.;
		for ( ; iir != m_allSources.end() ; iir++ ) {
			if( (*iir)->GetNBins() > maxrange ) maxrange = (double) ( (*iir)->GetNBins() );
		}

		m_gf_linear = new TF1("gf_linear", "gaus(0)", 0., maxrange);
		m_gf_linear->SetParameters(1, 1, m_bandwidth);
	
		m_gf_lowe = new TF1("gf_lowe", fitfunc_lowen, 0., maxrange, __fitfunc_lowen_npars);
		m_gf_lowe->SetParameters(1, 1, m_bandwidth, 1,1,1,1);

		m_gf_lowe_ZERO = new TF1("gf_lowe_ZERO", fitfunc_lowen_ZERO, 0., maxrange, __fitfunc_lowen_npars);
		m_gf_lowe_ZERO->SetParameters(1, 1, m_bandwidth, 1,1,1,1);
	}

	int GetMatrixWidth(){return __matrix_width;};
	int GetMatrixHeight(){return __matrix_height;};
	int GetMatrixSize(){return __matrix_size;};
	map<int, int> GetSurrogateStatusMap(){return m_surrogateStatus;};
	map<int, vector<double> > GetSurrogateParamMap(){return m_calibSurrogateConstants;};
	vector<TOTCalib *> GetSourcesVector(){return m_allSources;};

	map<int, vector<double>> GetMapCalibPointsConstants(){return m_calibPointsConstants;};
	map<int, vector< pair<double, double> > > GetMapCalibPoints(){return m_calibPoints_E_TOTfit;};
	map<int, vector<double> > GetMapCalibSigmas(){return m_calibPointsSigmas;};
	map<int, vector<double> > GetMapCalibIA(){return m_calibPoints_ia;};
	map<int, vector<double> > GetMapCalibIB(){return m_calibPoints_ib;};
	map<int, vector<double> > GetMapCalibIC(){return m_calibPoints_ic;};
	map<int, vector<double> > GetMapCalibIT(){return m_calibPoints_it;};

	void ParametersEstimation(int);
	bool GetParametersEstimationStatus(){return globalEstimationSuccess;};
	void SetThresholdBound(double e0){m_e0_bound = e0;};

	vector<double> GetParametersEstimation(){ 
		vector <double> v;   v.push_back(m_glob_const.first);   v.push_back(m_glob_sig.first);   v.push_back(m_glob_a.first);
		v.push_back(m_glob_b.first);   v.push_back(m_glob_c.first);   v.push_back(m_glob_t.first); v.push_back(m_glob_e0.first);    return v; };

	vector<double> GetParametersEstimationErrors(){ 
		vector <double> v;   v.push_back(m_glob_const.second);   v.push_back(m_glob_sig.second);   v.push_back(m_glob_a.second);
		v.push_back(m_glob_b.second);   v.push_back(m_glob_c.second);   v.push_back(m_glob_t.second); v.push_back(m_glob_e0.second);    return v; };

	void DumpParametersEstimation(vector< double> v, vector<double> verr){
		m_glob_const.first = v[0]; m_glob_const.second = verr[0];
		m_glob_sig.first = v[1]; m_glob_sig.second = verr[1];
		m_glob_a.first = v[2]; m_glob_a.second = verr[2];
		m_glob_b.first = v[3]; m_glob_b.second = verr[3];
		m_glob_c.first = v[4]; m_glob_c.second = verr[4];
		m_glob_t.first = v[5]; m_glob_t.second = verr[5];
		m_glob_e0.first = v[6]; m_glob_e0.second = verr[6];
		globalEstimationSuccess = true; };

    void WriteCalibToAsciiFiles(TString);
    void Choose2LinearPeaks(double p1, double p2){m_linearPeak1=p1;m_linearPeak2=p2;return;}
    void SetLowEnergyPeak(double p1){m_lowenPeak=p1;return;}    
    void DrawFullPixelCalib_coeff_histos();
    void SetLowEnFit_Params(double cons, double sig, double c, double t){m_lowen_fitParams[0]=cons; m_lowen_fitParams[1]=sig; m_lowen_fitParams[2]=c; m_lowen_fitParams[3]=t;return;}
    Double_t* GetLowEnFit_Params(){return m_lowen_fitParams;}
    
private:
	//////////////////////////////////////////////////////////////////
	// histograms and info for all pixels;
	//vector<TH1I *> m_calibhistos;
	vector<vector<double> > m_calibhistos;
	// key = pixel, values = points selected for calibration
	map<int, vector< pair<double, double> > > m_calibPoints_E_TOTfit;
	map<int, vector<int> > m_calibTOTPeaks;
	map<int, vector<double> > m_calibPointsSigmas;
	map<int, vector<double> > m_calibPointsConstants;
	map<int, vector<double> > m_calibPoints_ia;
	map<int, vector<double> > m_calibPoints_ib;
	map<int, vector<double> > m_calibPoints_ic;
	map<int, vector<double> > m_calibPoints_it;
	map<int, int > m_surrogateStatus;

	int m_nbins;
	int m_histoRebinning;
	Long64_t m_nFrames;
	int m_minpix, m_maxpix;
	TF1 ** m_kerneldensityfunctions;

	// key = pixel, val = vector of critical points
	unordered_map<int, vector<double> > m_critPointsMax;
    unordered_map<int, vector<double> > m_critPointsMax_amplitude; 
    // try if faster with vector
	vector < vector<double> > m_critPointsMax_vec;
    vector < vector<double> > m_critPointsMax_amplitude_vec;  
	// key = pixel, val = vector of critical points
	map<int, vector<double> > m_critPointsMin;
	// This map contains the final: a,b,c,t parameters of the calibration
	map<int, vector<double> > m_calibSurrogateConstants;
	// This map contains other final properties of the calibration like the prob TMath::Prob()
	map<int, vector<double> > m_calibSurrogateProperties;
	// This map contains only the linear pars of the calibration: a,b
	map<int, pair<double, double> > m_calibLinearMap;

	CalibHandler * m_calhandler;

	// Bandwidth for the Kernel density method
	double m_bandwidth;

	// TOT multiplier factor.  In case we have different clocks between sources
	double m_TOTMultiplierFactor;

	// Parameters histos
	vector<double> m_par_a_v;
	vector<double> m_par_b_v;
	vector<double> m_par_c_v;
	vector<double> m_par_t_v;
	vector<double> m_surr_prob_v;

//	TH1F * m_par_a;
//	TH1F * m_par_b;
//	TH1F * m_par_c;
//	TH1F * m_par_t;
//	TH1F * m_surr_prob;

	// Threshold information
	double m_thresholdEnergy;
	double m_thresholdEnergy_Err;

	// Low activity list
	vector<int> m_lowActivityList;
	vector<int> m_badPixelList;

	// output ROOT file
	TFile * m_output_root;
	vector<TH1 *> m_histosToSave;

	// If this is the Blender object store the other sources pointers
	vector<TOTCalib *> m_allSources;

	// Verbose level
	int m_verbose;

	// Functions for fitting peaks in the data
	TF1 * m_gf_linear;
	TF1 * m_gf_lowe;
	TF1 * m_gf_lowe_ZERO;
	vector<TF1 * > m_extra_tf1_to_erase;

	// Randon generators
	unsigned m_ranseed_time;
	TRandom1 * m_rand1;

	// Matrix size
	int __matrix_size;
	int __matrix_height;
	int __matrix_width;

	int m_peakMethod; // peak fit method
	int m_calMethod ; // calibration method

	vector<double> m_globalhisto;	
	vector<double> m_global_max; // Critical points of the whole map kernel
	vector<double> m_global_min; //

	TString	fp; // folder path (calib from image)
	TString fn;	// file name (calib from image)
	TString m_outputName; // output name

	bool globalEstimationSuccess;
	pair<double, double> m_glob_const;
	pair<double, double> m_glob_sig;
	pair<double, double> m_glob_a;
	pair<double, double> m_glob_b;
	pair<double, double> m_glob_c;
	pair<double, double> m_glob_t;
	double m_e0_bound = 0.;
	pair<double, double> m_glob_e0;

    double m_linearPeak1;
    double m_linearPeak2;
    double m_lowenPeak;  
    Double_t m_lowen_fitParams[4]; // constant, sigma, c, t
};

#endif
