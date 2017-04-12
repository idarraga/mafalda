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

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <map>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;


#define __npars_surrogate			    4
#define __npars_lowe_fitfunc            7
#define __min_tmathprobtest_val      0.90
#define __max_fit_tries                10
#define __fit_pars_randomization_max 1000
#define __fit_pars_randomization_min   20

#define __fraction_of_height_range_id 0.5

// Prototypes
double GausFuncAdd(double * x, double * par);
void printProgBar( int );
#define __fitfunc_lowen_npars  7
double fitfunc_lowen(double * x, double * par);
double surrogatefunc_calib(double * x, double * par);

// Calibration Handler for each source
class CalibHandler {

public:

	CalibHandler(string);
	~CalibHandler(){};
	string GetSourcename () { return m_sourceName; };
	map<int, double> GetCalibPoints(){ return m_calibPoints; };
	map<int, int> GetCalibPointsRegion(){ return m_calibPointsRegion; };
	double GetOneEnergyMatch(int i){return m_calibPoints[i]; };

	enum {
		__linear_reg = 0,
		__lowenergy_reg,
		__other_reg
	};

private:


	// TOT, Energy for the particular source or fluorescence data
	map<int, double> m_calibPoints;
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

	TOTCalib();
	TOTCalib(TString, TString, int minpix, int maxpix, int maxtot, Long64_t nFrames);
	TOTCalib(TString, TString, int minpix, int maxpix, int maxtot, Long64_t nFrames, TOTCalib *);

	void SetupJob(TString, TString, int minpix, int maxpix, int maxtot, Long64_t nFrames);

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

		vector< pair<double, double> > pointsSave;
		vector< double > pointsSaveSigmas;
		vector< double > pointsSaveConstants;
		vector< double > pointsSave_ia;
		vector< double > pointsSave_ib;
		vector< double > pointsSave_ic;
		vector< double > pointsSave_it;
		vector< int > calibTOTPeaks;
		vector<pair<double,double> > linearpairs;
		vector<int> peakFitStatus;
	};

	TH1I * GetHisto(int, TString extra = "");
	TF1 * GetKernelDensityFunction(int);
	int GetCriticalPoints(int i , vector<double> &, vector<double> &);
	int GetSign(double slope);
	map<int, vector<double> > GetMaxPeaksIdentified(){ return m_critPointsMax; };
	CalibHandler * GetCalibHandler(){ return m_calhandler; };
	vector< pair<double, double> > GetCalibPoints(int pix){ return m_calibPoints[pix]; };
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

	void Blender(TOTCalib * , TOTCalib *, TOTCalib *, TString, int = 0);
	void Blender(TOTCalib * , TOTCalib *, TString, int = 0);
	void Blender(TOTCalib * s2, TString outputName, int = 0);
	void Blender(TString, int = 0);
    
    // Methods for calibration
    enum {
        __standard = 0,
        __jakubek // Method proposed in J. Jakubek / Nuclear Instruments and Methods in Physics Research A 633 (2011) S262–S266
    };

	void ProcessOneSource(TOTCalib * s, store * sto, TGraphErrors * g, int pix, int & cntr);
	void ReorderSources();

	int PeakFit(TOTCalib *, int, int, TF1 *, TH1 *, store *);
    int PeakFit(TOTCalib *, int, int, TF1 *, TH1 *, store *, int);    
	TF1 * FittingFunctionSelector(double, TOTCalib *, int);
	void GetLinearFit(double & a, double & b, vector< pair<double,double> > p);
	void RandomFitParameters(TF1 * f, TH1 * h, int tot, TOTCalib*);

	vector<pair<double, double> > Extract_E_TOT_Points(int, TOTCalib * );
	int GetNumberOf_E_TOT_Points (TOTCalib * s);

	double DerivativeFivePointsStencil(TF1 *, double, double);
	double VectorSum(vector<double>);

	int XYtoX(pair<int, int> pix, int dimX);
	pair<int, int> XtoXY(int X, int dimX);
	TH2I * EntriesPlots(int);
	//TF1 * CreateKernelDensityFunction(TH1I *, double);
	TF1 * CreateKernelDensityFunction(int, vector<double>, double);

	void SetGlobalThresholdEnergyAndErr(double th, double th_err) { m_thresholdEnergy = th; m_thresholdEnergy_Err = th_err; };

	void FillHisto(int, int);

	void Finalize();
	TH1F * CreateParameterHistogram(vector<double>, TString);

	void DrawFullPixelCalib(int);
	void DrawFullPixelCalib(int, int);

	void SetVerboseLevel(int v) {m_verbose = v;};
	enum {
		__VER_DEBUG_LOOP = 0,
		__VER_DEBUG,
		__VER_INFO
	};

	double GetKernelBandWidth(){ return m_bandwidth; };
	double SetKernelBandWidth(double b) {return m_bandwidth = b; };

	int GetNBins() { return m_nbins; };

private:
	//////////////////////////////////////////////////////////////////
	// histograms and info for all pixels;
	//vector<TH1I *> m_calibhistos;
	vector<vector<double> > m_calibhistos;
	// key = pixel, values = points selected for calibration
	map<int, vector< pair<double, double> > > m_calibPoints;
	map<int, vector<int> > m_calibTOTPeaks;
	map<int, vector<double> > m_calibPointsSigmas;
	map<int, vector<double> > m_calibPointsConstants;
	map<int, vector<double> > m_calibPoints_ia;
	map<int, vector<double> > m_calibPoints_ib;
	map<int, vector<double> > m_calibPoints_ic;
	map<int, vector<double> > m_calibPoints_it;

	int m_nbins;
	int m_histoRebinning;
	Long64_t m_nFrames;
	int m_minpix, m_maxpix;
	TF1 ** m_kerneldensityfunctions;

	// key = pixel, val = vector of critical points
	map<int, vector<double> > m_critPointsMax;
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

	TH1F * m_par_a;
	TH1F * m_par_b;
	TH1F * m_par_c;
	TH1F * m_par_t;
	TH1F * m_surr_prob;

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
	vector<TF1 * > m_extra_tf1_to_erase;

	// Randon generators
	unsigned m_ranseed_time;
	TRandom1 * m_rand1;

	// Matrix size
	int __matrix_size;
	int __matrix_height;
	int __matrix_width;

};

#endif
