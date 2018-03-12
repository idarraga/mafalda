# TOTCalib

## How to install and run:

1) you must have cmake installed first.
2) in the TOTCalib folder, do:

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

3) write your macro and run it from the macro folder as:
```
$ root -l macro.C
```

## How to use macros:

### runTOTCalib.C 

Macro to run the TOT calibration (first implementation).
Run it with ROOT6 as: root -l runTOTCalib.C
Input source files must be MAFOutput root files from TOTCalibrationPreparation algo (in mafalda_framework).

### runTOTCalib2.C 

Macro to run the TOT calibration (2nd implementation, much faster than the first one).
Run it with ROOT6 as: root -l runTOTCalib2.C
Input source files must be MAFOutput root files from TOTCalibrationPreparation algo (in mafalda_framework).
It is an implementation of the method published in NIMPR A 633 (2011) S262-S266 (2nd method mentionned).
It requires 3 X-ray sources: 2 in the linear region and one close to the threshold (i.e. non-gaussian spectrum)
NOTES:
- There must be exactly 2 input root files for the peaks in the linear region and one for the low energy region
- Only one peak per input root file can be used (must be selected in the macro)
- Expected peak energies must be defined in the CalibHandler constructor (in TOTCalib.cpp), and must be sorted in ascending values.
- Needs the ROOT Spectrum library 
- Fitting range for low energy source is defined by the histogram range set in the macro
- Fitting range for linear peaks is calculated from peak amplitude (rebin if histo is noisy, i.e it has dips in the gaussian peak !). The fraction of height amplitude is set at 0.5 but can be modified in TOTCalib.h with the value (__fraction_of_height_range_id). See PeakFit2_gaussian() for details.
- Peak search optimization for the linear region sources: in case the user asks for optimization with the function OptimizeOnePeak(thl) in the macro, it removes too small found peaks (i.e with amplitude smaller than 'thl' x 'amplitude of the highest peak'). If only the highest remains but several peaks were defined, it is used as the peak to fit. This can help in the situations when the expected peaks are wrongly identified.
TO DO:
- adapt for more than 2 peaks in the linear regions
- speed up the calibration by processing every steps pixel by pixel (the use of the C++ type unordered_map to save the peaks of all the 65k pixels slows done the algorithm significantly) 

### runSavePixelResolution.C 

Macro to record spectra and fits for each pixel from one source file only.
Run it with ROOT6 as: root -l runSavePixelResolution.C
Source file must be a MAFOutput root file from TOTCalibrationPreparation algo (in mafalda_framework).
It generates a root file (output_pixelResolution.root) that contains histos, kernel and fit functions.
This root file is intended to be used with runExplorePixelResolution.C
If calibration files are given in argument when calling the function, it records calibrated histograms and fits as well (compatible with runExplorePixelResolution.C). The 4 text files (e.g. obtained by running runTOTCalib.C) must be given in the order a,b,c,t. 

### runExplorePixelResolution.C

Macro to display maps of values recorded by runSavePixelResolution.C.
Run it with ROOT6/CLING as:
root [0] TFile f("output_pixelResolution.root");
root [1] SingleHitFitMeans->Draw("colz"); // other histograms available
root [2] gStyle->SetOptStat(0); // remove stat box
root [3] c1->AddExec("",".x runExplorePixelTOTResolution.C");  
root [4] c1->ToggleEventStatus(); // to display pixel coordinate in the status bar (bottom line of canvas)
root [5] gStyle->SetPalette(52); // Better vizualisation with black & white.
When pointing a specific pixel with the mouse it opens a canvas with its histogram, kernel and fit function.
If calib files were used in runSavePixelResolution.C, a second canvas is opened which displays the calibrated spectrum with a gaussian fit.

