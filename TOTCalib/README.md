How to use macros:

-------------------- runTOTCalib.C --------------------
Macro to run the TOT calibration.
Run it with ROOT6 as: root -l runTOTCalib.C
Input sources file must be MAFOutput root files from TOTCalibrationPreparation algo (in mafalda_framework).

-------------- runSavePixelResolution.C ---------------
Macro to record spectra and fits for each pixel from one source file only.
Run it with ROOT6 as: root -l runSavePixelResolution.C
Source file must be a MAFOutput root file from TOTCalibrationPreparation algo (in mafalda_framework).
It generates a root file (output_pixelResolution.root) that contains histos, kernel and fit functions.
This root file is intended to be used with runExplorePixelResolution.C
If calibration files are given in argument when calling the function, it records calibrated histograms and fits as well (compatible with runExplorePixelResolution.C). The 4 text files (e.g. obtained by running runTOTCalib.C) must be given in the order a,b,c,t. 

------------- runExplorePixelResolution.C --------------
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
