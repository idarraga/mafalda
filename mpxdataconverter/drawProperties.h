#include <iostream>
#include "TROOT.h"
#include "TColor.h"
#include "draw_def.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TStyle.h"

class MPXPalette {

 private:

  Int_t           MaxColors;   
  Int_t           * palette;
  Int_t           index;
  Float_t         lightness, hue, r, g, b, rv, gv, bv;
  TColor          * color;
  UInt_t          failures;

 public:

  MPXPalette();
  ~MPXPalette();
  Int_t * generatePalette();
  
  Float_t Saturation;
  Float_t Maxlightness;
  Float_t Minlightness;
  Float_t MaxHue;
  Float_t MinHue;
  Float_t Lightness;


};

class HistoProperties {

 private:
  MPXPalette * histoPalette;

 public:
  HistoProperties();
  ~HistoProperties();
  Bool_t setMode(Int_t, TCanvas *);

};
