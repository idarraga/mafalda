#include "drawProperties.h"

HistoProperties::HistoProperties(){

  histoPalette = new MPXPalette();

}

HistoProperties::~HistoProperties(){
  delete histoPalette;
}


///////////////////////////////////////////////


Bool_t HistoProperties::setMode(Int_t visMode, TCanvas * c1){


  switch (visMode) {
    
  case VIS_MOD_GRAY:
    histoPalette->Saturation = GRAY_SATURATION_VAL;
    histoPalette->Maxlightness = GRAY_MAXLIGHTNESS_VAL;
    histoPalette->Minlightness = GRAY_MINLIGHTNESS_VAL;
    histoPalette->MaxHue = GRAY_MAXHUE_VAL;
    histoPalette->MinHue = GRAY_MINHUE_VAL;
    break;
    
  case VIS_MOD_JET:
    std::cout << "  [INFO] change to mode Jet" << std::endl;
    histoPalette->Saturation = JET_SATURATION_VAL;
    histoPalette->Maxlightness = JET_MAXLIGHTNESS_VAL;
    histoPalette->Minlightness = JET_MINLIGHTNESS_VAL;
    histoPalette->MaxHue = JET_MAXHUE_VAL;
    histoPalette->MinHue = JET_MINHUE_VAL;
    break;
    
  case VIS_MOD_HOT:
    histoPalette->Saturation = HOT_SATURATION_VAL;
    histoPalette->Maxlightness = HOT_MAXLIGHTNESS_VAL;
    histoPalette->Minlightness = HOT_MINLIGHTNESS_VAL;
    histoPalette->MaxHue = HOT_MAXHUE_VAL;
    histoPalette->MinHue = HOT_MINHUE_VAL;
    break;
    
  case VIS_MOD_COOL:
    histoPalette->Saturation = COOL_SATURATION_VAL;
    histoPalette->Maxlightness = COOL_MAXLIGHTNESS_VAL;
    histoPalette->Minlightness = COOL_MINLIGHTNESS_VAL;
    histoPalette->MaxHue = COOL_MAXHUE_VAL;
    histoPalette->MinHue = COOL_MINHUE_VAL;
    break;
    
  default:

    break;
  }

  Int_t * palette = histoPalette->generatePalette();
  gStyle->SetPalette(256, palette);
  gStyle->SetPalette(1, 0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameFillColor(palette[0]);
  //c1->SetFrameFillColor(1,0);
  c1->SetBorderMode(0);
  c1->ToggleEventStatus(); // see Event Status bar

  return true;
}


MPXPalette::MPXPalette(){

  Saturation = -1;
  Maxlightness = -1;
  Minlightness = -1;
  MaxHue = -1;
  MinHue = -1;
  Lightness = -1;

  MaxColors = 256;
  palette = new Int_t[MaxColors];
  index = -1;
  lightness = -1.0; hue = -1.0; r = -1.0; 
  g = -1.0; b = -1.0; rv = -1.0; gv = -1.0; 
  bv = -1.0;
  failures = 0;

}

MPXPalette::~MPXPalette(){
  delete palette;
  delete color;
}


Int_t * MPXPalette::generatePalette(){

  for (int i = 0 ; i < MaxColors ; i++)
    {
      index = palette[i] = MaxColors+1+i;     
      color = new TColor(index, 0, 0, 0);
      lightness = Maxlightness-(i+1)*((Maxlightness-Minlightness)/MaxColors);
      hue = MaxHue-(i+1)*((MaxHue-MinHue)/MaxColors);
      color->HLStoRGB(hue, lightness, Saturation, r, g, b);
      color->SetRGB(r, g, b);

      // check if the color can be defined
      //gVirtualX->GetRGB(index, rv, gv, bv);
      //gGXW->GetRGB(index, rv, gv, bv);
      
      if (r != rv || g != gv || b != bv) {
	failures++;
	//palette[i] =  i ? palette[i-1] : 1;
      }
    } 

  if (failures)
    std::cout << "  [WARNING] palette(): couldn't allocate " 
	      << failures << " of " << MaxColors << " colors" << std::endl;

  return palette;
}
