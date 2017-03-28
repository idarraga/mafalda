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

/* Steering the viewer */
#define SEEK_FORWARD     10
#define SEEK_BACKWARDS    5

/* these have to be negatives */
#define JUMP           -100
#define DONT_JUMP        -5
#define DONT_MOVE       -10

/* Look & Feel definitions */

#define MAX_COLORS 256

#define VIS_MOD_GRAY 1
#define VIS_MOD_JET 2
#define VIS_MOD_HOT 3
#define VIS_MOD_COOL 4

#define GRAY_SATURATION_VAL    0
#define GRAY_MAXLIGHTNESS_VAL  0
#define GRAY_MINLIGHTNESS_VAL  1
#define GRAY_MAXHUE_VAL        0
#define GRAY_MINHUE_VAL        0

#define JET_SATURATION_VAL     1
#define JET_MAXLIGHTNESS_VAL   0.5
#define JET_MINLIGHTNESS_VAL   0.6
#define JET_MAXHUE_VAL         230
#define JET_MINHUE_VAL         0

#define HOT_SATURATION_VAL     1
#define HOT_MAXLIGHTNESS_VAL   0.1
#define HOT_MINLIGHTNESS_VAL   1
#define HOT_MAXHUE_VAL         0
#define HOT_MINHUE_VAL         60

#define COOL_SATURATION_VAL    2
#define COOL_MAXLIGHTNESS_VAL  0.8
#define COOL_MINLIGHTNESS_VAL  0.8
#define COOL_MAXHUE_VAL        175
#define COOL_MINHUE_VAL        300

#define MAX_SCALE_VALUE
#define MIN_SCALE_VALUE        0

#define N_SIGNIFICANT_FIGURES  3
//graph->GetHistogram()->GetMinimum(); //That might get 


/*For Jet
  const float  Saturation = 1.0;
  const float  Maxlightness = 0.4;
  const float  Minlightness = 0.4;
  const float  MaxHue = 240; // from 0 to 360
  const float  MinHue = 0;
  const int    MaxColors = 256;   
  int          palette[MaxColors];
  int          index;
  float        saturation, lightness, hue, r, g, b, rv, gv, bv;
  TColor       *color;
  unsigned int failures = 0;
  
  for (int i=0 ; i<MaxColors ; i++) {
  index = palette[i] = MaxColors+1+i;     
  color = new TColor(index, 0, 0, 0);
  hue = MaxHue-(i+1)*((MaxHue-MinHue)/MaxColors);
  lightness = Maxlightness-(i+1)*((Maxlightness-Minlightness)/MaxColors);
  saturation = Saturation;
  //std::cout << hue << ", ";
  color->HLStoRGB(hue, lightness, saturation, r, g, b);
  color->SetRGB(r, g, b);
  
*/

/*
  For Gray
  const float  Saturation = 0;
  const float  Maxlightness = 0.0;
  const float  Minlightness = 1.0;
  const float  MaxHue = 0; // from 0 to 360
  const float  MinHue = 0;
  const int    MaxColors = 256;   
  int          palette[MaxColors];
  int          index;
  float        saturation, lightness, hue, r, g, b, rv, gv, bv;
  TColor       *color;
  unsigned int failures = 0;
 
  for (int i=0 ; i<MaxColors ; i++) {
    index = palette[i] = MaxColors+1+i;     
    color = new TColor(index, 0, 0, 0);
    hue = MaxHue-(i+1)*((MaxHue-MinHue)/MaxColors);
    lightness = Maxlightness-(i+1)*((Maxlightness-Minlightness)/MaxColors);
    saturation = Saturation;
    //std::cout << hue << ", ";
    color->HLStoRGB(hue, lightness, saturation, r, g, b);
    color->SetRGB(r, g, b);

*/

/*
For Hot
  const float  Saturation = 1;
  const float  Maxlightness = 0.0;
  const float  Minlightness = 1.0;
  const float  MaxHue = 0; // from 0 to 360
  const float  MinHue = 60;
  const int    MaxColors = 256;   
  int          palette[MaxColors];
  int          index;
  float        saturation, lightness, hue, r, g, b, rv, gv, bv;
  TColor       *color;
  unsigned int failures = 0;
 
  for (int i=0 ; i<MaxColors ; i++) {
    index = palette[i] = MaxColors+1+i;     
    color = new TColor(index, 0, 0, 0);
    hue = MaxHue-(i+1)*((MaxHue-MinHue)/MaxColors);
    lightness = Maxlightness-(i+1)*((Maxlightness-Minlightness)/MaxColors);
    saturation = Saturation;
    //std::cout << hue << ", ";
    color->HLStoRGB(hue, lightness, saturation, r, g, b);
    color->SetRGB(r, g, b);

*/

/*
For Cool
  const float  Saturation = 1;
  const float  Maxlightness = 0.0;
  const float  Minlightness = 1.0;
  const float  MaxHue = 0; // from 0 to 360
  const float  MinHue = 60;
  const int    MaxColors = 256;   
  int          palette[MaxColors];
  int          index;
  float        saturation, lightness, hue, r, g, b, rv, gv, bv;
  TColor       *color;
  unsigned int failures = 0;
 
  for (int i=0 ; i<MaxColors ; i++) {
    index = palette[i] = MaxColors+1+i;     
    color = new TColor(index, 0, 0, 0);
    hue = MaxHue-(i+1)*((MaxHue-MinHue)/MaxColors);
    lightness = Maxlightness-(i+1)*((Maxlightness-Minlightness)/MaxColors);
    saturation = Saturation;
    //std::cout << hue << ", ";
    color->HLStoRGB(hue, lightness, saturation, r, g, b);
    color->SetRGB(r, g, b);

*/
