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

/**
 *
 * This algorithm uses the output of BlobsFinder and classifies each
 *   blob under one of the following basic species:
 * 
 * - single hits
 * - double hits
 *        *        
 *          * , * *
 * - triple hits
 *        *
 *        * *
 * - quads
 *        * *
 *        * *
 * - long gamma
 *   gamma long-shapes due to scattered compton electrons, when
 *   gammas get into the Si at high angle.  Smallest angle between
 *   the momentum vector of the photon and 
 * - alphas --> uses parameters that may be adjusted !
 * - curly shapes
 *
 */

#ifndef PRBasicSpecies_h
#define PRBasicSpecies_h

#include "MPXAlgo/MediPixAlgo.h"
#include "MPXAlgo/Highlighter.h"

// This algo is using an object put in the StoreGate
//  by BlobsFinder, I need to include the header file
//  in order to know the object
#include "BlobsFinder/BlobsFinder.h"
#include "CalibrationLoader/CalibrationLoader.h"

#define MAX_PIECES_PER_FRAME 500

typedef struct {

  Int_t nSingleHits;
  Int_t nDoubleHits;
  Int_t nTripleHits;
  Int_t nQuadHits;
  Int_t nLongGamma;
  Int_t nHeavyBlobs;
  Int_t nHeavyTracks;
  Int_t nMip;
  Int_t nCurly;
  Int_t nNotBasicType;
  Int_t nAll;
} BasicSpeciesCntr_struct ;

typedef struct {
  Int_t nHits;
  Int_t nCounts;
} BasicSpeciesGeneralInfo_struct ;

//Float_t & mindistanceToLine, Float_t & fractionOfPixelsAtMindistance)

class PRBasicSpecies : public MediPixAlgo , public CalibrationLoader {

public:

  PRBasicSpecies();
  virtual ~PRBasicSpecies() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

  void SeparateBasicSpecies();
  bool SingleHit(vector<blob>::iterator);
  bool DoubleHit(vector<blob>::iterator);
  bool TripleHit(vector<blob>::iterator);
  bool QuadHit(vector<blob>::iterator);
  bool LongGamma(vector<blob>::iterator);
  bool Mip(vector<blob>::iterator);
  bool HeavyBlob(vector<blob>::iterator);
  bool HeavyTrack(vector<blob>::iterator);
  bool Curly(vector<blob>::iterator);

  void FillValuesForDisplay(Highlighter *, Int_t, blob);

  void ClearNtupleVars();

  void SetMIPOnOff(Bool_t flag){m_mipON = flag;};
  void SetHeavyBlobOnOff(Bool_t flag){m_HeavyBlobON = flag;};
  void SetHeavyTrackOnOff(Bool_t flag){m_HeavyTrackON = flag;};
  void SetCheckOverFlowOnOff(Bool_t flag){m_checkOverFlow = flag;};

  Bool_t FrameContainsOverFlowPixels();

  void DrawLines(blob);

private:

  AllBlobsContainer * m_aB;

  // to ntuple
  Int_t m_frameId;
  BasicSpeciesCntr_struct m_basicSpeciesCntr;
  BasicSpeciesGeneralInfo_struct m_basicGeneral;

  //HeavyBlob_struct m_heavyBlobs;
  vector<Int_t> m_clusterSize;
  vector<Float_t> m_clusterSizeWidth;
  vector<Float_t> m_clusterSizeEllipse;
  vector<Float_t> m_circleArea;
  vector<Float_t> m_counts;

  // selection parameters
  Float_t m_circleToEllipseMin;
  Float_t m_circleToEllipseMax;
  Int_t m_longGammaMax;
  Int_t m_nInnerPixelsCut;

  // Mips
  Float_t m_mipChisquareDof;
  Float_t m_mindistanceToLine;
  Float_t m_fractionOfPixelsAtMindistance;
  Int_t m_mipNPixels;
  // curly
  Int_t m_maxCurlyPixels;

  // switches
  Bool_t m_HeavyBlobON;
  Bool_t m_HeavyTrackON;
  Bool_t m_mipON;
  Bool_t m_checkOverFlow;

  // other vars
  Int_t m_overflowPixelsCntr;

  // highlighters
  Highlighter * m_highCurly;
  Highlighter * m_highHeavyBlobs;
  Highlighter * m_highHeavyTracks;
  Highlighter * m_highMips;
  Highlighter * m_highUnknown;

  ClassDef(PRBasicSpecies, 1)
};

#endif
