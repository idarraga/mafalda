 /*
 * 	Copyright 2010 John Idarraga, Mathieu Benoit
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

#ifndef BeamSpread_h
#define BeamSpread_h

#include "MPXAlgo/MediPixAlgo.h"

class AllBlobsContainer;

class BeamSpread : public MediPixAlgo {

public:

  BeamSpread();
  virtual ~BeamSpread() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

private:

  AllBlobsContainer * m_allBlobs;
  Double_t m_spreadDist;
  Double_t m_pixelSizeX;
  Double_t m_pixelSizeY;

  ClassDef(BeamSpread, 1)
};

#endif
