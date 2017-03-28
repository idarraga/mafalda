/*
 * 	Copyright 2012 John Idarraga
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

#ifndef __SimpleClustering_h
#define __SimpleClustering_h

#include "MPXAlgo/MediPixAlgo.h"

// inner single word
#define __inner_single_word  0x10
// low border --> Single 0x02
#define __low_border_word    0x02
// top border --> Single 0x10
#define __top_border_word    0x10
// left border --> Single 0x04
#define __left_border_word   0x04
// right border --> Single = 0x08
#define __right_border_word  0x08

// bottom left corner
#define __blc_word  0x01
// bottom right corner
#define __brc_word  0x02
// up left corner
#define __ulc_word  0x04
// up right corner
#define __urc_word  0x08

class SimpleClusterContainer : public CandidateContainer {

 public:

  SimpleClusterContainer(MediPixAlgo *);
  ~SimpleClusterContainer(){};

  void InsertSingleHit(int x, int y){ m_singleHits.push_back( make_pair(x, y) ); };
  void CleanUp(){ m_singleHits.clear(); };

  vector<pair<int, int> > GetSingleHits(){ return m_singleHits; };
  int GetNSingleHits(){return (int)m_singleHits.size();};

 public:

  // Should get to be private but I need to change things in the implementation
  vector< pair<int, int> > m_singleHits;

  ClassDef(SimpleClusterContainer, 1)
};


class SimpleClustering : public MediPixAlgo {

public:

	SimpleClustering();
	virtual ~SimpleClustering() { };

	// You ought to implement Init(), Execute() and Finalize()
	//  when you inherit from MediPixAlgo.  This model gives you
	//  direct access to data and services.
	void Init();
	void Execute();
	void Finalize();

private:

	// all blobs
	SimpleClusterContainer * m_scc;

	// config
	double m_maxOcc;

	ClassDef(SimpleClustering, 1)
};

#endif
