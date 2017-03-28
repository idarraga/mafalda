 /*
 * 	Copyright 2014 John Idarraga
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

#ifndef __DosSummary_h
#define __DosSummary_h

#include "MPXAlgo/MediPixAlgo.h"

#include "DosSummaryReport.h"
#include "BlobsFinder/BlobsFinder.h"

class DosSummary : public MediPixAlgo {

public:

  DosSummary();
  virtual ~DosSummary() { };

  // You ought to implement Init(), Execute() and Finalize()
  //  when you inherit from MediPixAlgo.  This model gives you
  //  direct access to data and services.
  void Init();
  void Execute();
  void Finalize();

private:

  AllBlobsContainer * m_aB;

  DosSummaryReport m_dos_summary_report;

  ClassDef(DosSummary, 1)
};

#endif
