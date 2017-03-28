/*
 * 	Copyright 2008 John Idarraga
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

#ifndef ConfigurationValue_h
#define ConfigurationValue_h

#include <string>

#include <TROOT.h>

#include "AnalysisCore/AnalysisCore_defs.h"
#include "MPXStoreGate/CandidateContainer.h"

using namespace std;

class MediPixAlgo;

/**
 * This class handles values to visualize and modify from the MediPix
 *  viewer.  Note the hyphen in the definitions in LinkDef.h
 *
 */

template<class T> 
class ConfigurationValue : public CandidateContainer {

 private:
  string m_valuesConfName;
  T * m_valuesConfVal; // pointer to actual values in the Algorithm
  
 public:
  ConfigurationValue();
  ConfigurationValue(MediPixAlgo *, T *, 
		     const char *, 
		     MPXDefs::SpecialObjs);
  ~ConfigurationValue();  // non need for virtual constructor, child class.
  
  T GetConfigValue(){ return * m_valuesConfVal; };
  string GetConfigName(){ return m_valuesConfName; };

  void SetConfigValue(T val);
  
  ClassDef(ConfigurationValue, 1)
};

#endif
