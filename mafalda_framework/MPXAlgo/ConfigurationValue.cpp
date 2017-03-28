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

#ifndef ConfigurationValue_cpp
#define ConfigurationValue_cpp

#include "ConfigurationValue.h"

templateClassImp(ConfigurationValue)

/**
 * This class handles values to visualize/tweak from the MediPix
 *  viewer.
 *
 */

// default constructor
template<class T> 
ConfigurationValue<T>::ConfigurationValue()
  : CandidateContainer() {
}

template<class T> 
ConfigurationValue<T>::ConfigurationValue(MediPixAlgo * algo, T * val, 
					  const char * valname, 
					  MPXDefs::SpecialObjs objType) 
  : CandidateContainer(algo, objType) {
  
  m_valuesConfName = valname;

  //m_valuesConfVal = new T; // <<----- check this out !!
  m_valuesConfVal = val;

}

template<class T> 
ConfigurationValue<T>::~ConfigurationValue()
{
  delete m_valuesConfVal;
}

template<class T> 
void ConfigurationValue<T>::SetConfigValue(T val){ 
  *m_valuesConfVal = val; 
}

// In a shared lib scenario (or if you define the template in a different file, like here) 
//  you have to tell what are you going to instantiate
template class ConfigurationValue<Int_t>;
template class ConfigurationValue<Float_t>;
template class ConfigurationValue<Double_t>;
template class ConfigurationValue<Bool_t>;

#endif
