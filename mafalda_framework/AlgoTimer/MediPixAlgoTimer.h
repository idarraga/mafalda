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

#ifndef MediPixAlgoTimer_h
#define MediPixAlgoTimer_h

#include <map>
#include <iostream>
#include <TROOT.h>
#include <TStopwatch.h>
#include <TString.h>

class MediPixAlgoTimer {

 public:
  MediPixAlgoTimer();
  virtual ~MediPixAlgoTimer() { };

  void SetupATimer(TString);
  void ContinueATimer(TString);
  void StopATimer(TString);
  Double_t GetRealTime(TString);
  Double_t GetCpuTime(TString);

  Double_t GetLastSlotRealTime(TString);
  Double_t GetLastSlotCPUTime(TString);

  void DumpTimers();

 private:

  Int_t nTimers;
  std::map<TString, TStopwatch> timers; // elapsed time
  std::map<TString, TStopwatch> timersInstant; // time between last and previous stop

};


#endif
