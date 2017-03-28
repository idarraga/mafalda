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

#ifndef MediPixAlgoTimer_cxx
#define MediPixAlgoTimer_cxx

#include "MediPixAlgoTimer.h"

MediPixAlgoTimer::MediPixAlgoTimer(){

  nTimers = 0;

}

void MediPixAlgoTimer::SetupATimer(TString timerName){

  timers[timerName] = TStopwatch();
  timers[timerName].Reset();

  timersInstant[timerName] = TStopwatch();
  timersInstant[timerName].Reset();

}

void MediPixAlgoTimer::ContinueATimer(TString timerName){

  // start
  timers[timerName].Start(kFALSE); /* start but don't reset */
  timersInstant[timerName].Start(kTRUE); /* start resetting */

}
void MediPixAlgoTimer::StopATimer(TString timerName){
  timers[timerName].Stop();
  timersInstant[timerName].Stop();
}
Double_t MediPixAlgoTimer::GetRealTime(TString timerName){
  return timers[timerName].RealTime();
}
Double_t MediPixAlgoTimer::GetCpuTime(TString timerName){
  return timers[timerName].CpuTime();
}
Double_t MediPixAlgoTimer::GetLastSlotRealTime(TString timerName){
  return timersInstant[timerName].RealTime();
}
Double_t MediPixAlgoTimer::GetLastSlotCPUTime(TString timerName){
  return timersInstant[timerName].CpuTime();
}

void MediPixAlgoTimer::DumpTimers(){

  std::map<TString, TStopwatch>::iterator timerItr = timers.begin();


  std::cout << "Timers ------------------------------------" << std::endl;
  std::cout << "  AlgoName:     realtime, cputime " << std::endl;
  for( ; timerItr != timers.end() ; timerItr++)
    {
      std::cout << "  " << (*timerItr).first << ": ";
      std::cout << GetRealTime((*timerItr).first) << ", "; 
      std::cout << GetCpuTime((*timerItr).first) << " [s]" << std::endl;
    }
  std::cout << "-------------------------------------------" << std::endl;

}
#endif
