/*
 * 	Copyright 2011 John Idarraga
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

#ifndef AnalysisCore_defs_h
#define AnalysisCore_defs_h

namespace MPXDefs {

  enum Flags {
    SERIALIZE_ME = 0,
    DO_NOT_SERIALIZE_ME,
  };

  enum SpecialObjs {
    NIL = 0,
    CONF_INT,
    CONF_FLOAT,
    CONF_DOUBLE,
    CONF_BOOL,
    CONF,  // Configuration objects.  They won't be erased over the entire run
    // Starting from here objects to erase on an event per event basis.
    VIS,       // Visualization drawable objects.
    VIS_SKIP,  // Visualization skip frame objects.
    REGULAR, // all regular objects erasable on an event per event basis
  };

}

#define LEFT_MSG_TAB 30

namespace MSG  {
  enum Level {
    NIL = 0,
    LOOP_DEBUG,
    DEBUG,
    INFO,
    WARNING,
    ERROR,
    FATAL,
    ALWAYS,
    NUM_LEVELS
  };

  typedef struct {
    Int_t endLevel;
  } Endreq ;

}

#define MAX_OCCUPANCY_CONF_LINE 1024

#endif
