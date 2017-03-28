# Modified by John Idarraga, 2008
# Medipix Group, Universite de Montreal	
#
# Makefile for the ROOT test programs.  This Makefile shows nicely how
# to compile and link applications using the ROOT libraries on all
# supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000
# -*- mode: makefile -*-
#
# Makefile containing platform dependencies for ROOT based projects.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ROOTCONFIG = $(shell which root-config)
ROOTCINT   = $(shell which rootcint)
ROOTCLING  = $(shell which rootcling)


shell         = /bin/bash

ARCH         := $(shell $(ROOTCONFIG) --prefix=$(ROOTSYS) --arch)
PLATFORM     := $(shell $(ROOTCONFIG) --prefix=$(ROOTSYS) --platform)

#ARCH         := $(shell $(ROOTSYS)/bin/root-config --prefix=$(ROOTSYS) --arch)
#PLATFORM     := $(shell $(ROOTSYS)/bin/root-config --prefix=$(ROOTSYS) --platform)


CXX           =
ObjSuf        = o
SrcSuf        = cpp
HdrSuf        = h
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
HASTHREAD    := $(shell $(ROOTCONFIG) --has-thread)

#ROOTCFLAGS   := $(shell $(ROOTSYS)/bin/root-config --prefix=$(ROOTSYS) --cflags)
#ROOTLDFLAGS  := $(shell $(ROOTSYS)/bin/root-config --prefix=$(ROOTSYS) --ldflags)
#ROOTLIBS     := $(shell $(ROOTSYS)/bin/root-config --prefix=$(ROOTSYS) --libs)
#ROOTGLIBS    := $(shell $(ROOTSYS)/bin/root-config --prefix=$(ROOTSYS) --glibs)
#HASTHREAD    := $(shell $(ROOTSYS)/bin/root-config --prefix=$(ROOTSYS) --has-thread)

# default optimized
EXTRACXXFLAGS = -O2 -W -Wall

opt:	EXTRACXXFLAGS=-O2 -W -Wall

dbg:	EXTRACXXFLAGS=-g3 -W -Wall

opt:	mafalda

dbg:	mafalda

ifeq ($(ARCH),linux)
# Linux with egcs, gcc 2.9x, gcc 3.x, gcc 4.x (check templates)
CXX           = g++ 
CXXFLAGS      = -W -Wall -fPIC $(EXTRACXXFLAGS)
LD            = g++
LDFLAGS       = -O2
SOFLAGS       = -shared
endif

ifeq ($(ARCH),linuxx8664gcc)
# Linux with egcs, gcc 2.9x, gcc 3.x, gcc 4.x (check templates)
CXX           = g++ 
CXXFLAGS      = -W -Wall -fPIC $(EXTRACXXFLAGS)
LD            = g++
LDFLAGS       = -O2
SOFLAGS       = -shared
endif

ifeq ($(ARCH),macosx)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++ 
CXXFLAGS      = -g3 -W -Wall
LD            = g++
LDFLAGS       = -O0
SOFLAGS       = -dynamiclib
endif

ifeq ($(ARCH),macosx64)
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++ 
CXXFLAGS      = -g3 -W -Wall -m64
LD            = g++
LDFLAGS       = -O0
SOFLAGS       = -dynamiclib
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

CXXFLAGS     += $(ROOTCFLAGS) -I./
LDFLAGS      += $(ROOTLDFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

EXPLLINKLIBS  = $(LIBS)
EXPLLINKLIBS += $(GLIBS)

#------------------------------------------------------------------------------
MAFO := $(shell ./buildList.scripts "ANO" $(ObjSuf) $(SrcSuf))
MAFO += MediPixAnalysisCoreDict.o

DICTDEP := $(shell ./buildList.scripts "DICT" $(HdrSuf) $(HdrSuf))

CINCDEP :=  $(shell ./buildList.scripts "CINCDEP" $(ObjSuf) $(SrcSuf))

MAFSO       = libMediPixAnalysisCore.$(DllSuf)
WZAN         = MediPixAnalysis$(ExeSuf)
WZANLIB      = $(MAFSO)
OBJS          = $(MAFO)
PROGRAMS      = $(WZAN)
#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
	
all:	mafalda

mafalda:            $(PROGRAMS)
	@echo "***                                MAFalda ***"
	@echo "*** Analysis Framework for pixel detectors ***"
	@echo "*** John Idarraga <idarraga@cern.ch>, 2013 ***"


$(MAFSO):	$(MAFO)
		@echo " ... linking ... "
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(EXPLLINKLIBS)

$(WZAN):       $(MAFSO)
#		$(LD) $(LDFLAGS)  $(WZANLIB) $(LIBS) \
#		   $(OutPutOpt)$(WZAN)
		@echo "done"

clean:
		@rm -f $(OBJS) core *Dict.* *.so *.o

dictclean:	
		@rm -f *Dict.*

distclean:      clean
		@rm -f $(PROGRAMS) $(MAFSO) $(WZANLIB) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d .def so_locations
		@rm -rf cxx_repository
		-@cd RootShower && $(MAKE) distclean

.SUFFIXES: .$(SrcSuf)


MediPixAnalysisCoreDict.$(SrcSuf): $(DICTDEP)
	@echo "Generating dictionary $@ ..."
	$(ROOTCLING) -f $@ -c $^
#	$(ROOTCINT) -f $@ -c $^
#	rootcint -f $@ -c $^

MediPixAnalysisCoreDict.$(SrcSuf).o:	
	@echo $(DICTDEP)

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) $(CINCDEP) -c $< -o $@


newalgo:	
	@$(shell) newalgo.sh
