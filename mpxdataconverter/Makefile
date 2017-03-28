#GLOBAL_CFLAGS= -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
#GLOBAL_CFLAGS= `getconf LFS_CFLAGS`
GLOBAL_LDFLAGS= `getconf LFS_LDFLAGS`
INCLUDE=
IBS=
OBJS = allpix_dm.o listHandler.o MediPixWriteToEntuple.o SelDict.o
SRCS = allpix_dm.cc listHandler.cpp MediPixWriteToEntuple.cc
HDRS = allpix_dm.h listHandler.h MediPixWriteToEntuple.h allpix_dm_consts.h

ifndef ROOTSYS
  DO_FLAG=FALSE  
else
  DO_FLAG=TRUE
  INCLUDE += $(shell $(ROOTSYS)/bin/root-config --cflags)
  LIBS += $(shell $(ROOTSYS)/bin/root-config --libs)
endif

GLOBAL_CFLAGS=-g3 -W -Wall
#GLOBAL_CFLAGS=-O2 -W -Wall
PROGNAME=mpxdataconverter

all:	$(PROGNAME) raw_to_pixeldm
	@echo "[DONE]"
	
ifeq ($(DO_FLAG),FALSE)
	@echo "ROOTSYS variable should by set to run make command"
else
	@rm -f check_file
	@rm -f touch_file
	@make $(PROGNAME)
endif
	@echo ""

########################################

$(PROGNAME):	msg1 MediPixROOTConverter.o ${OBJS}
	g++ -o $@ MediPixROOTConverter.o ${OBJS} ${LIBS}

msg1:	
	@echo "----------------------------------------------------"
	@echo "             Building converter: Pixelman to MAFalda"
	@echo "                                 John Idarraga, 2008"
	@echo "                                    idarraga@cern.ch"
	@echo "----------------------------------------------------"

########################################

raw_to_pixeldm:	msg2 allpix_dm.o MediPixWriteToEntuple.o SelDict.o raw_to_pixeldm.cpp
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c raw_to_pixeldm.cpp
	g++ -o $@ raw_to_pixeldm.o allpix_dm.o MediPixWriteToEntuple.o SelDict.o $(GLOBAL_CFLAGS) $(INCLUDE) ${LIBS}

msg2:	
	@echo "----------------------------------------------------"
	@echo "            Building converter: USBPixRaw to MAFalda"
	@echo "                                 John Idarraga, 2011"
	@echo "                                    idarraga@cern.ch"
	@echo "----------------------------------------------------"	
	
#########################################

MediPixROOTConverter.o:	MediPixROOTConverter.cpp
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c MediPixROOTConverter.cpp

allpix_dm.o:	allpix_dm.cc
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c allpix_dm.cc

listHandler.o:	listHandler.cpp
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c listHandler.cpp

control.o:	control.cpp
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c control.cpp

MediPixWriteToEntuple.o:	MediPixWriteToEntuple.cc
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c MediPixWriteToEntuple.cc

drawProperties.o:	drawProperties.cpp
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c drawProperties.cpp

SelDict.o:	
	@echo "Generating dictionary $@..."
	@rootcint -f SelDict.cpp -c MediPixWriteToEntuple.h allpix_dm.h LinkDef.h
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c SelDict.cpp

clean:
	rm -f *.o MediPixViewer MediPixROOTConverter
