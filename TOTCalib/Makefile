ifndef ROOTSYS
  DO_FLAG=FALSE  
else
  DO_FLAG=TRUE
  ROOTCONFIG = $(shell which root-config)
  ROOTCINT   = $(shell which rootcint)
  INCLUDE += $(shell $(ROOTSYS)/bin/root-config --prefix=$(ROOTSYS) --cflags)
  LIBS += $(shell $(ROOTSYS)/bin/root-config --prefix=$(ROOTSYS) --libs)
endif

# default optimized
GLOBAL_CFLAGS = -O2 -W -Wall

opt:	GLOBAL_CFLAGS=-O2 -W -Wall -fPIC
	
dbg:	GLOBAL_CFLAGS=-g3 -W -Wall -fPIC

opt:	TOTCalib

dbg:	TOTCalib
	
all:	TOTCalib

TOTCalib: TOTCalib.o TOTCalibDict.o
	g++ -shared -o libTOTCalib.so TOTCalib.o TOTCalibDict.o $(GLOBAL_CFLAGS) $(INCLUDE) ${LIBS} -lMathMore
	
TOTCalib.o:	TOTCalib.cpp TOTCalib.h
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c TOTCalib.cpp

TOTCalibDict.cpp: TOTCalib.h TOTCalib_LinkDef.h
	@echo "Generating dictionary $@ ..."
	$(ROOTCINT) -f $@ -c $^

TOTCalibDict.o:	TOTCalibDict.cpp
	g++ $(GLOBAL_CFLAGS) $(INCLUDE) -c TOTCalibDict.cpp

clean:	
	rm -f TOTCalib.o TOTCalibDict.o libTOTCalib.so TOTCalibDict.cpp
