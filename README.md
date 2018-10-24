# MAFalda framework
Code contributors: John Idarraga (idarraga)(2009+), Thomas Remy Victor Billoud (tbilloud), Oliver Michael Keller, jsroux, Navrit Bal (navrit)(2018).

The CERN hosted version on SVN is not maintained and is an old version, do not use it: https://twiki.cern.ch/twiki/bin/view/Main/MAFalda

## Compile the framework
* The only dependency is the [ROOT](http://root.cern.ch/) framework which should be installed and properly setup.
* This is being tested and run with ROOT v6.14.04, the most recent version at the time of writing.
* Note: Make sure the desired version of ROOT is 'sourced' before compiling this framework, otherwise it will complain about not being able to find ROOT.

```
cd mafalda_framework
make clean
make -j$(getconf _NPROCESSORS_ONLN)
```


## Run the software
* If the compilation finishes successfully (return code 0), you're ready to go.
* For a first test try:
```
cd mafalda_framework
root -l runAfterClusteringExample.C
```

## The top layer
* The file runAfterClusteringExample.C is a C++ macro meant to be interpreted using ROOT/CINT. When you issue the command: root -l runAfterClusteringExample.C the top layer is interpreted. The entire library including user defined algorithms is compiled for the particular architecture. The only piece of code which is interpreted is the Top layer which should not contain any exhausting computations or too complex syntax. The advantage of having an interpreted top layer is the use of the different CINT capabilities and the facility to make small changes which can propagate to the defined objects, and run without recompiling. 
* Note: For batch runs you will probably want to give root the '-q' flag in order to finish the root process at the end of the MAFalda run. 

In the Top layer you ought to take at least the following 3 actions:
1. Load the mafalda library
```
gSystem->Load("libMediPixAnalysisCore.so");
```
2. Create and instance of the AnalysisManager and load at least one input file
```
AnalysisManager mpxAnalysis("testdata/MPXNtuple_241Am_106Ru_137Cs_13keV_30deg_784_100V.root");
```
* Loading multiple input files can be done by calling the constructor of AnalysisManager with a different argument. The argument should be a clear text file containing a list of the root files to be loaded 
```
AnalysisManager mpxAnalysis("list_of_files.txt");
```
* Inside list_of_files.txt you will have something like:
```
pathtomydata/file1.root
pathtomydata/file2.root
pathtomydata/file3.root
```
3. Create an instance of a particular algorithm and connect it to the run
```
AfterClusteringExample * ac = new AfterClusteringExample;
mpxAnalysis.ConnectAlgo("AfterClusteringExample", ac);
```
* At this point you can call any member of the Algorithm class through the pointer that you just created.
4. Start the run
```
mpxAnalysis.Run();
```
A few examples of top layers, run*.C come with the package. 


# Updating to latest version
```
git pull
```
or if you want to have the exact same as the master branch:
```
git pull --reset hard origin/master
```


# Convert my data to the multiframe ROOT-based format (compatible with USBPix and Pixelman, FE-I3, FE-I4 and TimePix, Medipix2, 3)
* There is a small set of programs that are able to convert your data in different input formats (Pixelman ascii, binary, USBPix, TurboDAQ, etc) to the ROOT-based format used by the framework. It is available in the directory mpxdataconverter

* If you have been able to compile the framework, compilation of the mpxdataconverter should work as well.
```
cd mpxdataconverter
make clean
make -j$(getconf _NPROCESSORS_ONLN)
```

1. For Dexter and Pixelman data run:
```
./mpxdataconverter
```
2.  For USBPix data run:
```
./raw_to_pixeldm
```
* Arguments are necessary. Read the message you get when you run it without arguments. 
* Once the program is done, you get a single file (in ROOT-base Multiframe format) containing all the frames the program found in the directory you passed as an argument. You're ready to work now :)

*  Now let's suppose that, as output, you get a file named "MPXNtuple_myrun.root". To run MAFalda on this data you need to modify the runExample.C top layer macro. Open runExample.C (coming with the framework), and in the line 
```
AnalysisManager mpxAnalysis("testdata/MPXNtuple_myrun.root");
```
replace the existing line by your new filename "MPXNtuple_myrun.root" (relative or full path). Run the framework again. At this point you might be looking forward to write your own algorithms for your data. Use the following script:
```
./newalgo.sh
```
follow the instructions and have fun writing your own algo. To speed up your development check what exists already within the framework so you don't re-invent the wheel. The algorithm called BlobsFinder for example, identifies Clusters in a frame and is useful in a great deal of cases. You are welcome to use it with your data. It comes with MAFalda by default. For new algorithm implementation visit the section WriteAMPXAlgo.
