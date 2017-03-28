#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "TROOT.h"
#include "TString.h"

#include "MediPixWriteToEntuple.h"
#include "allpix_dm.h"


using namespace std;

void push_back_nbytes(unsigned int *, char *, int);

int main(int argc, char ** argv) {

	if(argc < 3) {
		cout << "use:" << endl;
		cout << argv[0] << "input_file(string) dataset(string)" << endl;
		exit(1);
	}

	//unsigned int nXYCWords = atoi(argv[1]);
	TString inputFile = argv[1];
	TString dataset = argv[2];

	WriteToNtuple MPXnTuple(dataset, "");
	FramesHandler frames(dataset);

	fstream filestr;
	filestr.open(inputFile, istream::in);

	unsigned int bytesRead = 0;
	char tempByte[4];
	unsigned int x = 0x0, y = 0x0, counts = 0x0;
	unsigned int prev_y = 0;
	//unsigned int zero1 = 0x0, zero2 = 0x0;
	unsigned int allcntr = 0;
	unsigned int nCountBytes = 4;

	while ( filestr.good() ) {

		// read X, 32 bits
		x = 0x0;
		filestr.get(tempByte[0]);
		filestr.get(tempByte[1]);
		filestr.get(tempByte[2]);
		filestr.get(tempByte[3]);
		push_back_nbytes(&x, tempByte, 4);

		//
		//zero1 = 0x0;
		//filestr.get(tempByte[0]);
		//filestr.get(tempByte[1]);
		//filestr.get(tempByte[2]);
		//filestr.get(tempByte[3]);
		//push_back_nbytes(&zero1, tempByte, 2);

		// read Y, 32 bits
		y = 0x0;
		filestr.get(tempByte[0]);
		filestr.get(tempByte[1]);
		filestr.get(tempByte[2]);
		filestr.get(tempByte[3]);
		push_back_nbytes(&y, tempByte, 4);

		//
		//zero2 = 0x0;
		//filestr.get(tempByte[0]);
		//filestr.get(tempByte[1]);
		//filestr.get(tempByte[2]);
		//filestr.get(tempByte[3]);
		//push_back_nbytes(&zero2, tempByte, 2);

		// Read Counts, 16 bits or 32 bits
		counts = 0x0;
		filestr.get(tempByte[0]);
		filestr.get(tempByte[1]);
		// if the next 16 bits are 0x0, the format is [X,Y,C]:=32,32,32 bits
		// if the next 16 bits are !0x0, then the format is [X,Y,C]:=32,32,16 bits
		filestr.get(tempByte[2]);
		filestr.get(tempByte[3]);
		nCountBytes = 4;

		if (tempByte[2] != 0x0 || tempByte[3] != 0x0) {
			// Counts format is 16 bits
			// rewind the last two lectures
			// cout << "C format is 16-bits ..." << endl;
			filestr.unget(); filestr.unget();
			nCountBytes = 2;
		}

		push_back_nbytes(&counts, tempByte, nCountBytes);

		//printf("%x, %x --> %x\n", x, y, counts);
		//printf("0x%x 0x%x, --> 0x%x | %d %d, --> %d\n", x, y, counts, x, y, counts);
		bytesRead += 10;

		/////////////////////////////////////////////////////////////////////////////////////////////
		// overflow recovery
		if(x > 255 || y > 255) {

			cout << "Coordinates overflow ! (possible data corruption) ..." <<
					"trying to recover by ignoring and searching next good position " << endl;

			// Rewind all last read (nCountBytes + 4 + 4) plus extra 32 bits
			for(unsigned int rew = 0 ; rew < nCountBytes + 4 + 4 + 4; rew++) { filestr.unget(); }

			// now search
			int zeroCount = 0;
			while ( zeroCount < 2 ) {
				filestr.get(tempByte[0]);
				filestr.get(tempByte[1]);
				printf("0x%x 0x%x  ", tempByte[0], tempByte[1]);
				if(tempByte[0] == 0x0 && tempByte[1] == 0x0) { zeroCount++; }
				else { zeroCount = 0; }
			}
			printf("\n");

			// searching next good position
			// tow consecutive > 0 pairs 0xXXXX and rewind 2 bytes
			bool searchN = true;
			unsigned int word1 = 0, word2 = 0;
			unsigned int last16Word;

			filestr.get(tempByte[0]);
			filestr.get(tempByte[1]);
			push_back_nbytes(&word1, tempByte, 2);
			last16Word = word1;

			while (searchN) {

				filestr.get(tempByte[0]);
				filestr.get(tempByte[1]);
				push_back_nbytes(&word2, tempByte, 2);

				if(last16Word != 0 && word2 != 0) { // good point rewind 2 and continue
					filestr.unget(); filestr.unget();
					searchN = false;
				}

				last16Word = word2;

			}
			//cout << "recovered ! attempt to continue reading after" << endl;
			//for(unsigned int rew = 0 ; rew < 6; rew++) filestr.get(tempByte[0]);

			continue; // jump to the beginning of the outer loop
		}
		/////////////////////////////////////////////////////////////////////////////////////////////

		// Before adding this pixel information check if it is time to store a frame
		if(y < prev_y) { // new frame
			frames.getFrameStructObject()->IncreaseId();
			frames.getFrameStructObject()->SetnX(256);
			frames.getFrameStructObject()->SetnY(256);
			MPXnTuple.fillVars(&frames);
		}
		prev_y = y;

		// Fill pixel information
		frames.getFrameStructObject()->FillOneElement(x, y, 256, counts);
		//printf("[%d] %d, %d --> %d\n", allcntr, x, y, counts);

		allcntr++;

	}

	MPXnTuple.closeNtuple();

	cout << "[ OK ] done with the conversion.  Check your output file --> "
		    << MPXnTuple.GetNtupleFileName() << endl;

}



void push_back_nbytes(unsigned int * val, char * bytes, int nbytes) {

	// indexes go like this
	// 0 --> lower byte  0x......XX
	// 3 --> higher byte 0xXX......

	*val &= 0x00000000;
	unsigned int tempVal;

	for(int idx = nbytes-1 ; idx >= 0 ; idx--) {

		// Get the byte
		tempVal = 0x0;
		tempVal ^= bytes[idx];
		// Clean up to have info for only one byte
		tempVal &= 0x000000ff;
		// Switch it to the right place
		for(int sw = 0 ; sw < idx ; sw++){
			tempVal = tempVal << 8;
		}
		// XOR the value
		*val ^= tempVal;

	}

}
