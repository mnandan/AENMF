#include <iostream>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "file_int.h"

using namespace std;

FileInt::FileInt(InpParams const &pars, AllData &dat_) :dat(dat_) {
	xIN.open(pars.xFile.c_str());
	wOUT.open(pars.wFile.c_str());
	hOUT.open(pars.hFile.c_str());;
	maxNNZ = 0;
	N = 0;
	maxD = 0;
// initialize delimiters to read data in LIBSVM format
	delims1 = new bool[128];
	for (int i = 0; i < 128; i++)
		delims1[i] = false;
	delims1[' '] = true;
	delims1['\t'] = true;

	delims2 = new bool[128];
	for (int i = 0; i < 128; i++)
		delims2[i] = false;
	delims2[' '] = true;
	delims2['\t'] = true;
	delims2[':'] = true;
}


void FileInt::readXFile() {
	if(! (xIN.is_open() && xIN.good())) {
		std::string msg = "Cannot open input X file.";
		errHand e(msg);
	}
	while (procLine() == true) {
		N++;
	}
}

bool FileInt:: procLine() {
	string line;
	if (xIN.good() && getline(xIN, line)) {
		UINT lineSz = (UINT)line.length();
		char * lineC = (char *)line.c_str();	// C array representation
		char *p = lineC;    // Pointer used to point to start of substrings
		char *endPtr;    // for error checking
		double fVal;
		UINT fNum, featInd = 0, cInd = 0;
		FeatType tempF;
		//remove initial spaces before label
		for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);
		//find label characters
		for(;!delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);
		++cInd; // skip label

		while (cInd < lineSz) {
			//remove initial spaces before feature number
			try {
				for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);
			} catch (...){
				errHand e("Error in input file read");
			}
			p = lineC + cInd;	//start of sub-string
			//find feature number
			try {
				for(;!delims2[(int)lineC[cInd]] && cInd < lineSz; ++cInd);
			} catch (...){
				errHand e("Error in input file read");
			}
			lineC[cInd] = 0;	//end of sub-string
			//convert feature number substring to int
			fNum = (UINT) strtoll(p, &endPtr, 10);
			if(maxD < fNum)
				maxD = fNum;
			// LIBSVM format features are numbered from 1, need it to be from 0
			fNum -= 1;
			if (endPtr == p || *endPtr != '\0')
				errHand e("Error in input file read");
			++cInd;
			//remove initial spaces before feature value
			try {
				for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);
			} catch (...){
				errHand e("Error in input file read");
			}
			p = lineC + cInd;
			//find feature value
			try {
				for(;!delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);
			} catch (...){
				errHand e("Error in input file read");
			}
			lineC[cInd] = 0;
			fVal = strtod(p, &endPtr);    //convert feature value to double
			if (endPtr == p || (*endPtr != '\0' && !isspace(*endPtr)))
				errHand e("Error in input file read");
			++cInd;
			if(featInd < maxNNZ) {
				Ftemp[featInd].fNum = fNum;
				Ftemp[featInd].fVal = fVal;
			} else {
				tempF.fNum = fNum;
				tempF.fVal = fVal;
				Ftemp.push_back(tempF);
				++maxNNZ;
			}
			++featInd;
		}
		dat.addXi(Ftemp, featInd);
	}

	if(xIN.good())
		return true;
	else
		return false;
}

void FileInt::writeW() {
	if (wOUT.is_open() && wOUT.good()) {
		wOUT <<setprecision(8);
		for(UINT ind = 0; ind < dat.R; ind++) {
			wOUT << '1';
			for (UINT k = 0; k < dat.D; k++)
					if(dat.WV[ind][k] != 0)
						wOUT <<" "<< k + 1<<':'<< dat.WV[ind][k];
			wOUT<<endl;
		}
	} else
		errHand e("Cannot write matrix factor W to file\n");
}

void FileInt::writeH() {
	if (hOUT.is_open() && hOUT.good()) {
		hOUT <<setprecision(8);
		for(UINT ind = 0; ind < N; ind++) {
			hOUT << '1';
			for (UINT k = 0; k < dat.R; k++)
				if(dat.HV[ind][k] != 0)
					hOUT <<" "<< k + 1<<':'<<dat.HV[ind][k];
			hOUT<<endl;
		}
	} else
		errHand e("Cannot write matrix factor H to file\n");
}
