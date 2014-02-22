#include "fileInt.h"

using namespace std;


bool fHand_T:: procLine() {
	string line;
	if (inpFH.good() && getline(inpFH, line)) {
		++lineNum;
		UINT lineSz = (UINT)line.length();
		char * lineC = (char *)line.c_str();
		char *p = lineC;

		UINT fNum;
		double fVal;
		feat_T tempF;
		char *endPtr;
		numFeats = 0;

		UINT cInd = 0;
		for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces

		p = lineC + cInd;
		for(;!delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find label characters

		lineC[cInd] = 0;
		LABEL_T label = (LABEL_T) strtol(p, &endPtr, 10);
		if(label <= 0)
			label = -1;
		else
			label = 1;
		if (endPtr == p || *endPtr != '\0') {
			cerr << "Error in input file read\n";
			exit(1);
		}
		++cInd;

		while (cInd < lineSz) {
			for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces

			p = lineC + cInd;
			for(;!delims2[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find fNum characters

			lineC[cInd] = 0;
			fNum = (UINT) strtoll(p, &endPtr, 10);
			if (endPtr == p || *endPtr != '\0') {
				cerr << "Error in input file read\n";
				exit(1);
			}
			++cInd;

			for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces

			p = lineC + cInd;
			for(;!delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find fVal characters
			lineC[cInd] = 0;
			fVal = strtod(p, &endPtr);
			if (endPtr == p || (*endPtr != '\0' && !isspace(*endPtr))) {
				cerr << "Error in input file read\n";
				exit(1);
			}
			if(numFeats < maxNNZ) {
				Ftemp[numFeats].fNum = fNum;
				Ftemp[numFeats].fVal = fVal;
			} else {
				tempF.fNum = fNum;
				tempF.fVal = fVal;
				Ftemp.push_back(tempF);
				++maxNNZ;
			}
			++numFeats;
			++cInd;
		}
		TrD->addVector(Ftemp, label, numFeats);
	}

	if(inpFH.good())
		return true;
	else
		return false;
}

bool fHand_T:: procLine2() {
	string line;
	if (inpFH.good() && getline(inpFH, line)) {
		++lineNum;
		UINT lineSz = (UINT)line.length();
		char * lineC = (char *)line.c_str();
		char *p = lineC;

		UINT fNum;
		double fVal;
		feat_T tempF;
		char *endPtr;
		numFeats = 0;

		UINT cInd = 0;
		for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces

		p = lineC + cInd;
		for(;!delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find label characters

		lineC[cInd] = 0;
		LABEL_T label = (LABEL_T) strtol(p, &endPtr, 10);
		if(label <= 0)
			label = -1;
		else
			label = 1;
		if (endPtr == p || *endPtr != '\0') {
			cerr << "Error in input file read\n";
			exit(1);
		}
		++cInd;

		for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces
		p = lineC + cInd;
		for(;!delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find B characters
		lineC[cInd] = 0;
		double BVal = strtod(p, &endPtr);
		if (endPtr == p || (*endPtr != '\0' && !isspace(*endPtr))) {
			cerr << "Error in input file read\n";
			exit(1);
		}
		++cInd;

		while (cInd < lineSz) {
			for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces

			p = lineC + cInd;
			for(;!delims2[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find fNum characters

			lineC[cInd] = 0;
			fNum = (UINT) strtoll(p, &endPtr, 10);
			if (endPtr == p || *endPtr != '\0') {
				cerr << "Error in input file read\n";
				exit(1);
			}
			++cInd;

			for(;delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//remove initial spaces

			p = lineC + cInd;
			for(;!delims1[(int)lineC[cInd]] && cInd < lineSz; ++cInd);		//find fVal characters
			lineC[cInd] = 0;
			fVal = strtod(p, &endPtr);
			if (endPtr == p || (*endPtr != '\0' && !isspace(*endPtr))) {
				cerr << "Error in input file read\n";
				exit(1);
			}
			if(numFeats < maxNNZ) {
				Ftemp[numFeats].fNum = fNum;
				Ftemp[numFeats].fVal = fVal;
			} else {
				tempF.fNum = fNum;
				tempF.fVal = fVal;
				Ftemp.push_back(tempF);
				++maxNNZ;
			}
			++numFeats;
			++cInd;
		}
		TrD->addVector2(Ftemp, label, numFeats, BVal);
	}

	if(inpFH.good())
		return true;
	else
		return false;
}
