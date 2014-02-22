/*
 * common.h
 *
 *  Created on: Aug 3, 2013
 *      Author: mnandan
 */

#ifndef FILE_INT_H_
#define FILE_INT_H_

#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "common.h"

using namespace std;

class fHand_T {
	ifstream inpFH;
	UINT lineNum;
	bool *delims1, *delims2;
	trainDat_T *TrD;
	UINT maxNNZ;
	vector <feat_T> Ftemp;
	UINT numFeats;
	UINT totVects;
	UINT epochSize;
public:
	fHand_T(char *fName, UINT totVects, UINT maxNNZ, trainDat_T *TrD) {
		inpFH.open(fName);
		if(! (inpFH.is_open() && inpFH.good())) {
			cerr << "Cannot open file " << fName << endl;
			exit(1);
		}
		lineNum = 0;
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
		this->TrD = TrD;
		this->maxNNZ = maxNNZ;
		this->totVects = totVects;
		this->epochSize = totVects;
		Ftemp.resize(maxNNZ);
	}

	~fHand_T() {
		inpFH.close();
		delete [] delims1;
		delete [] delims2;
		vector<feat_T>().swap(Ftemp);
	}

	UINT getTotVects() {
		return totVects;
	}

	UINT getEpochSize() {
		return epochSize;
	}

	UINT getMaxNNZ() {
		return maxNNZ;
	}

	void secondRead(){
		lineNum = 0;
		inpFH.seekg(0, ios::beg);
		TrD->secondRead();
	}

	void clearRead() {
		TrD->clearRead();
	}
	bool procLine();
	bool procLine2();
};
#endif /* FILE_INT_H_ */
