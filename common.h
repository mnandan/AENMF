/*
 * common.h
 *
 *  Created on: Aug 5, 2013
 *      Author: mnandan
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>

using namespace std;

typedef unsigned int UINT;
typedef char LABEL_T;
#define M_VAL (UINT)10
#define INF HUGE_VAL

typedef struct {
	double dist;
	UINT ind;
} DistDat;

struct FeatType {
	UINT fNum; // feature number
	double fVal;
};

struct DataVect {
	UINT numFeats;
	int index; // index of vector in input file
	double nrm; // norm of vector
	FeatType * F; // features
};

class TrainDat {
	DataVect * X;
	UINT D; // largest index of features among all vectors
	ofstream wOUT;
	ofstream hOUT;
	char aeFile[1024];
	UINT vectNum;

public:
	const DataVect * XC;
	UINT totNumVects;

	TrainDat(UINT D, UINT totNumVects, char * wFname, char *hFname) {
		this->D = D;
		this->totNumVects = totNumVects;
		X = new DataVect[totNumVects];
		for (UINT ind = 0; ind < totNumVects; ind++) {
			X[ind].F = NULL;
			X[ind].numFeats = 0;
		}
		XC = X;
		vectNum = 0;
		wOUT.open(wFname);
		hOUT.open(hFname);
	}

	~TrainDat() {
		wOUT.close();
		hOUT.close();
		for (UINT ind = 0; ind < totNumVects; ind++)
			if (X[ind].F != NULL)
				delete[] X[ind].F;
	}

	void addVector(vector<FeatType> &F, UINT numFeats) {
		X[vectNum].numFeats = numFeats;
		X[vectNum].index = vectNum;
		X[vectNum].nrm = 0;
		if (numFeats > 0) {
			X[vectNum].F = new FeatType[numFeats];
			for (UINT i = 0; i < numFeats; i++) {
				X[vectNum].F[i] = F[i];
				X[vectNum].nrm += (F[i].fVal) * (F[i].fVal);
			}
		} else {
			X[vectNum].F = NULL;
		}
		vectNum++;
	}

	void swapX(UINT i, UINT j) {
		swap(X[i], X[j]);
	}

	void putX(UINT si, UINT p, DataVect *X2, DistDat *distVals) {
		for (UINT k = 0; k < p; k++) //mass swap
			X2[k] = X[distVals[k].ind];
		for (UINT k = 0; k < p; k++) {
			UINT j = k + si;
			X[j] = X2[k];
		}
	}

	void writeWH(UINT R, UINT N, vector < vector< double> > &H) {
		for(UINT ind = 0; ind < R; ind++) {
			if (wOUT.is_open() && wOUT.good()) {
				UINT numF = X[ind].numFeats;
				wOUT<< 1;
				for (UINT fl = 0; fl < numF; fl++)
					wOUT <<" "<< X[ind].F[fl].fNum<<":"<< setprecision(8) << X[ind].F[fl].fVal;
				wOUT<<endl;
			} else {
				cerr << "Cannot write matrix factor W to file\n";
				exit(1);
			}
			writeH(ind, H);
		}
		for(UINT ind = R; ind < N; ind++)
			writeH(ind, H);
	}

	void writeH(UINT ind, vector < vector< double> > &H) {
		if (hOUT.is_open() && hOUT.good()) {
			vector <double>::iterator t = H[ind].begin();
			hOUT << '1';
			for (UINT i = 1;t != H[ind].end(); t++, i++)
				if(*t != 0)
					hOUT <<" "<< i<<':'<<setprecision(8) << *t;
			hOUT<<endl;
		} else {
			cerr << "Cannot write matrix factor H to file\n";
			exit(1);
		}
	}
};

#endif /* COMMON_H_ */
