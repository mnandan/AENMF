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
#include <unordered_map>
#include <algorithm>

using namespace std;

typedef unsigned int UINT;
typedef char LABEL_T;
#define M_VAL (UINT)2
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
	UINT index; // index of vector in input file
	double nrm; // norm of vector
	FeatType * F; // features
};

class TrainDat {
	DataVect * X;
	UINT D; // largest index of features among all vectors
	UINT R;
	ofstream hOUT, wOUT;
	char aeFile[1024];
	UINT vectNum;
	UINT totNumVects;
public:
	const DataVect * XC;
    vector <unordered_map < UINT, double> > W;	//key = fNum, val = fVal
	vector < vector< double> > H;	//factor H
	TrainDat(UINT D, UINT N, UINT R, char *hFname, char *wFname) {
		this->D = D;
		totNumVects = N;
		this->R = R;
		X = new DataVect[totNumVects];
		for (UINT ind = 0; ind < totNumVects; ind++) {
			X[ind].F = NULL;
			X[ind].numFeats = 0;
			H.push_back(vector <double> ());	//push empty vector
		}
		XC = X;
		vectNum = 0;

		for (UINT ind = 0; ind < R; ind++) {
			unordered_map < UINT, double> temp;
			W.push_back(temp);
		}

		hOUT.open(hFname);
		wOUT.open(wFname);
	}

	~TrainDat() {
		hOUT.close();
		wOUT.close();
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

	void copyXW() {	//initialize W matrix with the first numRp vector in X
		for(UINT i = 0; i<R; i++)
			for(UINT j = 0; j < X[i].numFeats; j++)
				W[i][X[i].F[j].fNum] = X[i].F[j].fVal;
		UINT i = 0;
		while(i < totNumVects) {	//Put X in original order
		    if(i == X[i].index)
		    	i++;
		    else
		    	swap(X[i], X[X[i].index]);
		}
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

	void writeWH() {
		for(UINT ind = 0; ind < totNumVects; ind++) {
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
		for(UINT ind = 0; ind < R; ind++) {
			if (wOUT.is_open() && wOUT.good()) {
				wOUT << '1';
				UINT size = W[ind].size();
				if(size > 0) {
					//sort each W[ind] hash before writing
					vector <UINT> fNums;
					auto t = W[ind].begin(), e = W[ind].end();
					for(; t!=e; t++)
						if(t->second > 0.0)
							fNums.push_back(t->first);
					sort(fNums.begin(),fNums.end());
					size = fNums.size();
					for (UINT i = 0; i< size;i++) {
						UINT fNum = fNums[i];
						wOUT <<" "<< fNum <<':'<<setprecision(8) << W[ind][fNum];
					}
				}
				wOUT<<endl;
			} else {
				cerr << "Cannot write matrix factor W to file\n";
				exit(1);
			}
		}
	}

	UINT getD() {
		return D;
	}

	UINT getN() {
		return totNumVects;
	}

	UINT getR() {
		return R;
	}
};

#endif /* COMMON_H_ */
