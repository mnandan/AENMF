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
} distDatType;

struct feat_T {
	UINT fNum; // feature number
	double fVal;
};

struct dataVect_T {
	UINT numFeats;
	LABEL_T label;
	double B; //beta
	double A; //alpha
	double nrm; //norm
	feat_T * F; //features
};

class trainDat_T {
	dataVect_T * X;
	UINT totRpNum;
	UINT maxF; // largest index of features among all vectors
	double *w;
	UINT epochSize;
	UINT vectsRead;
	UINT totNumVects;
	UINT currReadVect;
	ofstream AEF;
	ofstream LBD;
	ofstream XDAT;
	char aeFile[1024];
public:
	double C;
	const dataVect_T * XC;
	const double *wC;

	trainDat_T(UINT maxF, UINT epochSize, char * aeFile, char *wFile, double C1 = 1) {
		this->maxF = maxF;
		w = new double[maxF + 1];

		this->epochSize = epochSize;
		X = new dataVect_T[epochSize];
		for (UINT ind = 0; ind < epochSize; ind++) {
			X[ind].F = NULL;
			X[ind].numFeats = 0;
		}

		C = C1;

		XC = X;
		wC = w;
		this->totNumVects = epochSize;
		totRpNum = 0;
		vectsRead = 0;
		currReadVect = 0;
		AEF.open("ae.dat");
		LBD.open("lbd.dat");
		XDAT.open("X.dat");
	}

	~trainDat_T() {
		AEF.close();
		LBD.close();
		XDAT.close();

		delete[] w;
		for (UINT ind = 0; ind < epochSize; ind++)
			if (X[ind].F != NULL)
				delete[] X[ind].F;
	}

	void addVector(vector<feat_T> &F, LABEL_T label, UINT numFeats) {
		UINT currInd = vectsRead % (epochSize - totRpNum) + totRpNum;

		if (X[currInd].F != NULL)
			delete[] X[currInd].F;

		X[currInd].label = label;
		X[currInd].numFeats = numFeats;
		X[currInd].B = 1.0;
		X[currInd].A = 0;
		X[currInd].nrm = 0;
		if (numFeats > 0) {
			X[currInd].F = new feat_T[numFeats];
			for (UINT i = 0; i < numFeats; i++) {
				X[currInd].F[i] = F[i];
				X[currInd].nrm += (F[i].fVal) * (F[i].fVal);
			}
		} else {
			X[currInd].F = NULL;
		}
		currReadVect = currInd + 1;
		vectsRead++;
	}

	void addVector2(vector<feat_T> &F, LABEL_T label, UINT numFeats, double BVal) {
		UINT currInd = vectsRead % (epochSize - totRpNum) + totRpNum;

		if (X[currInd].F != NULL)
			delete[] X[currInd].F;

		X[currInd].label = label;
		X[currInd].numFeats = numFeats;
		X[currInd].B = BVal;
		X[currInd].A = 0;
		X[currInd].nrm = 0;
		if (numFeats > 0) {
			X[currInd].F = new feat_T[numFeats];
			for (UINT i = 0; i < numFeats; i++) {
				X[currInd].F[i] = F[i];
				X[currInd].nrm += (F[i].fVal) * (F[i].fVal);
			}
		} else {
			X[currInd].F = NULL;
		}
		currReadVect = currInd + 1;
		vectsRead++;
	}

	void secondRead(){
		totRpNum = 0;
		vectsRead = 0;
	}

	void clearRead() {
	}

	UINT getCurrVectNum () {
		return currReadVect;
	}
	void swapX(UINT i, UINT j) {
		///// critical section
		swap(X[i], X[j]);
	}

	void swapRp(UINT i, UINT j) {
		///// critical section
		//X[j].B = b;
		swap(X[i], X[j]);
		++totRpNum;
	}

	void clearTrp() {
		totRpNum = 0;
	}

	void putAW(UINT vI, double alpha_new) {
		double delta = (alpha_new - X[vI].A) * (double) X[vI].label;
		X[vI].A = alpha_new;
		for (UINT fI = 0; fI < X[vI].numFeats; fI++) //feature index
			w[X[vI].F[fI].fNum] += delta * X[vI].F[fI].fVal;
	}

	void putA(UINT ind, double a) {
		X[ind].A = a;
	}
	void putX(UINT si, UINT p, dataVect_T *X2, distDatType *distVals) {
		for (UINT k = 0; k < p; k++) //mass swap
			X2[k] = X[distVals[k].ind];
		for (UINT k = 0; k < p; k++) {
			UINT j = k + si;
			X[j] = X2[k];
		}
	}

	void clearW() {
		for (UINT fI = 0; fI <= maxF; fI++)
			w[fI] = 0;
	}

	void writeAE(UINT ind) {
		if (AEF.is_open() && AEF.good()) {
			UINT numF = X[ind].numFeats;
			AEF<< (int)X[ind].label;
			for (UINT fl = 0; fl < numF; fl++)
				AEF <<" "<< X[ind].F[fl].fNum<<":"<< setprecision(8) << X[ind].F[fl].fVal;
			AEF<<endl;
		} else {
			cerr << "Cannot write to AE file\n";
			exit(1);
		}
	}

	void writeLBD(UINT M, double * lambda) {
		if (LBD.is_open() && LBD.good()) {
			LBD<< 1;
			for (UINT fl = 0; fl < M; fl++)
				if(lambda[fl] > 0)
					LBD <<" "<< fl + 1<<":"<< setprecision(8) << lambda[fl];
			LBD<<endl;
		} else {
			cerr << "Cannot write to Lambda file\n";
			exit(1);
		}
	}

	void writeX(UINT ind) {
		if (XDAT.is_open() && XDAT.good()) {
			UINT numF = X[ind].numFeats;
			XDAT<< (int)X[ind].label;
			for (UINT fl = 0; fl < numF; fl++)
				XDAT <<" "<< X[ind].F[fl].fNum<<":"<< setprecision(8) << X[ind].F[fl].fVal;
			XDAT<<endl;
		} else {
			cerr << "Cannot write to XDAT file\n";
			exit(1);
		}
	}


	void closeAEF() {
		AEF.close();
	}

	void openAEF() {
		AEF.open(aeFile);
		if(!(AEF.is_open() && AEF.good())) {
			cerr << "Cannot open AE file"<< endl;
			exit(1);
		}
	}

	void writeW(char * modelFname) {
		ofstream outMdl(modelFname);
		if (outMdl.is_open() && outMdl.good()) {
			outMdl << "solver_type L2R_L1LOSS_SVC_DUAL\n" << "nr_class 2\nlabel 1 -1\nnr_feature "
					<< maxF << "\n" << "bias -1\nw" << endl;
			for (UINT fI = 1; fI <= maxF; fI++)
				outMdl << setprecision(8) << w[fI] << endl;
		} else {
			cerr << "Cannot write to model file\n";
			exit(1);
		}
		outMdl.close();
	}

	void readW(char * modelFname) {
		ifstream inMdl(modelFname);
		string line;
		for(int i = 0; i < 6; i++) {
			if (inMdl.is_open() && inMdl.good()) {
				getline(inMdl, line);
			} else {
				cerr << "Cannot read model file\n";
				exit(1);
			}
		}

		for(UINT i = 1; i <= maxF; i++) {
			if (inMdl.is_open() && inMdl.good()) {
				getline(inMdl, line);
				w[i] = strtod((char *)line.c_str(),NULL);
			} else {
				cerr << "Cannot read model file\n";
				exit(1);
			}
		}
		inMdl.close();
	}

	UINT getXsize() {
		return min(vectsRead, epochSize);
	}

	UINT getEphSz() {
		return epochSize;
	}

	UINT getRpNum() {
		return totRpNum;
	}

	void wCopy(const double *wC) {
		for(UINT fI = 0; fI <= maxF; fI++)
			w[fI] = wC[fI];
	}

	double * makeWCopy() {
		double *wCopy = new double[maxF + 1];
		for(UINT fI = 0; fI <= maxF; fI++)
			wCopy[fI] = w[fI];
		return wCopy;
	}

	double wDiff(double *wOrig) {
		double diffNorm = 0;
		double wNorm = 0;
		for(UINT fI = 0; fI <= maxF; fI++) {
			double diff = wOrig[fI] - w[fI];
			wNorm += wOrig[fI]*wOrig[fI];
			diffNorm += diff*diff;
		}
		diffNorm = diffNorm/wNorm;
		return diffNorm;
	}
	
	/*double wDiff(double *wOrig) {
		double diffNorm = 0;
		double wFcnt = 0;
		for(UINT fI = 0; fI <= maxF; fI++) {
			double diff = wOrig[fI] - w[fI];
			if(diff != 0.0) {
				wFcnt++;
			}
			diffNorm += diff*diff;
			//diffNorm += diff;
		}
		diffNorm = diffNorm/wFcnt;
		return diffNorm;
	}*/

	const feat_T* getF(UINT vI) {
		return (const feat_T*) X[vI].F;
	}

};

#endif /* COMMON_H_ */
