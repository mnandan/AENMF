/*
 * onePassDRS.h
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */

#ifndef ONEPASSDRS_H_
#define ONEPASSDRS_H_

#include "common.h"
#include "deriveAE.h"

class GetFact: public DeriveAE{
	const DataVect* X;
	TrainDat *trDat;
	UINT totVects,numRp,D;
	double * lbdArr;

	DataVect* X2;
	DistDat *distVals;
	UINT * clustEndInd;

	void segregateDataFLS2(UINT startIndex, UINT endIndex, UINT subSize);
    void segregateDataSLS(UINT startIndex, UINT endIndex, UINT V_VAL);
	double onePassDRSsub(UINT numRp, UINT startInd, UINT endInd);
	DistDat quickSelectDist(UINT n, UINT k);

	double getDotProductX(UINT i, UINT j) {
		double dotProduct = 0;
		UINT fn1 = 0, fn2 = 0;
		while (fn1 < X[i].numFeats && fn2 < X[j].numFeats) {
			if (X[i].F[fn1].fNum == X[j].F[fn2].fNum) {
				dotProduct += X[i].F[fn1].fVal * X[j].F[fn2].fVal;
				++fn2;
				++fn1;
			} else {
				if (X[i].F[fn1].fNum > X[j].F[fn2].fNum)
					++fn2;
				else
					++fn1;
			}
		}
		return dotProduct;
	}

	double getDotProductXW(UINT i, UINT j) {
		double dotProduct = 0;
		FeatType *F = X[i].F;
		for(UINT ind = 0; ind < X[i].numFeats; ind++) {
			auto t = trDat->W[j].find(F[ind].fNum);
			if ( t != trDat->W[j].end() )
				dotProduct += t->second *F[ind].fVal;
		}
		return dotProduct;
	}

	double getDotProductW(UINT i, UINT j) {
		double dotProduct = 0;
		if(trDat->W[i].size() <= trDat->W[j].size()) {
			auto t = trDat->W[i].begin(), e = trDat->W[i].end();
			for(; t!= e; t++) {
				auto wj = trDat->W[j].find(t->first);
				if ( wj != trDat->W[j].end() )
					dotProduct += wj->second * t->second;
			}
		}
		else {
			auto t = trDat->W[j].begin(), e = trDat->W[j].end();
			for(; t!= e; t++) {
				auto wi = trDat->W[i].find(t->first);
				if ( wi != trDat->W[i].end() )
					dotProduct += wi->second * t->second;
			}
		}
		return dotProduct;
	}

	void updateCacheW(UINT startInd, UINT indAdd) {
		rpCache[indAdd][indAdd] = getDotProductW(indAdd + startInd, indAdd + startInd);
		for (UINT ind = indAdd; ind > 0; ind--)
			rpCache[ind - 1][indAdd] = getDotProductW(indAdd + startInd, ind - 1 + startInd);
	}

	void updateCacheX(UINT startInd, UINT indAdd) {
		rpCache[indAdd][indAdd] = X[startInd + indAdd].nrm;
		for (UINT ind = indAdd; ind > 0; ind--)
			rpCache[ind - 1][indAdd] = getDotProductX(indAdd + startInd, ind - 1 + startInd);
	}

public:

	//double getFactors(UINT numRp, UINT startInd, UINT endInd, vector < vector< double> > &H);
	GetFact(TrainDat *trDat): DeriveAE(trDat->getR()) {
		this->trDat = trDat;
		X = trDat->XC;
		totVects = trDat->getN();
		numRp = trDat->getR();
		D = trDat->getD();
		X2 = NULL;
		distVals = NULL;
		clustEndInd = NULL;

		rpCache = new double*[numRp]; // size RxR
		for (UINT ind = 0; ind < numRp; ind++)
			rpCache[ind] = new double[numRp];
		xTz = new double[numRp];
		lbdArr = new double[numRp];
	}

	~GetFact() {
		for (UINT ind = 0; ind < numRp; ind++)
			delete[] rpCache[ind];
		delete[] rpCache;
		delete[] lbdArr;
		delete[] xTz;
	}

	void getWinit();
	void getH();
	void getW() { }
};


#endif /* ONEPASSDRS_H_ */
