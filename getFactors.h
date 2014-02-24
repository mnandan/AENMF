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
	DataVect* X2;
	DistDat *distVals;
	TrainDat *trDat;
	UINT * clustEndInd;
	double ** lbdArr;
	UINT totVects;

	void segregateDataFLS2(UINT startIndex, UINT endIndex, UINT subSize);
    void segregateDataSLS(UINT startIndex, UINT endIndex, UINT V_VAL);
	double onePassDRSsub(UINT numRp, UINT startInd, UINT endInd);
	DistDat quickSelectDist(UINT n, UINT k);

	double getDotProduct(UINT i, UINT j) {
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

	void updateCache(UINT startInd, UINT indAdd) {
		rpCache[indAdd][indAdd] = X[startInd + indAdd].nrm;
		for (UINT ind = indAdd; ind > 0; ind--)
			rpCache[ind - 1][indAdd] = getDotProduct(indAdd + startInd, ind - 1 + startInd);
	}

public:
	GetFact(TrainDat *trDat, UINT numRp): DeriveAE(numRp) {
		X2 = NULL;
		distVals = NULL;
		clustEndInd = NULL;
		lbdArr = NULL;
		this->trDat = trDat;
		X = trDat->XC;
		this->totVects = trDat->totNumVects;
	}

	double getFactors(UINT numRp, UINT startInd, UINT endInd, vector < vector< double> > &H);
	double getV(UINT numRp, UINT startInd, UINT endInd);
	double getH(UINT numRp, UINT startInd, UINT endInd, vector < vector< double> > &H);
};


#endif /* ONEPASSDRS_H_ */
