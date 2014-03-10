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
	DataVect *X2;
	DistDat * distVals;
	UINT *clustEndInd;

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

	void updateCacheX(UINT startInd, UINT indAdd) {
		rpCache[indAdd][indAdd] = X[startInd + indAdd].nrm;
		for (UINT ind = indAdd; ind > 0; ind--) {
			rpCache[ind - 1][indAdd] = getDotProductX(indAdd + startInd, ind - 1 + startInd);
			rpCache[indAdd][ind - 1] = rpCache[ind - 1][indAdd];
		}
	}
	DistDat quickSelectDist(UINT n, UINT k);
	void segregateDataFLS2(UINT startIndex, UINT endIndex, UINT P_VAL);
	void segregateDataSLS(UINT startIndex, UINT endIndex, UINT V_VAL);

public:
	GetFact(TrainDat *trDat): DeriveAE(trDat) {
		X2 = NULL;
		distVals = NULL;
		clustEndInd = NULL;
	}

	~GetFact() {

	}

	void getWinit();
	double getH();
	double getW();
};


#endif /* ONEPASSDRS_H_ */
