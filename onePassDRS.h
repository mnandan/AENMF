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
#include "aesvmSolver.h"

class onePDRS_T: public DeriveAE_T, public AesvmSolver_T {
	const dataVect_T* X;
	dataVect_T* X2;
	distDatType *distVals;
	UINT * clustEndInd;
	double ** lbdArr;
	UINT totVects;
	UINT epochSize;


	void segregateDataFLS2(UINT startIndex, UINT endIndex, UINT subSize);
    void segregateDataSLS(UINT startIndex, UINT endIndex, UINT V_VAL);
	double onePassDRSsub(UINT numRp, UINT startInd, UINT endInd);
	distDatType quickSelectDist(UINT n, UINT k);

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
	onePDRS_T( trainDat_T *trDat, UINT totVects, UINT numRp): DeriveAE_T(numRp), AesvmSolver_T(trDat, totVects) {
		X2 = NULL;
		distVals = NULL;
		clustEndInd = NULL;
		lbdArr = NULL;
		X = trDat->XC;
		this->epochSize = totVects;
		this->totVects = totVects;
        //this->M_VAL = M_VAL;
	}

	UINT getTotVects() {
		return totVects;
	}

	UINT getEpochSize() {
		return epochSize;
	}

	void clearWRp() {
		trDat->clearW();
		trDat->clearTrp();
	}

	double onePassDRS(UINT numRp, bool isFirstRun);
};


#endif /* ONEPASSDRS_H_ */
