/*
 * deriveAE.h
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */
#ifndef DERIVE_AE_H_
#define DERIVE_AE_H_

#include <iostream>
#include "common.h"
#include "all_data.h"
#include "data_handler.h"
#include "dense_mat.h"

struct DistDat {
	double dist;
	UINT ind;
};

enum {
	LOWER_BOUND, UPPER_BOUND, FREE
};

class DeriveAEW: public DataHandler {
	UINT* lambdaStat;
	UINT* origInd;
	DenseVect G;
	DenseVect lambda;
	UINT rpSize;
	// to sort distance in getWinit
	static bool dValComp (DistDat i,DistDat j) {
		return (i.dist - j.dist > FLOAT_ZERO);
	}
	bool isUpperBound(UINT i) {
		return (lambdaStat[i] == UPPER_BOUND);
	}
	bool isLowerBound(UINT i) {
		return (lambdaStat[i] == LOWER_BOUND);
	}
	void updateLambdaStat(UINT i);
	bool select_working_set(UINT &out_i, UINT &out_j);
	double getRepErr(double xNorm, UINT hInd);
public:
	DeriveAEW(AllData &dat_): DataHandler(dat_) {
// initialize data containers
		lambda = new DenseType[R];    // 1xR vector
		G = new DenseType[R];    // 1xR vector
		lambdaStat = new UINT[R];
		rpSize = 0;
		origInd = new UINT[N];
		for(UINT i = 0; i < N; i++) {
			origInd[i] = i;
		}
	}
	~DeriveAEW() {
		delete [] lambda;
		delete [] G;
		delete [] lambdaStat;
		delete [] origInd;
	}
	void getW();
};

#endif /* DERIVE_AE_H_ */
