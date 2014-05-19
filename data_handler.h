/*
 * data_handler.h
 *
 *  Created on: Apr 6, 2014
 *      Author: zarcon
 */
#ifndef DATA_HANDLER_H_
#define DATA_HANDLER_H_

#include <vector>
#include "sparse_mat.h"
#include "dense_mat.h"
#include "all_data.h"

class DataHandler {
protected:
// all data variables
	AllData &dat;
	UINT N,R,D;
	std::vector <SparseVect> const &X;
	DenseVect * const &H;
	DenseVect * const &W;
	DenseVect * rpCache;
//inline functions
	void updateCacheW();
	void updateCacheX(UINT indAdd);
	double getDotProduct(DenseVect const &d, DenseVect const &e);
	double getDotProduct(FeatVect const &d, UINT dSize, FeatVect const &e, UINT eSize);
	double getDotProduct(FeatVect const &d, UINT dSize, DenseVect const &e);
public:
	DataHandler(AllData &d1): dat(d1), N(d1.N), R(d1.R), D(d1.D), X(d1.XV), H(d1.HV), W(d1.WV) {
		assert(D > 0);
		rpCache = new DenseVect[R];
		for(UINT i = 0; i < R; i++) {
			rpCache[i] = new DenseType[R];    // RxR matrix
//			for(UINT j = 0; j < R; j++)
//				rpCache[i][j] = 0;
		}
	}
};

double inline DataHandler::getDotProduct(DenseVect const &d, DenseVect const &e) {
	double dotProduct = 0;
	UINT fSize = D;
	for(UINT fInd = 0; fInd < fSize; fInd++)
		if (d[fInd] != 0 && e[fInd] != 0)
			dotProduct += d[fInd]*e[fInd];
	return dotProduct;
}

double inline DataHandler::getDotProduct(FeatVect const &d, UINT dSize, FeatVect const &e, UINT eSize) {
	double dotProduct = 0;
	UINT df = 0, ef = 0;
	while (df < dSize && ef < eSize) {
		if (d[df].fNum == e[ef].fNum) {
			dotProduct += d[df].fVal*e[ef].fVal;
			df++;
			ef++;
		} else if (d[df].fNum > e[ef].fNum)
			ef++;
		else
			df++;
	}
	return dotProduct;
}

double inline DataHandler::getDotProduct(FeatVect const &d, UINT dSize, DenseVect const &e) {
	double dotProduct = 0;
	for(UINT df = 0; df < dSize; df++) {
		assert(d[df].fNum < D);
		double wfVal = e[d[df].fNum];
		if(wfVal != 0)
			dotProduct += wfVal*d[df].fVal;
	}
	return dotProduct;
}

//double inline DataHandler::getDotProduct(DenseVect const &e, FeatVect const &d) {
//}

void inline DataHandler::updateCacheX(UINT indAdd) {
	rpCache[indAdd][indAdd] = dat.XV[indAdd].nrm;
	for (UINT i = indAdd; i > 0; i--) {
		UINT ind = i;
		rpCache[ind][indAdd] = getDotProduct(dat.XV[indAdd].F, dat.XV[indAdd].fSize, dat.XV[ind].F, dat.XV[ind].fSize);
		rpCache[indAdd][ind] = rpCache[ind][indAdd];
	}
}

void inline DataHandler::updateCacheW() {
	for (UINT ind1 = 0; ind1 < R; ind1++) {
		rpCache[ind1][ind1] = getDotProduct(dat.WV[ind1], dat.WV[ind1]);
		for (UINT ind2 = ind1; ind2 > 0; ind2--) {
			UINT ind2u = ind2 - 1;
			rpCache[ind2u][ind1] = getDotProduct(dat.WV[ind1], dat.WV[ind2u]);
			rpCache[ind1][ind2u] = rpCache[ind2u][ind1];
		}
	}
}
#endif /* DATA_HANDLER_H_ */
