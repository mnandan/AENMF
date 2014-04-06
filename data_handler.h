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
	std::vector <DenseVect>  const &W;
	std::vector <DenseVect>  const &H;
	std::vector <DenseVect> rpCache;
//inline functions
	void updateCacheW();
	void updateCacheX(UINT indAdd);
	double getDotProduct(DenseVect const &d, DenseVect const &e);
	double getDotProduct(FeatVect const &d, FeatVect const &e);
	double getDotProduct(FeatVect const &d, DenseVect const &e);
	double getDotProduct(DenseVect const &e, FeatVect const &d);
public:
	DataHandler(AllData &d1): dat(d1), X(d1.XV), W(d1.WV), H(d1.HV){
		//std::assert(dat.D > 0);
// copy references from dat
		R = dat.R;
		N = dat.N;
		D = dat.D;
		DenseVect temp;
		temp.assign(R,0);
		rpCache.assign(R,temp);    // RxR matrix
	}
};

double inline DataHandler::getDotProduct(DenseVect const &d, DenseVect const &e) {
	double dotProduct = 0;
	UINT fSize = d.size();
//		std::assert(d.size() == e.size());
	for(UINT fInd = 0; fInd < fSize; fInd++)
		if (d[fInd] != 0 && e[fInd] != 0)
			dotProduct += d[fInd]*e[fInd];
	return dotProduct;
}

double inline DataHandler::getDotProduct(FeatVect const &d, FeatVect const &e) {
	double dotProduct = 0;
	UINT dSize = d.size(), eSize = e.size();
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

double inline DataHandler::getDotProduct(FeatVect const &d, DenseVect const &e) {
	double dotProduct = 0;
	UINT dSize = d.size();
	//std::assert(d[dSize].fNum < e.size());
	for(UINT df = 0; df < dSize; df++) {
		double wfVal = e[d[df].fNum];
		if(wfVal != 0)
			dotProduct += wfVal*d[df].fVal;
	}
	return dotProduct;
}

double inline DataHandler::getDotProduct(DenseVect const &e, FeatVect const &d) {
	return getDotProduct(d,e);
}

void inline DataHandler::updateCacheX(UINT indAdd) {
	rpCache[indAdd][indAdd] = dat.XV[indAdd].nrm;
	for (UINT i = indAdd; i > 0; i--) {
		UINT ind = i;
		rpCache[ind][indAdd] = getDotProduct(dat.XV[indAdd].F,dat.XV[ind].F);
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
