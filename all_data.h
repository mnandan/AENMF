/*
 * data.h
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#ifndef ALL_DATA_H_
#define ALL_DATA_H_

#include <vector>
#include <algorithm>
#include "common.h"
#include "sparse_mat.h"
#include "dense_mat.h"

class AllData {
//data
	SparseMat X_;
	DenseMat W_;
	DenseMat H_;
//parameters
	UINT N_, D_, R_;
public:
	UINT const &N;
	UINT const &D;
	UINT const &R;

	DenseVect* const &WV;
	DenseVect* const &HV;
	std::vector <SparseVect> const &XV;

	AllData(UINT Rval):
		R_(Rval), N(N_), D(D_), R(R_), WV(W_.M), HV(H_.M), XV(X_.M) {
		N_ = 0;
		D_ = 0;
	}
	void swapX(UINT i, UINT j) {
		std::swap(X_.M[i], X_.M[j]);
	}
	void addXi(std::vector<FeatType> const &F, UINT featInd) {
		N_++;
		X_.addVect(F, featInd);
	}

	void initWH(UINT Dval) {
		D_ = Dval;
		H_.init(N,R);
		W_.init(R,D);
	}

	void putWval(UINT i, UINT j, double val) {
		W_.put(i,j,val);
	}
	void putHval(UINT i, UINT j, double val) {
		H_.put(i,j,val);
	}
	void copyXW();
	void putHvect(UINT i, double * lambda);
	void putWvect(UINT i, double * lambda);
	void clearW();
	void clearH();
	double calcNorm();
};


#endif /* ALL_DATA_H_ */

