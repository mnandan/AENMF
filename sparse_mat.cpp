/*
 * xstore.cpp
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#include "sparse_mat.h"

void SparseMat::addVect(std::vector<FeatType> const &F, UINT featCnt) {
	M.push_back(SparseVect());
	M[size].index = size;
	M[size].nrm = 0;
	M[size].F = new FeatType[featCnt];
	M[size].fSize = featCnt;
	for (UINT i = 0; i < featCnt; i++) {
		M[size].F[i] = F[i];
		M[size].nrm += (F[i].fVal)*(F[i].fVal);
	}
	size++;
}
