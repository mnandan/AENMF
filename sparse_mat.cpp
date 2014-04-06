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
	for (UINT i = 0; i < featCnt; i++) {
		M[size].F.push_back(F[i]);
		double val = F[i].fVal;
		M[size].nrm += val*val;
	}
	size++;
}
