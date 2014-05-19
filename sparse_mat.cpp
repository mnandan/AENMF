/*
 * xstore.cpp
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#include "sparse_mat.h"

void SparseMat::addVect(std::vector<FeatType> const &F) {
	M.push_back(SparseVect());
	M[size].index = size;
	M[size].nrm = 0;
	UINT featCnt = F.size();
	M[size].F = new FeatType[featCnt];
	M[size].fSize = featCnt;
	std::vector<FeatType>::const_iterator t = F.begin(), e = F.end();
	UINT ind = 0;
	for (;t != e; t++) {
		M[size].F[ind].fNum = t->fNum;
		M[size].F[ind++].fVal = t->fVal;
		M[size].nrm += (t->fVal)*(t->fVal);
	}
	size++;
}
