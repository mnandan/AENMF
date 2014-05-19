/*
 * xstore.h
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#ifndef SPARSE_MAT_H_
#define SPARSE_MAT_H_

#include "common.h"
#include <vector>

struct FeatType {
	UINT fNum; // feature number
	double fVal;
};

typedef FeatType* FeatVect;

struct SparseVect {
	UINT index; // index of vector in input file
	UINT fSize;
	double nrm; // norm of vector
	FeatVect F; // features
	SparseVect(): index(0), fSize(0), nrm(0), F(NULL) {};
};

class SparseMat {
public:
	UINT size;
	std::vector <SparseVect> M;
	SparseMat(): size(0) {};
	~SparseMat() {
		for(UINT i = 0; i < size; i++)
			delete [] M[i].F;
	}
	void addVect(std::vector<FeatType> const &F, UINT featInd);
};

#endif /* SPARSE_MAT_H_ */
