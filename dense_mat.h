/*
 * xstore.h
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#ifndef DENSE_MAT_H_
#define DENSE_MAT_H_

#include "common.h"
#include <cassert>

typedef double DenseType;
typedef DenseType* DenseVect;

class DenseMat {
public:
	UINT size1, size2;
	DenseVect* M;
	DenseMat(): size1(0), size2(0), M(NULL) {};
	~DenseMat() {
		for(UINT i = 0; i < size1; i++)
			delete M[i];
		delete [] M;
	}
	void put(UINT i, UINT j, double val) {
		assert(((i < size1) && (j < size2)));
		M[i][j] = val;
	}
	void init(UINT i, UINT j);
};

#endif /* DENSE_MAT_H_ */
