/*
 * xstore.h
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#ifndef DENSE_MAT_H_
#define DENSE_MAT_H_

#include <vector>
#include "common.h"

typedef std::vector <double> DenseVect;

struct DenseMat {
	UINT size1, size2;
	std::vector <DenseVect>  M;
	DenseMat(): size1(0), size2(0) {}

	void put(UINT i, UINT j, double val) {
		//std::assert(size1 != 0 && size2 != 0);
		M[i][j] = val;
	}
	void init(UINT i, UINT j);
};

#endif /* DENSE_MAT_H_ */
