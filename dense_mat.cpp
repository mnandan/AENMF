/*
 * xstore.cpp
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#include "dense_mat.h"

void DenseMat::init(UINT i, UINT j) {
	M = new DenseVect [i];
	for(UINT i2 = 0; i2 < i; i2++) {
		M[i2] = new DenseType[j];
		for(UINT j2 = 0; j2 < j; j2++)
			M[i2][j2] = 0;
	}
	size1 = i;
	size2 = j;
}
