/*
 * xstore.cpp
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#include "dense_mat.h"

void DenseMat::init(UINT i, UINT j) {
	DenseVect F;
	F.assign(j,0);
	M.assign(i,F);	// ixj matrix
	size1 = i;
	size2 = j;
}
