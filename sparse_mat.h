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

typedef std::vector <FeatType> FeatVect;

struct SparseVect {
	UINT index; // index of vector in input file
	double nrm; // norm of vector
	std::vector <FeatType> F; // features
};

struct SparseMat {
	std::vector <SparseVect> M;
	UINT size;
	SparseMat(): size(0){}
	void addVect(std::vector<FeatType> const &F, UINT featCnt);
};

#endif /* SPARSE_MAT_H_ */
