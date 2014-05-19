/*
 * onePassDRS.h
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */
#ifndef GET_FACTORS_H_
#define GET_FACTORS_H_

#include "common.h"
#include "data_handler.h"
#include "all_data.h"
#include "dense_mat.h"

class GetFact: public DataHandler {
	// all data variables
	DenseVect* GM;
	DenseVect G, xTz;
	UINT* index;
	double deriveHi(UINT hInd);
public:
	GetFact(AllData &dat_): DataHandler(dat_) {
// initialize data containers
		G = new DenseType[R];    // 1xR vector
		xTz = new DenseType[R];    // 1xR vector
		index = new UINT[R];    // 1xR vector
		GM = new DenseVect[R];    // RxD matrix
		for(UINT i = 0; i < R; i++) {
			GM[i] = new DenseType[D];
//			G[i] = 0;
//			xTz[i] = 0;
//			index[i] = i;
//			for(UINT j = 0; j < D; j++)
//				GM[i][j] = 0;
		}
	}
	~GetFact() {
		delete [] G;
		delete [] xTz;
		delete [] index;
		for(UINT i = 0; i < R; i++)
			delete [] GM[i];
		delete [] GM;
	}
	void initWdata();
	double getH();
	double getW();
};


#endif /* GET_FACTORS_H_ */
