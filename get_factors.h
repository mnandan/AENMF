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
	std::vector <DenseVect> GM;
	DenseVect G, xTz, lambda;
	std::vector <UINT> index;
	double deriveHi(UINT hInd);
public:
	GetFact(AllData &dat_): DataHandler(dat_) {
// initialize data containers
		lambda.assign(D,0);    // 1xD vector
		GM.assign(R,lambda);    // RxD matrix
		G.assign(R,0);    // 1xR vector
		xTz.assign(R,0);    // 1xR vector
		index.assign(R,0);    // 1xR vector
	}
	double getH();
	double getW();
};


#endif /* GET_FACTORS_H_ */
