/*
 * common.h
 *
 *  Created on: Aug 5, 2013
 *      Author: mnandan
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <string>

#define TAU 1e-12
#define LAMBDA_MAX 1.0
#define LAMBDA_MIN 0.0
#define LAMBDA_MIN2 0.0
#define INF HUGE_VAL

#define REG_PARAM (double) 0.0
#define REG_PARAM2 (double) 0.0

typedef unsigned int UINT;

enum expVals {
	WRONG_PARAM,
	UNKNOWN
};

class errHand {
public:
	errHand(expVals e) {
		throw e;
	}

	errHand(std::string msg) {
		throw msg;
	}
//
//	~errHand() {
//	}
};

#endif /* COMMON_H_ */
