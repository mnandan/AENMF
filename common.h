/*
 * common.h
 *
 *  Created on: Aug 5, 2013
 *      Author: mnandan
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <string>

#define INF HUGE_VAL
#define FLOAT_ZERO 1e-9
#define LAMBDA_MAX 1.0
#define LAMBDA_MIN 0.0
#define LAMBDA_MIN_MF LAMBDA_MIN
#define DERIVE_AE_EPS 0.0001
#define GET_H_EPS 0.1
#define GET_W_EPS 0.1
#define REG_PARAM_H 0.0
#define REG_PARAM_W 0.0

typedef unsigned int UINT;

class errHand {
public:
	errHand(std::string msg) {
		throw msg;
	}
//
//	~errHand() {
//	}
};

#endif /* COMMON_H_ */
