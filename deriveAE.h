/*
 * deriveAE.h
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */

#ifndef DERIVEAE_H_
#define DERIVEAE_H_

#include "common.h"

#define TAU 1e-12
#define LAMBDA_MIN 0

class DeriveAE {
	double *G; // gradient of objective function
	UINT rpSize;
	bool *lambda_status;
	double *lambda;

	void update_lambda_status(int i) {
		if (lambda[i] <= LAMBDA_MIN)
			lambda_status[i] = false;
		else
			lambda_status[i] = true;
	}

	int select_working_set(int &i, int &j);

protected:
	double *xTz;
	double ** rpCache;

public:
	DeriveAE(UINT numRp) {
		G = new double[numRp];
		lambda_status = new bool[numRp];
		xTz = NULL;
		rpCache = NULL;
		lambda = NULL;
		rpSize = 0;
	}

	~DeriveAE() {
		delete[] G;
		delete[] lambda_status;
	}
	void deriveAE(UINT rpSize, double *lambda);
	double getRepErr(UINT rpSize, double xNorm, double *lambda);
};


#endif /* DERIVEAE_H_ */
