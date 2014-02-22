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
#define LAMBDA_MAX 1.0
#define LAMBDA_MIN 0

class DeriveAE_T {
	double *G; // gradient of objective function
	UINT rpSize;
	enum {
		LOWER_BOUND, UPPER_BOUND, FREE
	};
	char *alpha_status; // LOWER_BOUND, UPPER_BOUND, FREE
	double *alpha;
	UINT totVects;
    //UINT M_VAL;

	void update_alpha_status(int i) {
		if (alpha[i] >= LAMBDA_MAX)
			alpha_status[i] = UPPER_BOUND;
		else if (alpha[i] <= LAMBDA_MIN)
			alpha_status[i] = LOWER_BOUND;
		else
			alpha_status[i] = FREE;
	}

	bool is_upper_bound(int i) {
		return alpha_status[i] == UPPER_BOUND;
	}
	bool is_lower_bound(int i) {
		return alpha_status[i] == LOWER_BOUND;
	}
	bool is_free(int i) {
		return alpha_status[i] == FREE;
	}
	int select_working_set(int &i, int &j);

protected:
	double *xTz;
	double ** rpCache;

public:
	DeriveAE_T(UINT numRp) {
		G = new double[numRp];
		alpha_status = new char[numRp];
		xTz = NULL;
		rpCache = NULL;
		alpha = NULL;
        //this->M_VAL = M_VAL;
	}

	~DeriveAE_T() {
		delete[] G;
		delete[] alpha_status;
	}
	void deriveAE(UINT rpSize, double *lambda);
	double getRepErr(UINT rpSize, double xNorm, double *lambda);
};


#endif /* DERIVEAE_H_ */
