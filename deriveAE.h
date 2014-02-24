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

class DeriveAE {
	UINT rpSize;
	enum {
		LOWER_BOUND, UPPER_BOUND, FREE
	};
	char *lambdaStat; // LOWER_BOUND, UPPER_BOUND, FREE

	double *lambda;
	int *index;
	double *G;
	void update_lambdaStat(int i) {
		if (lambda[i] >= LAMBDA_MAX)
			lambdaStat[i] = UPPER_BOUND;
		else if (lambda[i] <= LAMBDA_MIN)
			lambdaStat[i] = LOWER_BOUND;
		else
			lambdaStat[i] = FREE;
	}

	bool is_upper_bound(int i) {
		return lambdaStat[i] == UPPER_BOUND;
	}
	bool is_lower_bound(int i) {
		return lambdaStat[i] == LOWER_BOUND;
	}
	bool is_free(int i) {
		return lambdaStat[i] == FREE;
	}
	int select_working_set(int &i, int &j);
protected:
	double *xTz;
	double ** rpCache;

public:
	DeriveAE(UINT numRp) {
		index = new int[numRp];
		G = new double[numRp];
		lambdaStat = new char[numRp];
		for(UINT i = 0; i< numRp; i++)
			index[i] = i;
		xTz = NULL;
		rpCache = NULL;
		lambda = NULL;
		rpSize = 0;
	}

	~DeriveAE() {
		delete[] index;
		delete[] lambdaStat;
		delete[]G;
	}
	void deriveAE(UINT rpSize, double *lambda);
	void deriveAE2(UINT rpSize, double *lambda);
	double getRepErr(UINT rpSize, double xNorm, double *lambda);
};


#endif /* DERIVEAE_H_ */
