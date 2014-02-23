/*
 * deriveAE.cpp
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "deriveAE.h"

using namespace std;
// refer LIBSVM for details. http://www.csie.ntu.edu.tw/~cjlin/libsvm/
// return 1 if already optimal, return 0 otherwise
int DeriveAE::select_working_set(int &out_i, int &out_j) {
	double Gmax = -INF;
	double Gmax2 = -INF;
	int Gmax_idx = -1;
	int Gmin_idx = -1;
	double obj_diff_min = INF;

	for (UINT t = 0; t < rpSize; t++)
		if (-G[t] >= Gmax) {
			Gmax = -G[t];
			Gmax_idx = t;
		}

	int i = Gmax_idx;

	for (UINT j = 0; j < rpSize; j++) {
		if (lambda_status[j]) {
			double grad_diff = Gmax + G[j];
			if (G[j] >= Gmax2)
				Gmax2 = G[j];
			if (grad_diff > 0) {
				double obj_diff;
				double quad_coef = rpCache[i][i] + rpCache[j][j]
						- 2 * rpCache[min(i, (int) j)][max(i, (int) j)];
				obj_diff = -(grad_diff * grad_diff) / quad_coef;
				if (obj_diff <= obj_diff_min) {
					Gmin_idx = j;
					obj_diff_min = obj_diff;
				}
			}
		}
	}

	if (Gmax + Gmax2 < 0.001)
		return 1;

	out_i = Gmax_idx;
	out_j = Gmin_idx;
	return 0;
}

void DeriveAE::deriveAE(UINT rpSize, double *lambda) {
	this->rpSize = rpSize;
	this->lambda = lambda;
	// initialize lambda_status
	{
		double lambdaInitVal = (double) 1.0 / rpSize; //should sum to 1
		for (UINT i = 0; i < rpSize; i++) {
			lambda[i] = lambdaInitVal;
			update_lambda_status(i);
		}
	}

	// initialize gradient
	{
		for (UINT i = 0; i < rpSize; i++) {
			G[i] = -xTz[i];
		}
		for (UINT i = 0; i < rpSize; i++) {
			if (lambda_status[i]) {
				for (UINT j = 0; j < rpSize; j++)
					G[j] += lambda[i] * rpCache[min(i, j)][max(i, j)];
			}
		}
	}
	// optimization step
	int iter = 0;
	while (iter < 1000) {
		int i, j;
		if (select_working_set(i, j) != 0)
			break;
		++iter;
		// update lambda[i] and lambda[j], handle bounds carefully
		double old_lambda_i = lambda[i];
		double old_lambda_j = lambda[j];
		double quad_coef = rpCache[i][i] + rpCache[j][j] - 2 * rpCache[min(i, j)][max(i, j)];
		if (quad_coef <= 0)
			quad_coef = TAU;
		double delta = (G[i] - G[j]) / quad_coef;
		double sum = lambda[i] + lambda[j];
		lambda[i] -= delta;
		lambda[j] += delta;
			if (lambda[j] < LAMBDA_MIN) {
				lambda[j] = LAMBDA_MIN;
				lambda[i] = sum - LAMBDA_MIN;
			}

			if (lambda[i] < LAMBDA_MIN) {
				lambda[i] = LAMBDA_MIN;
				lambda[j] = sum - LAMBDA_MIN;
			}
	
		// update G
		double delta_lambda_i = lambda[i] - old_lambda_i;
		double delta_lambda_j = lambda[j] - old_lambda_j;
		for (UINT k = 0; k < rpSize; k++)
			G[k] += (rpCache[min(i, (int) k)][max(i, (int) k)] * delta_lambda_i
					+ rpCache[min(j, (int) k)][max(j, (int) k)] * delta_lambda_j);
		update_lambda_status(i);
		update_lambda_status(j);
	}

}

double DeriveAE::getRepErr(UINT rpSize, double xNorm, double *lambda) {
	deriveAE(rpSize, lambda);
	double v1 = 0, v2 = 0, err;
	for (UINT i = 0; i < rpSize; i++) {
		if(lambda[i] != 0) {
			for (UINT j = 0; j < rpSize; j++)
				if(lambda[j] != 0)
					v1 += lambda[i] * lambda[j] * rpCache[min(i, j)][max(i, j)];
			v2 += lambda[i] * xTz[i];
		}
	}
	err = xNorm + v1 - 2 * v2;
	return err;
}
