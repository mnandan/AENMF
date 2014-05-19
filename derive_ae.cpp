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
#include "derive_ae.h"
#include <iostream>

using namespace std;
void DeriveAEW::getW() {
	UINT maxNormInd = 0;
	double maxDistVal = -INF;
	UINT rpInd = 0;
	for (UINT vI = 0; vI < N; vI++) {
		// find max norm vect: x1
		if (X[vI].nrm > maxDistVal) {
			maxNormInd = vI;
			maxDistVal = X[vI].nrm;
		}
	}
	dat.swapX(rpInd, maxNormInd);
	rpInd++;
	// find max norm from x1: x2
	std::vector<DistDat> dVals;
	dVals.resize(N);
	UINT maxNormInd2 = 0;
	double maxDistVal2 = -INF;
	for (UINT ind = rpInd; ind < N; ind++) {
		dVals[ind].ind = ind;
		if(X[ind].nrm > 0) {
			double curr_dist = X[ind].nrm + maxDistVal - \
					2 * getDotProduct(X[0].F,X[0].fSize, X[ind].F, X[ind].fSize);
			if (curr_dist > maxDistVal2) {
				maxNormInd2 = ind;
				maxDistVal2 = curr_dist ;
			}
			dVals[ind].dist = curr_dist;
		}
		else
			dVals[ind].dist = -INF;
	}
	dat.swapX(rpInd, maxNormInd2);
	updateCacheX(0);
	updateCacheX(1);
	dVals[maxNormInd2].dist = dVals[rpInd].dist;
	// find remaining R - rpInd vectors based on distance from polygon formed by
	// the selected rpInd vectors
	// using H matrix to store dot products temporarily
	for (UINT ind = 2; ind < N; ind++) {
		dat.putHval(ind, 0, getDotProduct(X[ind].F, X[ind].fSize, X[0].F, X[0].fSize));
		dat.putHval(ind, 1, getDotProduct(X[ind].F, X[ind].fSize, X[1].F, X[1].fSize));
	}
	sort(dVals.begin() + 2, dVals.end(), dValComp);
	double minDistVal;
	for (rpInd = 2; rpInd < std::min(N,R); rpInd++) {
		maxNormInd = 0;
		UINT dInd;
		rpSize = rpInd;
		maxDistVal = -INF;
		minDistVal = INF;
		for (dInd = rpInd; dInd < N; dInd++) {
			UINT ind = dVals[dInd].ind;
			if(ind < rpInd)
				continue;
			if(maxDistVal - dVals[dInd].dist > FLOAT_ZERO && maxNormInd != 0)
				break;
			//ind order is different due to sorting
			UINT hInd = origInd[ind];
			double rep_err = getRepErr(X[ind].nrm, hInd);
			dVals[dInd].dist = rep_err;
			if (rep_err > maxDistVal) {
				maxNormInd = ind;
				maxDistVal = rep_err;
				maxNormInd2 = dInd;
			} else if (rep_err < minDistVal) {
				minDistVal = rep_err;
			}
		}
		dat.swapX(rpInd, maxNormInd);
//		std::cout<<rpInd<<"\t"<<X[rpInd].index<<"\t"<<maxDistVal<<std::endl;
		dVals[maxNormInd2].dist = INF;
		std::swap(origInd[rpInd], origInd[maxNormInd]);
		if(dInd < N && minDistVal > FLOAT_ZERO) {
			std::vector<DistDat>::iterator t = dVals.begin() + dInd;
			std::vector<DistDat>::iterator e = dVals.end();
			for(; t != e; t++) {
				if(minDistVal - t->dist > FLOAT_ZERO)
					break;
			}
			sort(dVals.begin() + rpInd + 1, t, dValComp);
		} else {
			sort(dVals.begin() + rpInd + 1, dVals.end(), dValComp);
		}
		if(rpInd < std::min(N,R) - 1) {
			updateCacheX(rpInd);
			for (UINT ind = rpInd + 1; ind < N; ind++)
				dat.putHval(origInd[ind], rpInd, \
						    getDotProduct(X[rpInd].F, X[rpInd].fSize, X[ind].F, X[ind].fSize));
		}
	}
	dat.copyXW();
	dat.clearH();
}

void DeriveAEW::updateLambdaStat(UINT i) {
	if (lambda[i] >= LAMBDA_MAX)
		lambdaStat[i] = UPPER_BOUND;
	else if (lambda[i] <= LAMBDA_MIN)
		lambdaStat[i] = LOWER_BOUND;
	else
		lambdaStat[i] = FREE;
}

bool DeriveAEW::select_working_set(UINT &out_i, UINT &out_j) {
	double Gmax = -INF;
	double Gmax2 = -INF;
	int Gmax_idx = -1;
	int Gmin_idx = -1;
	double obj_diff_min = INF;

	for (UINT t = 0; t < rpSize; t++)
		if (!isUpperBound(t))
			if (-G[t] >= Gmax) {
				Gmax = -G[t];
				Gmax_idx = t;
			}

	int i = Gmax_idx;

	for (UINT j = 0; j < rpSize; j++) {
		if (!isLowerBound(j)) {
			double grad_diff = Gmax + G[j];
			if (G[j] >= Gmax2)
				Gmax2 = G[j];
			if (grad_diff > 0) {
				double obj_diff;
				double quad_coef = rpCache[i][i] + rpCache[j][j] \
				                    - 2*rpCache[i][j];
				obj_diff = -(grad_diff * grad_diff) / quad_coef;
				if (obj_diff <= obj_diff_min) {
					Gmin_idx = j;
					obj_diff_min = obj_diff;
				}
			}
		}
	}

	if (Gmax + Gmax2 < DERIVE_AE_EPS)		//converged
		return false;

	out_i = Gmax_idx;
	out_j = Gmin_idx;
	return true;
}

double DeriveAEW::getRepErr(double xNorm, UINT hInd) {
	//std::assert(rpSize != 0);
	double lInit = (double) 1.0 /(double) rpSize; //should sum to 1
	for (UINT i = 0; i < rpSize; i++) {
		G[i] = -H[hInd][i];
		lambda[i] = lInit;
		updateLambdaStat(i);
	}
	for (UINT i = 0; i < rpSize; i++) {
		for (UINT j = 0; j < rpSize; j++)
			G[j] += lambda[i] * rpCache[i][j];
	}
	// optimization step
	int iter = 0;
	while (iter < 1000) {
		UINT i, j;
		if (select_working_set(i, j) == false)
			break;
		++iter;
		// update lambda[i] and lambda[j], handle bounds carefully
		double oldLi = lambda[i];
		double oldLj = lambda[j];
		double quad_coef = rpCache[i][i] + rpCache[j][j] - 2 * rpCache[i][j];
		if (quad_coef <= 0)
			quad_coef = FLOAT_ZERO;
		double delta = (G[i] - G[j]) / quad_coef;
		double sum = lambda[i] + lambda[j];
		lambda[i] -= delta;
		lambda[j] += delta;
		if (lambda[i] > LAMBDA_MAX) {
			lambda[i] = LAMBDA_MAX;
			lambda[j] = sum - LAMBDA_MAX;
		} else if (lambda[j] > LAMBDA_MAX) {
			lambda[j] = LAMBDA_MAX;
			lambda[i] = sum - LAMBDA_MAX;
		}
		if (lambda[j] < LAMBDA_MIN) {
			lambda[j] = LAMBDA_MIN;
			lambda[i] = sum - LAMBDA_MIN;
		} else if (lambda[i] < LAMBDA_MIN) {
			lambda[i] = LAMBDA_MIN;
			lambda[j] = sum - LAMBDA_MIN;
		}
		// update G
		double delLi = lambda[i] - oldLi;
		double delLj = lambda[j] - oldLj;
		for (UINT k = 0; k < rpSize; k++)
			G[k] += (rpCache[k][i] * delLi + rpCache[k][j] * delLj);
		updateLambdaStat(i);
		updateLambdaStat(j);
	}
	double normSum = xNorm;
	for(UINT i = 0; i < rpSize; i++)
		if(lambda[i] != 0.0)
			normSum += (G[i] - H[hInd][i])*lambda[i];

	return normSum;
}
