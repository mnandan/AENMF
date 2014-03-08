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

//void DeriveAE::deriveAE(UINT rpSize, double *lambda, const DataVect * XC, int sI) {
//	this->rpSize = rpSize;
//	this->lambda = lambda;
//	double gradient = 0;
//	// initialize
//	double G, Gmax, Gmin;
//	double GmaxOld = INF, GminOld = -INF;
//	int activeSize = rpSize;
//	for (UINT i = 0; i < rpSize; i++)
//		lambda[i] = 0;
//
//	for(UINT i = 0; i <= D; i++)
//		w[i] = 0;
//	// optimization step
//	int iter = 0;
//	while (iter++ < 100) {
//		for (UINT i = 0; i < activeSize; i++) {
//			UINT j = i + rand() % (activeSize - i);
//			swap(index[i], index[j]);
//		}
//		Gmax = -INF;
//		Gmin = INF;
//		for (UINT s = 0; s < activeSize; s++) {
//			UINT vI = index[s]; //relative vector index
//			UINT aI = vI + sI; //actual vector index
//			gradient = -xTz[vI];
//			FeatType * F = XC[aI].F;
//			G = 0;
//			for (UINT fI = 0; fI < XC[aI].numFeats; fI++) //feature index
//				gradient += w[F[fI].fNum] * F[fI].fVal;
//
//			if (lambda[vI] == 0) {
//				if (gradient > GmaxOld) {
//					activeSize--;
//					swap(index[s], index[activeSize]);
//					s--;
//					continue;
//				} else if (gradient < 0)
//					G = gradient;
//			} else
//				G = gradient;
//
//			if (fabs(G) > TAU) {
//				double lambdaNew = max(lambda[vI] - G / rpCache[vI][vI], 0.0); //constraint
//				double delta = lambdaNew - lambda[vI];
//				lambda[vI] = lambdaNew;
//				for (UINT fI = 0; fI < XC[aI].numFeats; fI++) //feature index
//					w[F[fI].fNum] += delta * F[fI].fVal;
//			}
//			Gmax = max(Gmax, G);
//			Gmin = min(Gmin, G);
//		}
//		if (Gmax - Gmin <= 0.1) {
//			if (activeSize == rpSize)
//				break;
//			else {
//				activeSize = rpSize;
//				//cout << "*";
//				GmaxOld = INF;
//				GminOld = -INF;
//				continue;
//			}
//		}
//		GmaxOld = Gmax;
//		GminOld = Gmin;
//		if (GmaxOld <= 0)
//			GmaxOld = INF;
//		if (GminOld >= 0)
//			GminOld = -INF;
//	}
//}

void DeriveAE::deriveAE(UINT rpSize, double *lambda) {
	this->rpSize = rpSize;
	this->lambda = lambda;
	// initialize
	double PG, Gmax, Gmin;
	double GmaxOld = INF, GminOld = -INF;
	int activeSize = rpSize;
	for (UINT i = 0; i < rpSize; i++) {
		lambda[i] = 0;
		G[i] = -xTz[i];
	}

	// optimization step
	int iter = 0;
	while (iter++ < 1000) {
		for (UINT i = 0; i < activeSize; i++) {
			UINT j = i + rand() % (activeSize - i);
			swap(index[i], index[j]);
		}
		Gmax = -INF;
		Gmin = INF;
		for (UINT s = 0; s < activeSize; s++) {
			UINT vI = index[s]; //relative vector index
			PG = 0;
			if (lambda[vI] == 0) {
				if (G[vI] > GmaxOld) {
					activeSize--;
					swap(index[s], index[activeSize]);
					s--;
					continue;
				} else if (G[vI] < 0)
					PG = G[vI];
			} else
				PG = G[vI];

			if (fabs(PG) > TAU) {
				double lambdaNew = max(lambda[vI] - PG / rpCache[vI][vI], 0.0); //constraint
				double delta = lambdaNew - lambda[vI];
				lambda[vI] = lambdaNew;
				for (UINT k = 0; k < rpSize; k++)
					G[k] += rpCache[min(vI, (UINT) k)][max(vI, (UINT) k)] * delta;

			}
			Gmax = max(Gmax, PG);
			Gmin = min(Gmin, PG);
		}
		if (Gmax - Gmin <= 0.1) {
			if (activeSize == rpSize)
				break;
			else {
				activeSize = rpSize;
				//cout << "*";
				GmaxOld = INF;
				GminOld = -INF;
				continue;
			}
		}
		GmaxOld = Gmax;
		GminOld = Gmin;
		if (GmaxOld <= 0)
			GmaxOld = INF;
		if (GminOld >= 0)
			GminOld = -INF;
	}
}

int DeriveAE::select_working_set(int &out_i, int &out_j) {
	double Gmax = -INF;
	double Gmax2 = -INF;
	int Gmax_idx = -1;
	int Gmin_idx = -1;
	double obj_diff_min = INF;

	for (UINT t = 0; t < rpSize; t++)
		if (!is_upper_bound(t))
			if (-G[t] >= Gmax) {
				Gmax = -G[t];
				Gmax_idx = t;
			}

	int i = Gmax_idx;

	for (UINT j = 0; j < rpSize; j++) {
		if (!is_lower_bound(j)) {
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

void DeriveAE::deriveAE2(UINT rpSize, double *lambda) {
	this->rpSize = rpSize;
	// initialize lambdaStat
	this->lambda = lambda;
	{
		double alphaInitVal = (double) 1.0 / rpSize; //should sum to 1
		for (UINT i = 0; i < rpSize; i++) {
			lambda[i] = alphaInitVal;
			update_lambdaStat(i);
		}
	}

	// initialize gradient
	{
		for (UINT i = 0; i < rpSize; i++) {
			G[i] = -xTz[i];
		}
		for (UINT i = 0; i < rpSize; i++) {
			if (!is_lower_bound(i)) {
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
		double old_alpha_i = lambda[i];
		double old_alpha_j = lambda[j];
		double quad_coef = rpCache[i][i] + rpCache[j][j] - 2 * rpCache[min(i, j)][max(i, j)];
		if (quad_coef <= 0)
			quad_coef = TAU;
		double delta = (G[i] - G[j]) / quad_coef;
		double sum = lambda[i] + lambda[j];
		lambda[i] -= delta;
		lambda[j] += delta;
		if (sum > LAMBDA_MAX) {
			if (lambda[i] > LAMBDA_MAX) {
				lambda[i] = LAMBDA_MAX;
				lambda[j] = sum - LAMBDA_MAX;
			}
		} else {
			if (lambda[j] < LAMBDA_MIN) {
				lambda[j] = LAMBDA_MIN;
				lambda[i] = sum - LAMBDA_MIN;
			}
		}
		if (sum > LAMBDA_MAX) {
			if (lambda[j] > LAMBDA_MAX) {
				lambda[j] = LAMBDA_MAX;
				lambda[i] = sum - LAMBDA_MAX;
			}
		} else {
			if (lambda[i] < LAMBDA_MIN) {
				lambda[i] = LAMBDA_MIN;
				lambda[j] = sum - LAMBDA_MIN;
			}
		}
		// update G
		double delta_alpha_i = lambda[i] - old_alpha_i;
		double delta_alpha_j = lambda[j] - old_alpha_j;
		for (UINT k = 0; k < rpSize; k++)
			G[k] += (rpCache[min(i, (int) k)][max(i, (int) k)] * delta_alpha_i
					+ rpCache[min(j, (int) k)][max(j, (int) k)] * delta_alpha_j);
		update_lambdaStat(i);
		update_lambdaStat(j);
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
