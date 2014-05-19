/*
 * onePassDRS.cpp
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include "get_factors.h"


void GetFact::initWdata() {
	for (UINT ind = 0; ind < R; ind++) {
		index[ind] = ind;
		for (UINT ind2 = 0; ind2 < R; ind2++)
			rpCache[ind][ind2] = 0;
		for (UINT ind2 = 0; ind2 < D; ind2++)
			GM[ind][ind2] = 0;
	}
	UINT xFnum;
	double tempH;
	for (UINT ind = 0; ind < N; ind++) {
		for (UINT ind2 = 0; ind2 < R; ind2++) {
			double hVal = H[ind][ind2];
			rpCache[ind2][ind2] += hVal*hVal;
			for (UINT ind3 = ind2; ind3 > 0; ind3--) {
				UINT ind3u = ind3 - 1;
				tempH = hVal*H[ind][ind3u];
				rpCache[ind2][ind3u] += tempH;
				rpCache[ind3u][ind2] += tempH;
			}
		}
		FeatVect const &F = X[ind].F;
		xFnum = X[ind].fSize;
		if(ind < R) {
			for (UINT fInd = 0; fInd < xFnum; fInd++) {
				GM[ind][F[fInd].fNum] -= REG_PARAM_W*F[fInd].fVal;
				for (UINT ind2 = 0; ind2 < R; ind2++) {
					GM[ind2][F[fInd].fNum] -= H[ind][ind2]*F[fInd].fVal;
				}
			}
			if(REG_PARAM_W > 0)
				for (UINT fInd = 0; fInd < D; fInd++)
					GM[ind][fInd] += REG_PARAM_W*W[ind][fInd];
		} else {
			for (UINT fInd = 0; fInd < xFnum; fInd++)
				for (UINT ind2 = 0; ind2 < R; ind2++) {
					GM[ind2][F[fInd].fNum] -=H[ind][ind2]*F[fInd].fVal;
				}
		}
	}

	for (UINT ind = 0; ind < R; ind++)
		for (UINT fInd = 0; fInd < D; fInd++)
			for (UINT ind2 = 0; ind2 < R; ind2++)
				GM[ind][fInd] += rpCache[ind][ind2]*W[ind2][fInd];
}

double GetFact::getW( ) {
	initWdata();
	// optimization step
	double gPod, lambdaNew, delta, gradSum, prevGradSum = INF;
	UINT iter = 0;
	while (iter++ < 1000) {
		for (UINT i = 0; i < R; i++) {
			UINT j = i + rand() % (R - i);
			std::swap(index[i], index[j]);
		}
		gradSum = 0;
		for (UINT s = 0; s < R; s++) {
			UINT vI = index[s]; //vector index
			if(rpCache[vI][vI] != 0) {
				for (UINT fInd = 0; fInd < D; fInd++) {
					double denom = REG_PARAM_W + rpCache[vI][vI];
					gPod = GM[vI][fInd]/denom;
					lambdaNew = std::max(W[vI][fInd] - gPod, LAMBDA_MIN_MF);
					delta = lambdaNew - W[vI][fInd];
					if(fabs(delta) > FLOAT_ZERO) {
						dat.putWval(vI,fInd,lambdaNew);		//W[vI][fInd] = lambdaNew;
						GM[vI][fInd] += delta*REG_PARAM_W;
						for (UINT gI = 0; gI < R; gI++)
							GM[gI][fInd] += delta*rpCache[gI][vI];
					}
				}
			}
		}
		for (UINT s = 0; s < R; s++)
			for (UINT fInd = 0; fInd < D; fInd++)
				gradSum += GM[s][fInd]*GM[s][fInd];
		gradSum = sqrt(gradSum);
		if(fabs(prevGradSum - gradSum) < GET_W_EPS) {
			//cout << "deltaSum = "<< deltaSum <<endl;
			break;
		}
		prevGradSum = gradSum;
	}
	return gradSum;
}


double GetFact::getH() {
	updateCacheW();	//update rpCache
	for (UINT i = 0; i < R; i++)
		index[i] = i;
// derive lambda for other vectors
	double frobNormSq = 0;
	for (UINT ind = 0; ind < N; ind++) {
		for (UINT rpInd2 = 0; rpInd2 < R; rpInd2++)
			xTz[rpInd2] = getDotProduct(X[ind].F, X[ind].fSize, W[rpInd2]);
		double norm = deriveHi(ind);
//		std::cout<<norm<<"\t"<<X[ind].nrm<<"\n";
		frobNormSq += norm;
	}
	return sqrt(frobNormSq);
}

double GetFact::deriveHi(UINT hInd) {
	// initialize
	double PG, Gmax, Gmin;
	for (UINT i = 0; i < R; i++) {
		G[i] = -xTz[i];
		for (UINT k = 0; k < R; k++)
			G[i] += rpCache[i][k] * H[hInd][k];
	}
	// optimization step
	int iter = 0;
	while (iter++ < 1000) {
		for (UINT i = 0; i < R; i++) {
			UINT j = i + rand() % (R - i);
			std::swap(index[i], index[j]);
		}
		Gmax = -INF;
		Gmin = INF;
		for (UINT s = 0; s < R; s++) {
			UINT vI = index[s]; //relative vector index
			if (!(H[hInd][vI] == 0 && G[vI] >= 0)) {
				PG = G[vI];
				if (fabs(PG) > FLOAT_ZERO) {
					double gVal = PG/(REG_PARAM_H + rpCache[vI][vI]);
					double lambda2 = std::max(H[hInd][vI] - gVal, LAMBDA_MIN_MF);
					double delta = lambda2 - H[hInd][vI];
					if(delta != 0) {
						dat.putHval(hInd,vI,lambda2);
						for (UINT k = 0; k < R; k++)
							G[k] += rpCache[vI][k] * delta;
					}
				}
				Gmax = std::max(Gmax, PG);
				Gmin = std::min(Gmin, PG);
			}
		}
		if (Gmax - Gmin <= GET_H_EPS)
			break;
	}
	double normSum = X[hInd].nrm;
	for(UINT i = 0; i < R; i++) {
		if(H[hInd][i] != 0) {
			normSum += (G[i] - xTz[i])*H[hInd][i];
		}
	}
	return normSum;
}

//double GetFact::deriveHi(UINT hInd) {
//	// initialize
//	double PG, Gmax, Gmin;
//	double GmaxOld = INF, GminOld = -INF;
//	UINT activeSize = R;
//	for (UINT i = 0; i < R; i++) {
//		G[i] = -xTz[i];
//		for (UINT k = 0; k < R; k++)
//			G[i] += rpCache[i][k] * H[hInd][k];
//	}
//	// optimization step
//	int iter = 0;
//	while (iter++ < 1000) {
//		for (UINT i = 0; i < activeSize; i++) {
//			UINT j = i + rand() % (activeSize - i);
//			std::swap(index[i], index[j]);
//		}
//		Gmax = -INF;
//		Gmin = INF;
//		for (UINT s = 0; s < activeSize; s++) {
//			UINT vI = index[s]; //relative vector index
//			PG = 0;
//			if (H[hInd][vI] == 0) {
//				if (G[vI] > GmaxOld) {
//					activeSize--;
//					std::swap(index[s], index[activeSize]);
//					s--;
//					continue;
//				} else if (G[vI] < 0)
//					PG = G[vI];
//			} else
//				PG = G[vI];
//
//			if (fabs(PG) > TAU) {
//				double gVal = PG/(REG_PARAM + rpCache[vI][vI]);
//				double lambda2 = std::max(H[hInd][vI] - gVal, LAMBDA_MIN_MF);
//				double delta = lambda2 - H[hInd][vI];
//				if(delta != 0) {
//					dat.putHval(hInd,vI,lambda2);
//					for (UINT k = 0; k < R; k++)
//						G[k] += rpCache[vI][k] * delta;
//				}
//			}
//			Gmax = std::max(Gmax, PG);
//			Gmin = std::min(Gmin, PG);
//		}
//		if (Gmax - Gmin <= GET_H_EPS) {
//			if (activeSize == R)
//				break;
//			else {
//				activeSize = R;
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
//	double normSum = X[hInd].nrm;
//	for(UINT i = 0; i < R; i++) {
//		if(H[hInd][i] != 0) {
//			normSum += (G[i] - xTz[i])*H[hInd][i];
//		}
//	}
//	return normSum;
//}
