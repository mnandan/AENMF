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
#include "getFactors.h"
#include <fstream>
#include <exception>

using namespace std;

//void GetFact::getWinit() {
//	UINT numSubs = (UINT)ceil((double)numRp/(double)M_VAL);
//	UINT subSize = (UINT)ceil((double)totVects/(double)numSubs);
//	X2 = new DataVect[totVects];
//	distVals = new DistDat[totVects];
//	clustEndInd = new UINT[totVects];
//
//	UINT cSInd = 0, totBatchRp = 0, cEInd = 0;
//	UINT flsSubSize = min(subSize*50,totVects);
//
//	segregateDataFLS2(0, totVects, flsSubSize);
//	segregateDataSLS(0, totVects, subSize);
//
//	while (cSInd < totVects) {
//		cEInd = clustEndInd[cSInd];
//		UINT clustSize = cEInd - cSInd;
//		if (clustSize > M_VAL) {
//			// find max norm vect: x1
//			UINT maxNormInd = 0;
//			double maxDistVal = -INF;
//			UINT rpInd = 0;
//			for (UINT vI = cSInd; vI < cEInd; vI++) {
//				if (X[vI].nrm > maxDistVal) {
//					maxNormInd = vI;
//					maxDistVal = X[vI].nrm;
//				}
//			}
//			trDat->swapX(cSInd +rpInd, maxNormInd);
//			rpInd++;
//			// find max norm from x1: x2
//			UINT maxNormInd2 = 0;
//			double maxDistVal2 = -INF;
//			for (UINT ind = cSInd + rpInd; ind < cEInd; ind++) {
//				double curr_dist = X[ind].nrm + X[cSInd].nrm - 2 * getDotProductX(cSInd, ind);
//				if (curr_dist > maxDistVal2) {
//					maxNormInd2 = ind;
//					maxDistVal2 = curr_dist ;
//				}
//			}
//			trDat->swapX(cSInd + rpInd, maxNormInd2);
//			rpCache[0][0] = X[cSInd].nrm;
//			rpCache[1][1] = X[cSInd + 1].nrm;
//			rpCache[0][1] = (rpCache[0][0] + rpCache[1][1] - maxDistVal2 ) / 2.0;
//			for (rpInd = 2; rpInd < M_VAL; rpInd++) {
//				maxNormInd = 0;
//				maxDistVal = -INF;
//				for (UINT ind = cSInd + rpInd; ind < cEInd; ind++) {
//					for (UINT rpInd2 = 0; rpInd2 < rpInd; rpInd2++)
//						xTz[rpInd2] = getDotProductX(ind, rpInd2 + cSInd);
//					double rep_err = getRepErr(rpInd, X[ind].nrm);
//
//					if (rep_err > maxDistVal) {
//						maxNormInd = ind;
//						maxDistVal = rep_err;
//					}
//				}
//				trDat->swapX(cSInd + rpInd, maxNormInd);
//				updateCacheX(cSInd, rpInd);
//			}
//		}
//		UINT rpSz = min(clustSize, (UINT) M_VAL);
//		for (UINT ind = 0; ind < rpSz; ind++) {
//			UINT origInd = ind + cSInd;
//			trDat->swapX(totBatchRp, origInd);
//			totBatchRp++;
//		}
//		cSInd = cEInd;
//	}
//
//	delete[] X2;
//	delete[] distVals;
//	delete[] clustEndInd;
//	trDat->copyXW();
//}

void GetFact::getWinit() {
	double *lambda = new double[R];
	UINT *origInd = new UINT [N];
	UINT maxNormInd = 0;
	double maxDistVal = -INF;
	UINT rpInd = 0;
	// find max norm vect: x1
	for (UINT vI = 0; vI < N; vI++) {
		origInd[vI] = vI;
		if (X[vI].nrm > maxDistVal) {
			maxNormInd = vI;
			maxDistVal = X[vI].nrm;
		}
	}
	trDat->swapX(rpInd, maxNormInd);
	swap(origInd[0], origInd[maxNormInd]);
	rpInd++;
	// find max norm from x1: x2
	UINT maxNormInd2 = 0;
	double maxDistVal2 = -INF;
	for (UINT ind = rpInd; ind < N; ind++) {
		double curr_dist = X[ind].nrm + maxDistVal - 2 * getDotProductX(0, ind);
		if (curr_dist > maxDistVal2) {
			maxNormInd2 = ind;
			maxDistVal2 = curr_dist ;
		}
	}
	trDat->swapX(rpInd, maxNormInd2);
	rpCache[0][0] = X[0].nrm;
	rpCache[1][1] = X[1].nrm;
	rpCache[0][1] = (rpCache[0][0] + rpCache[1][1] - maxDistVal2 ) / 2.0;
	swap(origInd[1], origInd[maxNormInd2]);
	// find remaining R - rpInd vectors based on distance from polygon formed by
	// the selected rpInd vectors
	// using H matrix to store dot products temporarily
	for (UINT ind = 2; ind < N; ind++) {
		for (UINT ind2 = 0; ind2 < 2; ind2++)
			trDat->putHval(origInd[ind],ind2,getDotProductX(ind2, ind));
	}
	for (rpInd = 2; rpInd < min(N,R); rpInd++) {
		maxNormInd = 0;
		maxDistVal = -INF;
		for (UINT ind = rpInd; ind < N; ind++) {
			double rep_err = getRepErr(rpInd, X[ind].nrm, lambda, origInd[ind]);
			if (rep_err > maxDistVal) {
				maxNormInd = ind;
				maxDistVal = rep_err;
			}
		}
		trDat->swapX(rpInd, maxNormInd);
		swap(origInd[rpInd],origInd[maxNormInd]);
		if(rpInd < min(N,R) - 1) {
			updateCacheX(0, rpInd);
			for (UINT ind = rpInd + 1; ind < N; ind++)
				trDat->putHval(origInd[ind],rpInd,getDotProductX(rpInd, ind));
		}
	}
	trDat->copyXW();
	trDat->clearH();
	delete [] origInd;
	delete [] lambda;
}


double GetFact::getH() {
	return deriveH();
}

double GetFact::getW() {
	return deriveW();
}
