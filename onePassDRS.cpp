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
#include "onePassDRS.h"
#include <fstream>

using namespace std;

double onePDRS_T::onePassDRS(UINT numRp, bool isFirstRun) {
	UINT numVects = trDat->getXsize();
	return onePassDRSsub(numRp, 0, numVects);
}

double onePDRS_T::onePassDRSsub(UINT numRp, UINT startInd, UINT endInd) {
	UINT numVects = endInd - startInd;
	UINT numSubs = ceil((double)numRp/(double)M_VAL);
	UINT subSize = ceil((double)numVects/(double)numSubs);
	X2 = new dataVect_T[numVects];
	distVals = new distDatType[numVects];
	clustEndInd = new UINT[endInd];
	lbdArr = new double*[subSize]; // size VxM
	rpCache = new double*[numRp]; // size MxM
	for (UINT ind = 0; ind < numRp; ind++) {
		lbdArr[ind] = new double[numRp];
		rpCache[ind] = new double[numRp];
	}
	for (UINT ind = numRp; ind < subSize; ind++)
		lbdArr[ind] = new double[numRp];
	xTz = new double[numRp];

	UINT prevCSI = 0;
	UINT prevCS = 0;
	UINT cSInd = startInd; //clustStartIndex
	UINT cEInd; //clustEndIndex
	UINT totBatchRp = 0;

	UINT flsSubSize = min(subSize*50,numVects);
	segregateDataFLS2(startInd, endInd, flsSubSize);
	segregateDataSLS(startInd, endInd, subSize);

	while (cSInd < endInd) {
		cEInd = clustEndInd[cSInd];
		UINT clustSize = cEInd - cSInd;
		for (UINT ind = 0; ind < M_VAL; ind++) {
			for (UINT rpInd = 0; rpInd < M_VAL; rpInd++)
				lbdArr[ind][rpInd] = 0;
			lbdArr[ind][ind] = 1;
		}

		if (clustSize > M_VAL) {
			UINT rpInd = 0;
			// find max norm vect: x1
			UINT maxNormInd = 0;
			double maxDistVal = -INF;
			for (UINT vI = cSInd; vI < cEInd; vI++) {
				if (X[vI].nrm > maxDistVal) {
					maxNormInd = vI;
					maxDistVal = X[vI].nrm;
				}
			}
			trDat->swapX(cSInd + rpInd, maxNormInd);
			rpInd++;

			// find max norm from x1: x2
			UINT maxNormInd2 = 0;
			double maxDistVal2 = -INF;
			for (UINT ind = cSInd + rpInd; ind < cEInd; ind++) {
				double curr_dist = X[ind].nrm + X[cSInd].nrm - 2 * getDotProduct(cSInd, ind);
				if (curr_dist > maxDistVal2) {
					maxNormInd2 = ind;
					maxDistVal2 = curr_dist ;
				}
			}
			trDat->swapX(cSInd + rpInd, maxNormInd2);
			rpCache[0][0] = X[cSInd].nrm;
			rpCache[1][1] = X[cSInd + 1].nrm;
			rpCache[0][1] = (rpCache[0][0] + rpCache[1][1] - maxDistVal2 ) / 2.0;
			rpInd++;
			// find x_i that gives max f value: x3 ... xM
			for (rpInd = 2; rpInd < M_VAL; rpInd++) {
				maxNormInd = 0;
				maxDistVal = -INF;
				for (UINT ind = cSInd + rpInd; ind < cEInd; ind++) {
					//currEndIndex <- startIndex + rpInd
					for (UINT rpInd2 = 0; rpInd2 < rpInd; rpInd2++)
						xTz[rpInd2] = getDotProduct(ind, rpInd2 + cSInd);
					double rep_err = getRepErr(rpInd, X[ind].nrm, lbdArr[subSize - 1]);

					if (rep_err > maxDistVal) {
						maxNormInd = ind;
						maxDistVal = rep_err;
					}
				}
				trDat->swapX(cSInd + rpInd, maxNormInd);
				updateCache(cSInd, rpInd);
				//cout<< rpInd<<endl;
			}
		}

		UINT rpSz = min(clustSize, (UINT) M_VAL);

		for (UINT ind = 0; ind < rpSz; ind++) {
			UINT origInd = ind + cSInd;
			trDat->swapRp(totBatchRp, origInd);
			trDat->writeAE(totBatchRp);
			totBatchRp++;
		}

		prevCSI = cSInd;
		prevCS = min(clustSize, (UINT) M_VAL);
		cSInd = cEInd;
	}
    for (UINT ind = 0; ind < numRp; ind++) {
        for (UINT rpInd = 0; rpInd < numRp; rpInd++)
            lbdArr[ind][rpInd] = 0;
        lbdArr[ind][ind] = 1;
    }
    //cout<<totBatchRp<<endl;
    for (UINT ind = 0; ind < numRp; ind++) {
        updateCache(0, ind);
        trDat->writeLBD(numRp, lbdArr[ind]);
        trDat->writeX(ind);
    }
// derive lambda for other vectors
    for (UINT ind = numRp; ind < numVects; ind++) {
        for (UINT rpInd2 = 0; rpInd2 < numRp; rpInd2++)
            xTz[rpInd2] = getDotProduct(ind, rpInd2);
        deriveAE(numRp, lbdArr[numRp]);
        trDat->writeLBD(numRp, lbdArr[numRp]);
        trDat->writeX(ind);        
    }

	for (UINT ind = 0; ind < numRp; ind++) {
		delete[] rpCache[ind];
		delete[] lbdArr[ind];
	}
	for (UINT ind = numRp; ind < subSize; ind++)
		delete[] lbdArr[ind];
	delete[] rpCache;
	delete[] lbdArr;
	delete[] xTz;
	delete[] X2;
	delete[] distVals;
	delete[] clustEndInd;

	return 0;
}

// C code for the quickselect algorithm (Hoare's selection algorithm), Ref. Numerical Recipes in C
// selects element in array (arr) that is kth in largest among the first n elements
// quickselect has an expected complexity of O(N) and worst case complexity of O(N^2), though in practice it is fast

distDatType onePDRS_T::quickSelectDist(UINT n, UINT k) {
	UINT beforeK, afterK, mid;
	UINT i, j;
	distDatType a; // a is the partition element

	beforeK = 0;
	afterK = n - 1;
	while (1) {
		if (afterK <= beforeK + 1) { //Active partition contains 1 or 2 elements
			if (afterK == beforeK + 1 && distVals[afterK].dist < distVals[beforeK].dist) //Case of 2 elements
				swap(distVals[beforeK], distVals[afterK]);
			return distVals[k];
		} else {
			mid = (beforeK + afterK) / 2;
			swap(distVals[mid], distVals[beforeK + 1]);
// rearrange distVals so that distVals[afterK] >= distVals[beforeK + 1] >= distVals[beforeK]
			if (distVals[beforeK].dist > distVals[afterK].dist)
				swap(distVals[beforeK], distVals[afterK]);
			if (distVals[beforeK + 1].dist > distVals[afterK].dist)
				swap(distVals[beforeK + 1], distVals[afterK]);
			if (distVals[beforeK].dist > distVals[beforeK + 1].dist)
				swap(distVals[beforeK], distVals[beforeK + 1]);
// distVals[beforeK + 1] is selected as the partition element, a, with sentinels, distVals[afterK] and distVals[beforeK]
			a = distVals[beforeK + 1];
// find all {(i,j): j>= i and distVals[i] < a < distVals[j]} and swap
// at the end of the while loop, we have distVals[beforeK]....distVals[i] <= a <= distVals[j]....distVals[afterK]
			i = beforeK + 1;
			j = afterK;
			while (1) {
				do
					i++;
				while (distVals[i].dist < a.dist);
				do
					j--;
				while (distVals[j].dist > a.dist);
				if (j < i)
					break;
				swap(distVals[i], distVals[j]);
			}
// put partition element, a, in correct position in rearranged array
			distVals[beforeK + 1] = distVals[j];
			distVals[j] = a;
// if the position of a, j, is ahead of k, search only till j - 1 in the next iteration
			if (j >= k)
				afterK = j - 1;
// if the position of a, j, is before k, search only from i
			if (j <= k)
				beforeK = i;
		}
	}
}

// Implements FLS2 of DeriveRS. Data is divided into sets of size less than P
void onePDRS_T::segregateDataFLS2(UINT startIndex, UINT endIndex, UINT P_VAL) {
	UINT anchorInd = -1;
	double maxDist = -INF;
	for (UINT k = startIndex; k < endIndex; k++) {	//find max as anchor
		if(maxDist < X[k].nrm) {
			maxDist = X[k].nrm;
			anchorInd = k;
		}
	}

	if (endIndex - startIndex >= 2 * P_VAL) { //steps 1,2,3, and 5
		UINT p = 0;
		UINT k;
		UINT midVal;
		for (k = startIndex; k < endIndex; k++) { // step 1
			distVals[p].dist = getDotProduct(k, anchorInd);
			distVals[p++].ind = k;
		}
		midVal = p / 2;
		quickSelectDist(p, midVal); // step 2
// steps 3 and 5
		trDat->putX(startIndex, p, X2, distVals);

		midVal += startIndex;
		segregateDataFLS2(startIndex, midVal, P_VAL);
		segregateDataFLS2(midVal, endIndex, P_VAL);
	} else if (endIndex - startIndex > P_VAL) { //steps 1,2,3, and 4
		UINT p = 0;
		UINT k;
		UINT midVal;
		for (k = startIndex; k < endIndex; k++) {
			distVals[p].dist = getDotProduct(k, anchorInd);
			distVals[p++].ind = k;
		}
		midVal = p / 2;
		quickSelectDist(p, midVal);

		trDat->putX(startIndex, p, X2, distVals);

		for (k = 0; k < midVal; k++)
			clustEndInd[k + startIndex] = midVal + startIndex;
		midVal += startIndex;
		for (k = midVal; k < endIndex; k++) //step 4
			clustEndInd[k] = endIndex;
	} else { //step 4
		for (UINT p = startIndex; p < endIndex; p++)
			clustEndInd[p] = endIndex;
	}
}

/// Implements SLS of DeriveRS. Data is divided into sets of size less than V
void onePDRS_T::segregateDataSLS(UINT startIndex, UINT endIndex, UINT V_VAL) {
	UINT pSubsetStart = startIndex;
	UINT pSubsetEnd = clustEndInd[pSubsetStart];
	while(1) {
		UINT vSubsetStart = pSubsetStart;
		UINT vSubsetEnd = vSubsetStart + V_VAL;
		while(1) {
			if(vSubsetEnd >= pSubsetEnd) {
				for (UINT k = vSubsetStart; k < pSubsetEnd; k++)
					clustEndInd[k] = pSubsetEnd;
				break;
			}
			UINT k, p = 0;
			UINT anchorInd = -1;
			double maxDist = -INF;
			for (k = vSubsetStart; k < pSubsetEnd; k++) {	//find max as anchor
				if(maxDist < X[k].nrm) {
					maxDist = X[k].nrm;
					anchorInd = k;
				}
			}
			for (k = vSubsetStart; k < pSubsetEnd; k++) {
				distVals[p].dist = X[k].nrm - 2*getDotProduct(k, anchorInd);
				distVals[p++].ind = k;
			}
			quickSelectDist(p, V_VAL); // select closest V vectors
			trDat->putX(vSubsetStart, p, X2, distVals);
			for (k = vSubsetStart; k < pSubsetEnd; k++)
				clustEndInd[k] = vSubsetEnd;
			vSubsetStart = vSubsetEnd;
			vSubsetEnd += V_VAL;
		}
		if(pSubsetEnd < endIndex) {		// step 2
			pSubsetStart = pSubsetEnd;
			pSubsetEnd = clustEndInd[pSubsetStart];
		}
		else
			break;
	}
}
