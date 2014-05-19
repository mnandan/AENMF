/*
 * data.cpp
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#include "all_data.h"

void AllData::copyXW() {	//initialize W matrix with the first R vectors in X
	for(UINT i = 0; i< R; i++) {
		const FeatVect &F = XV[i].F;
		for(UINT j = 0; j < XV[i].fSize; j++)
			putWval(i, F[j].fNum, F[j].fVal);
	}
}
void AllData::putHvect(UINT i, double * lambda) {
	for (UINT k = 0; k < R; k++)
		H_.put(i,k,lambda[k]);
}
void AllData::putWvect(UINT i, double * lambda) {
	for (UINT k = 0; k < D; k++)
		W_.put(i,k,lambda[k]);
}

void AllData::clearW() {
	for (UINT ind1 = 0; ind1 < R; ind1++)
		for (UINT ind2 = 0; ind2 < D; ind2++)
			W_.put(ind1,ind2,0);
}

void AllData::clearH() {
	for (UINT ind1 = 0; ind1 < N; ind1++)
		for (UINT ind2 = 0; ind2 < R; ind2++)
			H_.put(ind1,ind2,0);
}

double AllData::calcNorm() {
	double frobNorm = 0;
	std::vector <double> hwVect;
	for(UINT j = 0;j < D;j++)
		hwVect.push_back(0);

	for(UINT i = 0;i < N;i++) {
		for(UINT j = 0;j < D;j++)
			hwVect[j] = 0;
		for(UINT j = 0;j < R;j++) {
			if(HV[i][j] != 0 )
				for(UINT k = 0;k < D;k++)
					hwVect[k] += HV[i][j]*WV[j][k];
		}
		UINT fSize = XV[i].fSize;
		UINT prevFnum = 0;
		for(UINT j = 0; j < fSize; j++) {
			for(;prevFnum < XV[i].F[j].fNum; prevFnum++)
				frobNorm += hwVect[prevFnum] * hwVect[prevFnum];
			double val = XV[i].F[j].fVal - hwVect[XV[i].F[j].fNum];
			frobNorm += val*val;
			prevFnum++;
		}
	}
	return frobNorm;
}
