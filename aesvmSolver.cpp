/*
 * aesvmSolver.cpp
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "aesvmSolver.h"
#include <ctime>

using namespace std;

void AesvmSolver_T::solveAesvm() {
	UINT epochSize = trDat->getEphSz();
	const dataVect_T *X = trDat->XC;
	const double* w = trDat->wC;
	const double C = trDat->C;
	srand((UINT)time(NULL));

	// PG: projected gradient, for shrinking and stopping
	double PG;
	PGmax_old = INF;
	PGmin_old = -INF;
	double PGmax_new, PGmin_new;
	double alpha_new;

	index = new UINT[epochSize];
	for (UINT vI = 0; vI < epochSize; vI++) {
		index[vI] = vI;
	}
	int maxIterNum = 1000;

	int iterNum = 0;
	UINT numVects = epochSize;//trDat->getXsize() - trDat->getRpNum();
	unsigned active_size = epochSize;

	while (iterNum < maxIterNum) {
		UINT i, j;
			PGmax_new = -INF;
			PGmin_new = INF;

			for (i = 0; i < active_size; i++) {
				j = i + rand() % (active_size - i);
				swap(index[i], index[j]);
			}

			for (UINT s = 0; s < active_size; s++) {
				UINT vI = index[s]; // vector index
				double Gradient = 0;
				double Cnow = C;
				const feat_T *F = trDat->getF(vI);
				int y = (int) X[vI].label;
				for (UINT fI = 0; fI < X[vI].numFeats; fI++) //feature index
					Gradient += w[F[fI].fNum] * F[fI].fVal;
				Gradient = Gradient * y - 1;

				PG = 0;
				if (X[vI].A == 0) {
					if (Gradient > PGmax_old) {
						active_size--;
						swap(index[s], index[active_size]);
						s--;
						continue;
					} else if (Gradient < 0)
						PG = Gradient;
				} else if (X[vI].A == Cnow) {
					if (Gradient < PGmin_old) {
						active_size--;
						swap(index[s], index[active_size]);
						s--;
						continue;
					} else if (Gradient > 0)
						PG = Gradient;
				} else
					PG = Gradient;

				PGmax_new = max(PGmax_new, PG);
				PGmin_new = min(PGmin_new, PG);

				if (fabs(PG) > TAU) {
					alpha_new = min(max(X[vI].A - Gradient / X[vI].nrm, 0.0), Cnow); // box constraints
					trDat->putAW(vI, alpha_new);
				}
				//unlock vi
			}

			iterNum++;
			//if (iterNum % 10 == 0)
				//cout << ".";

			if (PGmax_new - PGmin_new <= 0.1) {
				if (active_size == numVects)
					break;
				else {
					active_size = numVects;
					//cout << "*";
					PGmax_old = INF;
					PGmin_old = -INF;
					continue;
				}
			}
			PGmax_old = PGmax_new;
			PGmin_old = PGmin_new;
			if (PGmax_old <= 0)
				PGmax_old = INF;
			if (PGmin_old >= 0)
				PGmin_old = -INF;
		}

	delete[] index;
}

void AesvmSolver_T::solveAesvm2() {
	UINT epochSize = trDat->getEphSz();
	const dataVect_T *X = trDat->XC;
	const double* w = trDat->wC;
	const double C = trDat->C;
	srand((UINT)time(NULL));

	// PG: projected gradient, for shrinking and stopping
	double PG;
	PGmax_old = INF;
	PGmin_old = -INF;
	double PGmax_new, PGmin_new;
	double alpha_new;

	index = new UINT[epochSize];
	for (UINT vI = 0; vI < epochSize; vI++) {
		index[vI] = vI;
		//trDat->putA(vI,0);
	}
	int maxIterNum = 1000;

	int iterNum = 0;
	UINT numVects = epochSize;
	unsigned active_size = epochSize;

	while (iterNum < maxIterNum) {
		UINT i, j;
			PGmax_new = -INF;
			PGmin_new = INF;

			for (i = 0; i < active_size; i++) {
				j = i + rand() % (active_size - i);
				swap(index[i], index[j]);
			}

			for (UINT s = 0; s < active_size; s++) {
				UINT vI = index[s]; // vector index
				double Gradient = 0;
				double Cnow = C * X[vI].B;
				const feat_T *F = trDat->getF(vI);
				int y = (int) X[vI].label;
				for (UINT fI = 0; fI < X[vI].numFeats; fI++) //feature index
					Gradient += w[F[fI].fNum] * F[fI].fVal;
				Gradient = Gradient * y - 1;

				PG = 0;
				if (X[vI].A == 0) {
					if (Gradient > PGmax_old) {
						active_size--;
						swap(index[s], index[active_size]);
						s--;
						continue;
					} else if (Gradient < 0)
						PG = Gradient;
				} else if (X[vI].A >= C) {
					if (Gradient < PGmin_old* X[vI].A/C) {
						active_size--;
						swap(index[s], index[active_size]);
						s--;
						continue;
					} else if (Gradient > 0)
						PG = Gradient*C/X[vI].A;
				} else
					PG = Gradient;

				PGmax_new = max(PGmax_new, PG);
				PGmin_new = min(PGmin_new, PG);

				if (fabs(PG) > TAU) {
					alpha_new = min(max(X[vI].A - Gradient / X[vI].nrm, 0.0), Cnow); // box constraints
					trDat->putAW(vI, alpha_new);
				}				
			}

			iterNum++;
			//if (iterNum % 50 == 0)
				//cout << ".";

			if (PGmax_new - PGmin_new <= 0.1) {
				if (active_size == numVects)
					break;
				else {
					active_size = numVects;
					//cout << "*";
					PGmax_old = INF;
					PGmin_old = -INF;
					continue;
				}
			}
			PGmax_old = PGmax_new;
			PGmin_old = PGmin_new;
			/*if (PGmax_old <= 0)
				PGmax_old = INF;
			if (PGmin_old >= 0)
				PGmin_old = -INF;*/
		}

	delete[] index;
}

