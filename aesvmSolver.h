/*
 * aesvmSolver.h
 *
 *  Created on: Aug 6, 2013
 *      Author: mnandan
 */

#ifndef AESVMSOLVER_H_
#define AESVMSOLVER_H_

#include "deriveAE.h"
#include <vector>

class AesvmSolver_T {
private:
	UINT *index;
	UINT totVects;

protected:
	trainDat_T *trDat;
	double PGmax_old;
	double PGmin_old;

public:
	AesvmSolver_T(trainDat_T *trDat, UINT totVects) {
		this->trDat = trDat;
		this->totVects = totVects;
		PGmax_old = INF;
		PGmin_old = -INF;
	}

	void solveAesvm();
	void solveAesvm2();
};

#endif /* AESVMSOLVER_H_ */
