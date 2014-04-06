/*
 * inp_params.h
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#ifndef INP_PARAMS_H_
#define INP_PARAMS_H_

#include "common.h"
#include <string>

class InpParams {
	UINT R_;
	std::string xF_, wF_, hF_;
	bool vrbs_;
public:
	const UINT &R;
	const std::string &xFile, &wFile, &hFile;
	const bool &verbose;
	// bind const public references to private variables
	InpParams(): R(R_), xFile(xF_), wFile(wF_), hFile(hF_), verbose(vrbs_) {
		R_ = 20;
		xF_ = "";
		wF_ = "wMat.dat";
		hF_ = "hMat.dat";
		vrbs_ = true;
	}
	bool init(int argc, const char* const params[]);
	void dispParList();
};


#endif /* INP_PARAMS_H_ */
