/*
 * common.h
 *
 *  Created on: Aug 3, 2013
 *      Author: mnandan
 */

#ifndef FILE_INT_H_
#define FILE_INT_H_

#include <fstream>
#include <vector>
#include <string>

#include "common.h"
#include "inp_params.h"
#include "all_data.h"

class FileInt {
// file handles
	std::ifstream xIN;
	std::ofstream wOUT, hOUT;
// references to data
	AllData &dat;
// data attributes
	UINT N, maxD, maxNNZ;
// variables
	bool *delims1, *delims2;
	std::vector <FeatType> Ftemp;
public:
	FileInt(InpParams const &pars, AllData &dat_);
	void init(InpParams const &pars);
	~FileInt() {
		xIN.close();
		wOUT.close();
		hOUT.close();
		delete [] delims1;
		delete [] delims2;
	}
	UINT getMaxD() {
		return maxD;
	}
	bool procLine();
	void readXFile();
	void writeW();
	void writeH();
};
#endif /* FILE_INT_H_ */
