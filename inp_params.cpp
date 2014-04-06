/*
 * userInp.cpp
 *
 *  Created on: Apr 5, 2014
 *      Author: zarcon
 */

#include "inp_params.h"
#include <string>
#include <cstdlib>
#include <iostream>

bool InpParams::init(int argc, const char* const params[]) {
	switch(argc) {
	case 6:
		if(strtol(params[5], NULL, 10) == 0)
			vrbs_ = false;
	case 5:
		hF_ = params[4];
	case 4:
		wF_ = params[3];
	case 3:
		R_ = strtol(params[2], NULL, 10);
		if(R_ < 1)
			return false;
	case 2:
		xF_ = params[1];
		return false;
	default:
		return true;
	}
}

void InpParams::dispParList() {
	std::cout << "Please call executable as: ";
	std::cout << "./AENMF <input X file> <optional parameters>\n";
	std::cout<<"optional parameter 1: <R default=20>\n";
	std::cout<<"optional parameter 2: <output W file default=wMat.dat>\n";
	std::cout<<"optional parameter 3: <output H file default=hMat.dat>\n";
}


