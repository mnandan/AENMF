/*
 * onePAE.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: mn
 */

#include <iostream>
#include <ctime>
#include <cmath>
#include "common.h"
#include "file_int.h"
#include "get_factors.h"
#include "all_data.h"
#include "derive_ae.h"

int main(int argc, char *argv[]) {
// get user input parameters and initialize DS
	InpParams pars;
	if(pars.init(argc, argv)) {
		pars.dispParList();
		return 1;
	}
	// Initialize data matrices X, W and H in dat
    AllData dat(pars.R);
	FileInt files(pars, dat);    // to read and write files
	try {	// read input data from pars.xFile
		files.readXFile();
	} catch (std::string &msg) {
		std::cerr<<msg<<std::endl;
		return 1;
	}
	if(pars.verbose)
		std::cout<< "Input file read successfully"<<std::endl;
	dat.initWH(files.getMaxD());
	// DeriveAEW performs initialization of W
	DeriveAEW winit(dat);
	// GetFact performs factorization of dat.X to W and H
	GetFact fact(dat);

// Compute initial W and H matrices
	clock_t begin = clock();
	winit.getW();
	double frobNorm = fact.getH();
	if(pars.verbose)
		std::cout << "Initial reconstruction error = "<<frobNorm<<std::endl;
	else
		std::cout << "Computing factors ";
// Perform ALS till convergence criteria is met
	double prevNorm;
	UINT iter = 1;
	while(iter <= 100) {
		prevNorm = frobNorm;
		fact.getW();    // optimize over W
		frobNorm = fact.getH();    // optimize over H
		if(pars.verbose)
			std::cout << "\titeration "<< iter <<": reconstruction error = " \
			          <<frobNorm<<std::endl;
		else if(iter%5 == 0)
			std::cout<<'.';
		iter++;
		// check if convergence criteria is met
		if(fabs(prevNorm - frobNorm) < 0.01)
			break;
	}
	if(!pars.verbose && iter > 5)
		std::cout << std::endl;    // for pretty display after '.'
// Compute time taken and write output files
	double cpuTimeTaken = difftime(clock(), begin) / CLOCKS_PER_SEC;
	std::cout << "Final reconstruction error = " <<frobNorm<<", time = " \
			  << cpuTimeTaken <<std::endl;
	files.writeW();
	files.writeH();
	return 0;
}

