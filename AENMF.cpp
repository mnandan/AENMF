/*
 * onePAE.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: mn
 */

#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <ctime>
#include <algorithm>
#include "fileInt.h"
#include "getFactors.h"
#include "common.h"
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {
	char inpFname[1024];
	char wFname[1024];
	char hFname[1024];

// get user input parameters
	if(argc < 7) {
		cerr<<"Insufficient number of arguments. ";
		cerr<<"AENMF N R D inputFile wFile hFile\n";
		exit(1);
	}
	UINT N = (UINT)strtol(argv[1], NULL, 10);
	UINT R = (UINT)strtol(argv[2], NULL, 10);
	UINT D = (UINT)strtol(argv[3], NULL, 10);
	strcpy(inpFname, argv[4]);
	strcpy(wFname, argv[5]);
	strcpy(hFname, argv[6]);

// initialize data structures
	clock_t begin;
	TrainDat trDat(D, N, R, hFname, wFname);
	fHand_T inpF(inpFname, &trDat);
// read input data
	if(inpF.readFile())
		cout<< "Input file read successfully\n";
	else
		cout<< "Input file has only "<< inpF.getLineNum()<<" lines. Expected "<< N <<" lines\n";
	begin = clock();
	GetFact fact(&trDat);
// Compute initial W and H matrices
	fact.getWinit();
	double frobNorm = fact.getH();
	double prevNorm = frobNorm;
	cout << "init err = "<<frobNorm << ", ";
	for(UINT i = 1; i <= 100; i++) {
		fact.getW();
		frobNorm = fact.getH();
//		cout << frobNorm << " iter "<<i<<endl;
		if(prevNorm - frobNorm < 0.01)
			break;
		prevNorm = frobNorm;
	}
// Compute time taken and write output files
	double cpuTimeTaken = difftime(clock(), begin) / CLOCKS_PER_SEC;
	cout << "final err = " <<frobNorm<<", time = " << cpuTimeTaken <<endl;
	trDat.writeWH();
	return 1;
}

