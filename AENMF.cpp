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
	UINT N = strtol(argv[1], NULL, 10);
	UINT R = strtol(argv[2], NULL, 10);
	UINT D = strtol(argv[3], NULL, 10);
	strcpy(inpFname, argv[4]);
	strcpy(wFname, argv[5]);
	strcpy(hFname, argv[6]);

// initialize data structures
	clock_t begin;
	vector < vector< double> > H;	//factor H, W is stored in X of trDat
	for(UINT i = 0; i < N; i++)
		H.push_back(vector <double> ());	//push empty vector
	TrainDat trDat(D, N, wFname, hFname);
	fHand_T inpF(D, N, inpFname, &trDat);
// read input data
	begin = clock();
	if(inpF.readFile())
		cout<< "Input file read successfully\n";
	else
		cout<< "Input file has only "<< inpF.getLineNum()<<" lines. Expected "<< N <<" lines\n";

#if 0
	{
		UINT R1 = 100*R;
		GetFact fact1(&trDat, R1);
	// Compute V a Dx100R matrix
		fact1.getV(R1,0,N);
		cout<< "Computed V successfully\n";
	// Compute W a DxR matrix from V
		fact1.getV(R,0,R1);
		cout<< "Computed W successfully\n";
	}
	GetFact fact2(&trDat, R);
// Compute H a RxN matrix using W
	fact2.getH(R,0,N,H);
	cout<< "Computed H successfully\n";
#else
	GetFact fact(&trDat, R);
// Compute V a Dx100R matrix
	fact.getFactors(R,0,N,H);
#endif
// Compute time taken and write output files
	double cpuTimeTaken = difftime(clock(), begin) / CLOCKS_PER_SEC;
	cout << "CPU time taken = " << cpuTimeTaken <<endl;
	trDat.writeWH(R,N,H);
	return 1;
}

