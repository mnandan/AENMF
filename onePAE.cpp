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
#include "onePassDRS.h"
#include "common.h"

using namespace std;



void reader(fHand_T &inpF, onePDRS_T &onePD, UINT aeSize, UINT &allRpNum) {
	int blockNum = 1;
	UINT vectRead = 0;
	UINT totVects = inpF.getTotVects();
	UINT eSz = inpF.getEpochSize();
	while (vectRead < totVects) {
		UINT currES = min(vectRead + eSz, totVects) - vectRead;
		UINT linesRead = 0;
		UINT epochRpSz = aeSize;
		while (linesRead < currES && inpF.procLine()) {
			++linesRead;
		}
		vectRead += linesRead;

		if (linesRead == 0)
			break;

		cout << "vectRead = "<<vectRead<<" ER = ," <<onePD.onePassDRS(epochRpSz, true)<<endl;
		blockNum++;
	}
}

int main(int argc, char *argv[]) {

	char inpFname[1024];
	char aeFname[1024];
	char modelFname[1024];
	UINT aeSize = 0;
	UINT maxF;
	UINT epochSize;

	if(argc < 5) {
		cerr<<"Insufficient number of arguments\n";
		exit(1);
	}

	epochSize = strtol(argv[1], NULL, 10);
	maxF = strtol(argv[2], NULL, 10);
	aeSize = strtol(argv[3], NULL, 10);
	strcpy(inpFname, argv[4]);

	clock_t begin = clock();
	trainDat_T TrD(maxF, epochSize, aeFname, modelFname, 1);
	fHand_T inpF(inpFname, epochSize, maxF, &TrD);
	onePDRS_T onePD(&TrD, epochSize, aeSize);
	int retCode = 0;

	UINT allRpNum = 0;
	reader(inpF, onePD, aeSize, allRpNum);
	double cpuTimeTaken = difftime(clock(), begin) / CLOCKS_PER_SEC;
	cout << "CPU time taken = ,0," << cpuTimeTaken <<endl;
	return retCode;
}

