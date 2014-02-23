all: DRSpl
CXX = g++
CFLAGS = -Wall -Wconversion -O3 -fPIC

fileInt.o: fileInt.cpp
	$(CXX) $(CFLAGS)  -c -o fileInt.o fileInt.cpp 

deriveAE.o: deriveAE.cpp
	$(CXX) $(CFLAGS)  -c -o deriveAE.o deriveAE.cpp 

onePassDRS.o: onePassDRS.cpp
	$(CXX) $(CFLAGS)  -c -o onePassDRS.o onePassDRS.cpp 

aesvmSolver.o: aesvmSolver.cpp
	$(CXX) $(CFLAGS)  -c -o aesvmSolver.o aesvmSolver.cpp 

DRSpl: onePAE.cpp deriveAE.o onePassDRS.o aesvmSolver.o fileInt.o
	$(CXX) $(CFLAGS)  -o DRSpl onePAE.cpp deriveAE.o onePassDRS.o aesvmSolver.o fileInt.o

clean:
	rm -f *~ *.o DRSpl
