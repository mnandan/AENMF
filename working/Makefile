all: AENMF
CXX = g++
CFLAGS = -Wall -Wconversion -O3 -fPIC

fileInt.o: fileInt.cpp
	$(CXX) $(CFLAGS)  -c -o fileInt.o fileInt.cpp 

deriveAE.o: deriveAE.cpp
	$(CXX) $(CFLAGS)  -c -o deriveAE.o deriveAE.cpp 

getFactors.o: getFactors.cpp
	$(CXX) $(CFLAGS)  -c -o getFactors.o getFactors.cpp 

AENMF: AENMF.cpp deriveAE.o getFactors.o fileInt.o
	$(CXX) $(CFLAGS)  -o AENMF AENMF.cpp deriveAE.o getFactors.o fileInt.o

clean:
	rm -f *~ *.o *.dat AENMF
