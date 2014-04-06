all: AENMF
CXX = g++
CFLAGS = -Wall -Wconversion -O3 -fPIC

inp_params.o: inp_params.cpp
	$(CXX) $(CFLAGS)  -c -o inp_params.o inp_params.cpp 

dense_mat.o: dense_mat.cpp
	$(CXX) $(CFLAGS)  -c -o dense_mat.o dense_mat.cpp 

sparse_mat.o: sparse_mat.cpp
	$(CXX) $(CFLAGS)  -c -o sparse_mat.o sparse_mat.cpp 

all_data.o: all_data.cpp
	$(CXX) $(CFLAGS)  -c -o all_data.o all_data.cpp 

file_int.o: file_int.cpp
	$(CXX) $(CFLAGS)  -c -o file_int.o file_int.cpp 
	
derive_ae.o: derive_ae.cpp
	$(CXX) $(CFLAGS)  -c -o derive_ae.o derive_ae.cpp 

get_factors.o: get_factors.cpp
	$(CXX) $(CFLAGS)  -c -o get_factors.o get_factors.cpp 

AENMF: AENMF.cpp inp_params.o sparse_mat.o dense_mat.o all_data.o file_int.o derive_ae.o get_factors.o
	$(CXX) $(CFLAGS)  -o AENMF AENMF.cpp inp_params.o sparse_mat.o dense_mat.o all_data.o file_int.o derive_ae.o get_factors.o

clean:
	rm -f *~ *.o *.dat AENMF
