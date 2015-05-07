
#GSL_L	=	/usr/local/lib/
#GSL_I	=	/usr/local/include/


GSL_L	=	/usr/lib/
GSL_I	=	/usr/include/



clean: rm *.o

shell_iso_hw: work.o HW_iso.o
	g++  -I$(GSL_I) -L$(GSL_L) -fopenmp -lgsl -lgslcblas -lm  work.o HW_iso.o -o shell_iso_hw

work.o: work.cpp HW_iso.h
	g++ -c work.cpp

HW_iso.o: HW_iso.cpp HW_iso.h
	g++ -c HW_iso.cpp 
