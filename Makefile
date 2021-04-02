# -----------------------------------------------------------------
#   Supported platforms: Unix / Linux
# ---------------------------------------------------------------------

CXX = g++
CXXFLAGS = -Wall -g -O2

smr_linux: bfile.o CommFunc.o dcdflib.o SMR.o SMR_data.o SMR_data_p1.o SMR_data_p2.o SMR_data_p3.o SMR_plot.o StatFunc.o StrFunc.o
	${CXX} ${CXXFLAGS} bfile.o CommFunc.o SMR.o SMR_data.o SMR_data_p1.o SMR_data_p2.o SMR_data_p3.o SMR_plot.o StatFunc.o StrFunc.o dcdflib.o -lm -lz -lomp -o smr

bfile.o: bfile.cpp bfile.hpp
	${CXX} ${CXXFLAGS} -c bfile.cpp 

CommFunc.o: CommFunc.cpp CommFunc.h
	${CXX} ${CXXFLAGS} -c CommFunc.cpp

dcdflib.o: dcdflib.cpp dcdflib.h cdflib.h ipmpar.h
	${CXX} ${CXXFLAGS} -c dcdflib.cpp

SMR.o: SMR.cpp SMR.h bfile.hpp StrFunc.h SMR_data_p1.h SMR_data_p2.h SMR_data_p3.h SMR_plot.h 
	${CXX} ${CXXFLAGS} -c SMR.cpp

SMR_data.o: SMR_data.cpp SMR_data.h
	${CXX} ${CXXFLAGS} -c SMR_data.cpp

SMR_data_p1.o: SMR_data_p1.cpp SMR_data_p1.h
	${CXX} ${CXXFLAGS} -c SMR_data_p1.cpp

SMR_data_p2.o: SMR_data_p2.cpp SMR_data_p2.h
	${CXX} ${CXXFLAGS} -c SMR_data_p2.cpp

SMR_data_p3.o: SMR_data_p3.cpp SMR_data_p3.h
	${CXX} ${CXXFLAGS} -c SMR_data_p3.cpp

SMR_plot.o: SMR_plot.cpp SMR_plot.h
	${CXX} ${CXXFLAGS} -c SMR_plot.cpp

StatFunc.o: StatFunc.cpp StatFunc.h
	${CXX} ${CXXFLAGS} -c StatFunc.cpp

StrFunc.o: StrFunc.cpp StrFunc.h
	${CXX} ${CXXFLAGS} -c StrFunc.cpp


clean:
	rm *.o
