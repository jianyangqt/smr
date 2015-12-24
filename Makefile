
# -----------------------------------------------------------------
#   Makefile for GCTA 
#   
#   Supported platforms: Unix / Linux
# ---------------------------------------------------------------------

# Directory of the target
OUTPUT = smr_linux

# Compiler
CXX = g++

# EIGEN library
EIGEN_PATH = ../Lib/eigen

# Intel MKL library
#MKL_PATH = /opt/intel/mkl

# Compiler flags
CXXFLAGS = -w -O3 -m64 -fopenmp -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG 
LIB += -static -lz -Wl,-lm -ldl
#LIB += -lz -Wl, -lm -ldl

HDR += CommFunc.h \
	   cdflib.h \
	   dcdflib.h \
           SMR.h \
	   ipmpar.h \
           StatFunc.h \
           StrFunc.h \
            SMR_data.h \
            eData.h 
SRC = SMR.cpp \
           CommFunc.cpp \
           SMR_data.cpp \
	   dcdflib.cpp \
           StatFunc.cpp \
           StrFunc.cpp	\
           eData.cpp 
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean: 
	rm -f *.o
