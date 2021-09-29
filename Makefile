#Makfile for linux
CC = gcc
CFLAGS = -O3
CXX = g++
CXXFLAGS = -O3

CPPFLAGS = 
LDFLAGS = 
LIBS =  -lm -lz -lomp

all: smr smr_static

smr: CPP/bfile.o CPP/CommFunc.o CPP/dcdflib.o CPP/SMR.o CPP/SMR_data.o \
	 CPP/SMR_data_p1.o CPP/SMR_data_p2.o CPP/SMR_data_p3.o CPP/SMR_plot.o \
	 CPP/StatFunc.o CPP/StrFunc.o C/calmt.o C/file_parser.o CPP/main.o
	$(CXX) $(CXXFLAGS)  $?  $(LIBS) -o $@

smr_static: CPP/bfile.o CPP/CommFunc.o CPP/dcdflib.o CPP/SMR.o CPP/SMR_data.o \
	 CPP/SMR_data_p1.o CPP/SMR_data_p2.o CPP/SMR_data_p3.o CPP/SMR_plot.o \
	 CPP/StatFunc.o CPP/StrFunc.o C/calmt.o C/file_parser.o CPP/main.o
	$(CXX) $(CXXFLAGS) $? $(LIBS) -static -ldl -lpthread -o $@

bfile.o: CPP/bfile.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

CommFunc.o: CPP/CommFunc.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

dcdflib.o: CPP/dcdflib.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR.o: CPP/SMR.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_data.o: CPP/SMR_data.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_data_p1.o: CPP/SMR_data_p1.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_data_p2.o: CPP/SMR_data_p2.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_data_p3.o: CPP/SMR_data_p3.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_plot.o: CPP/SMR_plot.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

StatFunc.o: CPP/StatFunc.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

StrFunc.o: CPP/StrFunc.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

calmt.o: C/calmt.c
	${CC} ${CFLAGS} ${CPPFLAGS} -c $?

file_parser: C/file_parser.c
	${CC} ${CFLAGS} ${CPPFLAGS} -c $?

main.o: CPP/main.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -c $?


clean:
	@rm C/*.o
	@rm CPP/*.o
	@rm smr


install:
	@echo The binary is under building directory named smr.



