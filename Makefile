#Makfile for linux

CXX = g++
CXXFLAGS = -g -O2
CPPFLAGS ?= 
LDLIBS = -lm -lz -lomp

smr: bfile.o CommFunc.o dcdflib.o SMR.o SMR_data.o SMR_data_p1.o SMR_data_p2.o SMR_data_p3.o SMR_plot.o StatFunc.o StrFunc.o
	$(CXX) $(CXXFLAGS)  $?  $(LDLIBS) -o $@

bfile.o: bfile.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

CommFunc.o: CommFunc.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

dcdflib.o: dcdflib.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR.o: SMR.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_data.o: SMR_data.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_data_p1.o: SMR_data_p1.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_data_p2.o: SMR_data_p2.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_data_p3.o: SMR_data_p3.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

SMR_plot.o: SMR_plot.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

StatFunc.o: StatFunc.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?

StrFunc.o: StrFunc.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $?


clean:
	@rm *.o smr


install:
	@echo The binary is under building directory named smr.
