#Makfile for linux

CXX = g++
CXXFLAGS = -g -O2
CPPFLAGS ?= 
LDLIBS = -lm -lz -lomp

smr: CPP/bfile.o CPP/CommFunc.o CPP/dcdflib.o CPP/SMR.o CPP/SMR_data.o \
	 CPP/SMR_data_p1.o CPP/SMR_data_p2.o CPP/SMR_data_p3.o CPP/SMR_plot.o \
	 CPP/StatFunc.o CPP/StrFunc.o CPP/calmt.o CPP/main.o
	$(CXX) $(CXXFLAGS)  $?  $(LDLIBS) -o $@

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

calmt.o: CPP/calmt.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -c $?

main.o: CPP/main.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -c $?

clean:
	@rm C/*.o
	@rm CPP/*.o
	@rm smr


install:
	@echo The binary is under building directory named smr.



